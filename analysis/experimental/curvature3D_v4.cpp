//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <deque>
#include <list>
#include <chrono>

#define sThreshold 20
#define pSig 0.3

struct edgeType {
	int x,y;
};

//for sorting keys, logarithmic search
bool compByX(edgeType a, edgeType b)
{
	//if(a.x<b.x)
	//	return true;
	//return false;
	return a.x<b.x;
}

//for sorting keys, logarithmic search
bool equByXY(edgeType a, edgeType b)
{
	return (a.x==b.x) && (a.y==b.y);
}


//for simple linear search
class equByY
{
	public:
		bool operator () (edgeType a)
		{
			return a.y==y;
		};
		int x;
		int y;
};

//for simple linear search
class equByX
{
	public:
		bool operator () (edgeType a)
		{
			return a.x==x;
		};
		int x;
		int y;
};

//bool equByY(edgeType a, edgeType b)
//{
//	if(a.y==b.y)
//		true;
//	false;
//}

//to transform into unique range
//bool equByX(edgeType a, edgeType b)
//{
//	if(a.x==b.x)
//		true;
//	false;
//}

template <typename T>
int hash(T c, T nFactor)
{
	//threeVector<int> c;
	while(c.x>nFactor.x)c.x-=nFactor.x;
	while(c.y>nFactor.y)c.y-=nFactor.y;
	while(c.z>nFactor.z)c.z-=nFactor.z;
	while(c.x<0)c.x+=nFactor.x;
	while(c.y<0)c.y+=nFactor.y;
	while(c.z<0)c.z+=nFactor.z;
	
	return c.x+c.y*nFactor.x+c.z*nFactor.x*nFactor.y;
}

template <typename T>
T unhash(int v, T nFactor)
{
	T c;
	c.x=v%nFactor.x;
	c.y=static_cast<int>(v/nFactor.x)%nFactor.y;
	c.z=static_cast<int>(v/(nFactor.x*nFactor.y))%nFactor.z;
	return c;
}

//I should double check this
//I double checked it. The radius '(dotProduct(a,a)*dotProduct(b,b)*dotProduct(c,c))/twoCrossSqr'
// is wrong or undersized. Can get radius from positions near circle.
template <typename T>
threeVector<double> circumscribe(T a, T b, T c)
{
	//std::cerr << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``" << std::endl;
	//std::cerr << a.x << '\t' << a.y << '\t' << a.z << std::endl;
	//std::cerr << b.x << '\t' << b.y << '\t' << b.z << std::endl;
	
	T crossProd=crossProduct(a,b);
	
	//std::cerr << crossProd.x << '\t' << crossProd.y << '\t' << crossProd.z << std::endl;
	
	
	double twoCrossSqr=2.0*(dotProduct(crossProd,crossProd));
	//double twoCrossSqr=2.0*(magnitude(crossProd));
	
	//std::cerr << "twoCrossSqr: " << twoCrossSqr << std::endl;
	
	threeVector<double> rFact;
	
	//radius
	//rFact.t=(dotProduct(a,a)*dotProduct(b,b)*dotProduct(c,c))/twoCrossSqr;
	//rFact.t=magnitude(a)*magnitude(b)*magnitude(c)/twoCrossSqr;
	
	//std::cerr << "twoCrossSqr: " << twoCrossSqr << std::endl;
	
	T aNeg=a;
	T bNeg=b;
	T cNeg=c;
	
	aNeg.x=-aNeg.x;
	aNeg.y=-aNeg.y;
	aNeg.z=-aNeg.z;
	
	bNeg.x=-bNeg.x;
	bNeg.y=-bNeg.y;
	bNeg.z=-bNeg.z;
	
	cNeg.x=-cNeg.x;
	cNeg.y=-cNeg.y;
	cNeg.z=-cNeg.z;
	
	//the offset scalers (threeVector r=rFact.x*a+rFact.y*b+rFact.z*c;)
	//rFact.x=dotProduct(b,b)*(dotProduct(a,cNeg))/twoCrossSqr;
	//rFact.y=dotProduct(cNeg,cNeg)*(dotProduct(aNeg,b))/twoCrossSqr;
	//rFact.z=dotProduct(a,a)*(dotProduct(c,bNeg))/twoCrossSqr;
	
	rFact.x=dotProduct(b,b)*(dotProduct(a,cNeg))/twoCrossSqr;
	rFact.y=dotProduct(cNeg,cNeg)*(dotProduct(aNeg,b))/twoCrossSqr;
	rFact.z=dotProduct(a,a)*(dotProduct(c,bNeg))/twoCrossSqr;
	return rFact;
	//std::cerr << "=====================================" << std::endl;
	//std::cin.get();
}



int main(int argc, char **argv)
{
	if(argc<8)
	{
		std::cerr << argv[0] << " name surfSpacing noTraverseRatio vertShellRatio vertSpacing dry type1 cType1 ...\n";
		std::cerr << "\tdry is for no output, if the surface spacing needs to be optimized. 0=false 1=true\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double surfSpacing, noTraverseRatio,vertSpacing,vertShellRatio;
	int dry;
	cmdArg >> surfSpacing >> noTraverseRatio >> vertShellRatio >> vertSpacing >> dry;
	double surfSpacingSqr=surfSpacing*surfSpacing;
	double vertSpacingSqr=vertSpacing*vertSpacing;
	double vertSpacingSqrMin=vertSpacing*vertSpacing*vertShellRatio*vertShellRatio;
	double vertSpacingSqrMax=vertSpacing*vertSpacing*(1.0+vertShellRatio)*(1.0+vertShellRatio);
	double acosSixty=acos(M_PI/3.0);
	
	double orthFactor=0.8019235/(vertSpacing*0.5*2.0);
	double crossFactor=0.8019235/(sqrt(vertSpacing*vertSpacing*0.25)*2.0);
	//double maxSphere=sqrt(4.0*2.0*surfSpacingSqr);//sqrt((2.0*surfSpacing)^2*dimension)
	//double maxSphereSqr=maxRad*maxRad;
	//double minSphereSqr=minRad*minRad;
	
	//surfSkip=1;//floor(2.0*maxRad/surfSpacing);
	
	//std::cout << surfSkip << std::endl;
	
	std::vector<int> types;
	std::vector<int> cTypes;
	
	for(int i=7;i<argc;i++)
	{
		int type, cType;
		cmdArg >> type >> cType;
		types.push_back(type);
		cTypes.push_back(cType);
	}
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//set time
	double time=0;

	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string newName("size_");
	newName+=argv[1];
	newName+=".dat";
		
	std::fstream sizeFile;
	sizeFile.open(newName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//grab particles, by type, for surfaces
	std::vector<twoVector<int> > pIndex;
	
	for(int i=0;i<System.readNParticles();i++)
	{
		bool match=false;
		for(int tIndex=0;!match && tIndex<types.size();tIndex++)
			match=(p[i].type==types[tIndex]);
		
		if(match)
		{
			twoVector<int> buf;
			buf.x=i;
			buf.y=-1;
			pIndex.push_back(buf);
		}
	}
	
	//molecule chains end to end if needed
	for(int i=0;i<cTypes.size();i++)
	{
		if(cTypes[i]<System.readNMolecules() && cTypes[i]>-1)
		{
			molecule<double,fourVector<int> > m=System.getMolecule()[cTypes[i]];
			if(m.readType()==CHAIN)
			{
				fourVector<int> *bond=m.getBonds();
				
				for(int j=0;j<m.readNBond();j++)
				{
					int start=bond[j].s[START];
					int nChains=bond[j].s[NCHAINS];
					int length=bond[j].s[CHAINLENGTH];
					
					for(int k=0;k<pIndex.size();k++)
					{
						int chain=(pIndex[k].x-start)/length;
						int seg=(pIndex[k].x-start)%length;
						if(seg==0)//its at the 0 end, use the length-1 end
							pIndex[k].y=(length-1)+chain*length;
						else if(seg==length-1)//its at the length-1 end, use the 0 end
							pIndex[k].y=chain*length;
						else//its in the middle, just use the 0 end
							pIndex[k].y=chain*length;
						if(pIndex[k].x<start || pIndex[k].x>start+nChains*length)
						{
							std::cerr << "Bad molecule!" << std::endl;
							throw 0;
						}
					}
				}
			}
		}
		else
		{
			std::cerr << "Bad molecule! #" << cTypes[i] << std::endl;
			throw 0;
		}
	}
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<int> sFlagOrigin;
	int nSurfaces=0;
	int frame=0;
	double avgNeighbors=0;
	
	while(xyzFile.load())
	{
		time+=System.readStoreInterval();
		std::cerr << time << '\t';
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		//size=System.readSize();
		//set up a quick search using the minimum surface spacing
		threeVector<int> surfN;
		surfN.x=static_cast<int>(size.x/surfSpacing);
		surfN.y=static_cast<int>(size.y/surfSpacing);
		surfN.z=static_cast<int>(size.z/surfSpacing);
		threeVector<double> surfR;
		surfR.x=size.x/static_cast<double>(surfN.x);
		surfR.y=size.y/static_cast<double>(surfN.y);
		surfR.z=size.z/static_cast<double>(surfN.z);
		
		typedef std::vector<position<double> > cellType;//fastest is vector, deque is a little slower
		
		std::map<int,std::vector<int> > pFilt;
		
		std::map<int,cellType > cells,etoe;
		
		threeVector<double> zero=0;
		
		#define unusedFlag -1
		#define noTraverseFlag -2
		#define traverseFlag 0
		
		for(int j=0;j<pIndex.size();j++)
		{
			int i=pIndex[j].x;
			threeVector<int> c;
			c.x=static_cast<int>(p[i].x/surfR.x)%surfN.x;
			c.y=static_cast<int>(p[i].y/surfR.y)%surfN.y;
			c.z=static_cast<int>(p[i].z/surfR.z)%surfN.z;
			int hash=c.x+c.y*surfN.x+c.z*surfN.x*surfN.y;
			//int hash=cellHash(p[i]);
			std::map<int,cellType >::iterator it=cells.lower_bound(hash);
			if(it!=cells.end() && it->first==hash)//!(cells.key_comp()(it->first, hash)))
			{
				//std::cerr << it->first << '\t' << hash
				
				//push operation
				it->second.push_back(p[i]);
				etoe[hash].push_back(p[pIndex[j].y]);
				pFilt[hash].push_back(unusedFlag);
			}
			else
			{
				//initialize operation
				cellType buf;
				buf.push_back(p[i]);
				
				cells.insert(it, std::map<int,cellType >::value_type(hash,buf));
				
				buf[0]=p[pIndex[j].y];
				etoe[hash]=buf;
				
				std::vector<int> sFlags;
				sFlags.push_back(unusedFlag);
				pFilt[hash]=sFlags;
				
			}
		}
		
		std::cerr << "using " << pIndex.size() << " points for surfaces" << std::endl;
		
		if(frame==0)
		{
			std::cerr << "Determining average number of neighbors per particle..." << std::endl;
			long int nNeighbors=0;
			//obtaining surfaces
			for(auto &cell:cells)
			{
				int currentHash=cell.first;
				cellType &cellPoints=cell.second;
				
				threeVector<int> iCoor=unhash(currentHash,surfN);
				
				for(auto &aP:cellPoints)
				{
					//if(oldHist.size()>0)
					std::vector<double> dist;
					double minDist=0;//size.x*size.x+size.y*size.y+size.z*size.z;
					//range search
					for(int k=0;k<27;k++)
					{
						threeVector<int> nCoor=iCoor;
						nCoor.x+=k%3-1;
						nCoor.y+=int(k/3)%3-1;
						nCoor.z+=int(k/9)-1;
						int neighborHash=hash(nCoor,surfN);
						auto it=cells.lower_bound(neighborHash);
						if(it!=cells.end() && it->first==neighborHash)
						{
							cellType &neigPoints=it->second;
							for(auto &neP:neigPoints)
							{
								threeVector<double> d;
								d.x=neP.x-aP.x;
								d.y=neP.y-aP.y;
								d.z=neP.z-aP.z;
								
								while(d.x>size.x/2.0)d.x-=size.x;
								while(d.y>size.y/2.0)d.y-=size.y;
								while(d.z>size.z/2.0)d.z-=size.z;
								while(d.x<-size.x/2.0)d.x+=size.x;
								while(d.y<-size.y/2.0)d.y+=size.y;
								while(d.z<-size.z/2.0)d.z+=size.z;
								
								double dr=d.x*d.x+d.y*d.y+d.z*d.z;
								if(dr<surfSpacingSqr && dr!=0)// && dotProduct(normalSum,otherSum)>sqrt(2.0)/2.0)
									nNeighbors++;
							}
						}
					}
				}
			}
			avgNeighbors=nNeighbors/pIndex.size();
			std::cerr << "Average number of neighbors: " << avgNeighbors << std::endl;
			std::cerr << "Using " << avgNeighbors*noTraverseRatio << " neighbors for no traverse condition!" << std::endl;
		}
		
		std::cerr << "using " << pIndex.size() << " points for surface discrimination" << std::endl;
		
		std::vector< std::vector< position<double> > > surfaces, surfacesETOE;
		
		
		//obtaining surfaces
		for(auto pFIT=pFilt.begin();pFIT!=pFilt.end();pFIT++)//a cell
		for(auto pTest=pFIT->second.begin();pTest!=pFIT->second.end();pTest++)//a particle in cell
		{
			//unaccounted surface
			if(*pTest==unusedFlag)
			{
				//check if it already has too many neighbors
				int surfaceHash=pFIT->first;
				
				threeVector<int> iCoor=unhash(surfaceHash,surfN);
				
				int cellOffset=std::distance(pFIT->second.begin(),pTest);
				
				position<double> sPoint=cells[surfaceHash][cellOffset];
				
				//integrate nearby points
				int nNeighbors=0;
				for(int k=0;k<27;k++)
				{
					threeVector<int> nCoor=iCoor;
					nCoor.x+=k%3-1;
					nCoor.y+=int(k/3)%3-1;
					nCoor.z+=int(k/9)-1;
					int neighborHash=hash(nCoor,surfN);
					
					auto it=cells.lower_bound(neighborHash);
					if(it!=cells.end() && it->first==neighborHash)
					{
						for(auto &neP:it->second)
						{
							//position<double> neP=it->second;
							threeVector<double> d;
							d.x=neP.x-sPoint.x;
							d.y=neP.y-sPoint.y;
							d.z=neP.z-sPoint.z;
							
							while(d.x>size.x/2.0)d.x-=size.x;
							while(d.y>size.y/2.0)d.y-=size.y;
							while(d.z>size.z/2.0)d.z-=size.z;
							while(d.x<-size.x/2.0)d.x+=size.x;
							while(d.y<-size.y/2.0)d.y+=size.y;
							while(d.z<-size.z/2.0)d.z+=size.z;
							
							if(d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
								nNeighbors++;
						}
					}
				}
				
				if(nNeighbors>avgNeighbors*noTraverseRatio)
					*pTest=noTraverseFlag;
				else
					*pTest=surfaces.size();
				
				std::deque<twoVector<int> > stack;
				
				twoVector<int> currentParticle;
				currentParticle.x=surfaceHash;
				currentParticle.y=cellOffset;
				
				if(*pTest!=noTraverseFlag)
					stack.push_back(currentParticle);
				
				//std::vector<std::map<int, position<double> >::iterator> surface;
				cellType surface, surfaceETOE;
				
				while(stack.size()>0)
				{
					currentParticle=stack.back();
					
					stack.pop_back();
					
					surfaceHash=currentParticle.x;
					cellOffset=currentParticle.y;
					
					position<double> aP=cells[surfaceHash][cellOffset];
					position<double> aETOE=etoe[surfaceHash][cellOffset];
					
					surface.push_back(aP);
					surfaceETOE.push_back(aETOE);
					
					iCoor=unhash(surfaceHash,surfN);
					
					std::vector< twoVector<int> > candidates;
					
					//surface search while integrating nearby points
					nNeighbors=0;
					for(int k=0;k<27;k++)
					{
						threeVector<int> nCoor=iCoor;
						nCoor.x+=k%3-1;
						nCoor.y+=int(k/3)%3-1;
						nCoor.z+=int(k/9)-1;
						int neighborHash=hash(nCoor,surfN);
						auto it=pFilt.lower_bound(neighborHash);
						if(it!=pFilt.end() && it->first==neighborHash)
						{
							for(auto neFilt=it->second.begin();neFilt!=it->second.end();neFilt++)
							{
								int neOffset=std::distance(it->second.begin(),neFilt);
								position<double> neP=cells[neighborHash][neOffset];
								threeVector<double> d;
								d.x=neP.x-aP.x;
								d.y=neP.y-aP.y;
								d.z=neP.z-aP.z;
								
								while(d.x>size.x/2.0)d.x-=size.x;
								while(d.y>size.y/2.0)d.y-=size.y;
								while(d.z>size.z/2.0)d.z-=size.z;
								while(d.x<-size.x/2.0)d.x+=size.x;
								while(d.y<-size.y/2.0)d.y+=size.y;
								while(d.z<-size.z/2.0)d.z+=size.z;
								
								double dr=d.x*d.x+d.y*d.y+d.z*d.z;
								
								if(dr<surfSpacingSqr && (!(neOffset==cellOffset) != !(neighborHash==surfaceHash)))
								{
									
									
									
									nNeighbors++;
									if(*neFilt==unusedFlag)
									{
										twoVector<int> neParticle;
										neParticle.x=it->first;
										neParticle.y=neOffset;
										candidates.push_back(neParticle);
										
										
										
										//std::cout << neP.x << ' ' << neP.y << ' ' << neP.z << std::endl;
									}
								}
							}
						}
					}
					//std::cerr << stack.size() << ' ' << candidates.size() << std::endl;
					if(nNeighbors<avgNeighbors*noTraverseRatio)
					{
						for(auto &candidate:candidates)
						{
							stack.push_back(candidate);
							pFilt[candidate.x][candidate.y]=surfaces.size();
						}
					}
					else if(surface.size()>1)
					{
						pFilt[surfaceHash][cellOffset]=surfaces.size();
					}
					//std::cerr << stack.size() << std::endl;
					//std::cin.get();
				}
				if(surface.size()>1)
				{
					surfaces.push_back(surface);
					surfacesETOE.push_back(surfaceETOE);
				}
				//nSurfaces++;
			}
		}
		std::cerr << "Found " << surfaces.size() << " surfaces." << std::endl;
		
		/*
		if(surfaces.size()==1)
		{
			std::fstream dumpFile;
			dumpFile.open("typeProfile.xyz", std::ios::app | std::ios::out);
			dumpFile << surfaces[0].size() << "\ntest\n";
			for(auto &aP:surfaces[0])
				dumpFile << aP.type << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << '\n';
			dumpFile.close();
			throw 0;
		}*/
		//set up a quick search using the minimum surface spacing
		threeVector<int> vertN;
		vertN.x=static_cast<int>(size.x/(vertSpacing));
		vertN.y=static_cast<int>(size.y/(vertSpacing));
		vertN.z=static_cast<int>(size.z/(vertSpacing));
		threeVector<double> vertR;
		vertR.x=size.x/static_cast<double>(vertN.x);
		vertR.y=size.y/static_cast<double>(vertN.y);
		vertR.z=size.z/static_cast<double>(vertN.z);
		
		typedef std::vector<threeVector<double> > normalType;
		
		std::vector<std::map<int,cellType > > surfCells;
		std::vector<std::map<int,normalType > > surfNormals;
		std::vector<threeVector<double> > centerOfMass;
		
		auto sETOE=surfacesETOE.begin();
		
		for(auto& surface:surfaces)
		{
			threeVector<double> com;
			com=0.0;
			//if(surface.size()>10)
			std::cerr << "surface size: " << surface.size() << std::endl;
			std::map<int,cellType > sCells;
			std::map<int,normalType > normals;
			std::vector<position<double> >::iterator itETOE=sETOE->begin();
			//std::sort(surface.begin(),surface.end());
			for(auto& pI:surface)
			{
				threeVector<double> normal;
				normal.x=itETOE->x-pI.x;
				normal.y=itETOE->y-pI.y;
				normal.z=itETOE->z-pI.z;
				normal.x-=(normal.x>size.x/2.0)?size.x:0;
				normal.y-=(normal.y>size.y/2.0)?size.y:0;
				normal.z-=(normal.z>size.z/2.0)?size.z:0;
				normal.x-=(normal.x<-size.x/2.0)?size.x:0;
				normal.y-=(normal.y<-size.y/2.0)?size.y:0;
				normal.z-=(normal.z<-size.z/2.0)?size.z:0;
				normal=unitVector(normal);
				if(normal.x!=normal.x)
				{
					std::cerr << itETOE->x << ' ' << pI.x << std::endl;
					throw 0;
				}
				
				threeVector<int> c;
				c.x=static_cast<int>(pI.x/vertR.x)%vertN.x;
				c.y=static_cast<int>(pI.y/vertR.y)%vertN.y;
				c.z=static_cast<int>(pI.z/vertR.z)%vertN.z;
				int hash=c.x+c.y*vertN.x+c.z*vertN.x*vertN.y;
				
				std::map<int,cellType >::iterator it=sCells.lower_bound(hash);
				if(it!=sCells.end() && it->first==hash)//!(cells.key_comp()(it->first, hash)))
				{
					//push operation
					it->second.push_back(pI);
					
					normals[hash].push_back(normal);
				}
				else
				{
					//initialize operation
					cellType buf;
					buf.push_back(pI);
					//cells[hash]=buf;
					sCells.insert(it, std::map<int,cellType >::value_type(hash,buf));
					
					
					normalType normalBuf;
					normalBuf.push_back(normal);
					normals[hash]=normalBuf;
				}
				com.x+=pI.x;
				com.y+=pI.y;
				com.z+=pI.z;
				
				itETOE++;
			}
			//nSurface++;
			surfCells.push_back(sCells);
			surfNormals.push_back(normals);
			sETOE++;
			com.x/=surface.size();
			com.y/=surface.size();
			com.z/=surface.size();
			centerOfMass.push_back(com);
		}
		
		//if(surfaces.size()<4)
		{
			for(auto &sCells:surfCells)
				if(sCells.size()>10)
					std::cerr << "using " << sCells.size() << " points for surface curvature" << std::endl;
		}
		
		//examined surfaces
		typedef std::map<int, position<double> >  filtType;
		std::vector<filtType> vFilt;
		
		//typedef std::map<int,threeVector<double> > filtNormalType;
		//std::vector<filtNormalType> vFiltNormals;
		
		#define unusedPoint -1
		#define usedPoint 0;
		int surf=0;
		for(auto& sCells:surfCells)
		{
			filtType filt;
			//filtNormalType filtNormal;
			//subsample our pIndex
			//c++11 syntax reduces this:
			//std::map<int,cellType >::iterator it
			for(auto& it: sCells)
			{
				int gridPos=it.first;
				
				position<double> pAvg;
				pAvg.x=0;
				pAvg.y=0;
				pAvg.z=0;
				
				//std::vector<int>::iterator pI
				for(auto& pI: it.second)
				{
					pAvg.x+=pI.x;
					pAvg.y+=pI.y;
					pAvg.z+=pI.z;
					//std::cout << surf << ' ' << pI.x << ' ' << pI.y << ' ' << pI.z << std::endl;
				}
				pAvg.x/=static_cast<double>(it.second.size());
				pAvg.y/=static_cast<double>(it.second.size());
				pAvg.z/=static_cast<double>(it.second.size());
				pAvg.type=-1;
				
				//filtNormal[gridPos]=zero;
				
				filt[gridPos]=pAvg;
			}
			vFilt.push_back(filt);
			//vFiltNormals.push_back(filtNormal);
			surf++;
		}
		
		//throw 0;
		std::vector< double > surfaceCurvature;
		
		if(vFilt.size()>1)
			std::cout << time;
		
		//try every surface to get surface normals
		for(int s=0;s<vFilt.size() && vFilt.size()>1;s++)
		{
			
			int nSpacing=100;
		
			std::vector<threeVector<double> > histogram;
			//std::vector<twoVector<double> > radDist;
			std::vector<std::vector<double> > gradSquareValues;
			//std::vector<double > gradSquareValues;
			//double dGrad=0.0001;//smallest distinguishible gradSquare
			//double dRad=0.1;
			double patchInc=vertSpacing*sqrt(2.0)/(static_cast<double>(nSpacing+1)*2.0);
			//double rSpacingScalar=surfSpacing/M_PI;
			
			for(int i=0;i<nSpacing+1;i++)
			{
				threeVector<double> hist;
				hist.x=static_cast<double>(i)*patchInc+patchInc*0.5;
				hist.y=0;
				hist.z=0;
				histogram.push_back(hist);
				
				gradSquareValues.push_back(std::vector<double>());
			}
			
			int nFound=0;
		
			double gradSquared=0;
			#define vertexPoint 0
			#define edgePoint 1
			#define ignorePoint 2
			auto& filt=vFilt[s];//for nearest neighbor traversal
			auto& sCells=surfCells[s];//for vertex candidates
			auto& normals=surfNormals[s];
			auto& surface=surfaces[s];
			//auto& filtNormal=vFiltNormals[s];
			
			int iSurface=0;
			
			int patch=0;
			int nElements=0;
			//unaccounted surface
			for(auto element=filt.begin();element!=filt.end() && filt.size()>5;element++)
			{
				if((nElements++)%10==0)
					std::cerr << element->first << std::endl;
				//unaccounted surface element?
				if(element->second.type==unusedPoint)
				{
					threeVector<double> oneZ=0;
					oneZ.z=1.0;
					
					element->second.type=usedPoint;
					
					//stack.push_back(element);
					
					
					//rundown:
					//get the local alignment
					//	gather nearby normal vectors (filtNormal->nList)
					//	generate local normal (avg(nList))
					//check that the local alignment is continuous
					//keep up with local alignment
					//	gather local center's of mass (filt->bcP)
					//	maintain stack (stack.push_back)
					//	generate new local normal (avg(sum(crossProduct)))
					//now align the clustered group
					//	gather nearby particles (sCells->nearbyGroups)
					//	generate per particle normal (avg(sum(crossProduct)))
					//get the projection matrix
					//get the curvature
					
					
					
					auto currentPoint=element;
					
					int currentHash=currentPoint->first;
					
					threeVector<int> iCoor=unhash(currentHash,vertN);
					
					currentPoint->second.type=usedPoint;
					
					
					auto& aGroup=sCells[currentHash];
					
					//filtNormal[currentHash]=normalSum;
					
					//now align the clustered group
					
					//auto normal=normals[currentHash].begin();//iterator
					std::vector<position<double> > nearbyGroups;
					std::vector<threeVector<double> > nearbyNormals;
					
					for(int k=0;k<27;k++)
					{
						threeVector<int> nCoor=iCoor;
						nCoor.x+=k%3-1;
						nCoor.y+=int(k/3)%3-1;
						nCoor.z+=int(k/9)-1;
						int neighborHash=hash(nCoor,vertN);
						
						auto bGroup=sCells.lower_bound(neighborHash);
						if(bGroup!=sCells.end() && bGroup->first==neighborHash && neighborHash!=currentHash)
						{
							for(auto& point:bGroup->second)
								nearbyGroups.push_back(point);
							for(auto& normal:normals[neighborHash])
								nearbyNormals.push_back(normal);
						}
					}
					
					auto normalGroup=normals[currentHash];
					
					auto myNormal=normalGroup.begin();
					
					//double patchSize=vertSpacing*sqrt(2)/2.0;
					
					//std::vector<double> inDistance,outDistance;
					
					
					
					//std::cerr << nearbyGroups.size() << std::endl;
					int subpatch=0;
					//do this for every particle in aGroup
					for(auto& aP:aGroup)
					{
						double minGradSquare=1.0;
						double minHistIndex=0;
						for(double patchSize=vertSpacing*sqrt(2)/2.0-patchInc;patchSize>1.0;patchSize-=patchInc)
						{
							threeVector<double> normalSum=0;
							int nNormals=0;
							for(int i=aGroup.size()-1;i>=0;--i)
							{
								threeVector<double> d;
								d.x=aGroup[i].x-aP.x;
								d.y=aGroup[i].y-aP.y;
								d.z=aGroup[i].z-aP.z;
								d.x-=(d.x>size.x/2.0)?size.x:0;
								d.y-=(d.y>size.y/2.0)?size.y:0;
								d.z-=(d.z>size.z/2.0)?size.z:0;
								d.x+=(d.x<-size.x/2.0)?size.x:0;
								d.y+=(d.y<-size.y/2.0)?size.y:0;
								d.z+=(d.z<-size.z/2.0)?size.z:0;
								if(d.x*d.x+d.y*d.y+d.z*d.z<patchSize*patchSize*4.0 && dotProduct(normalGroup[i],*myNormal)>0)
								{
									normalSum+=normalGroup[i];
									nNormals++;
								}
							}
							for(int i=nearbyGroups.size()-1;i>=0;--i)
							{
								threeVector<double> d;
								d.x=nearbyGroups[i].x-aP.x;
								d.y=nearbyGroups[i].y-aP.y;
								d.z=nearbyGroups[i].z-aP.z;
								d.x-=(d.x>size.x/2.0)?size.x:0;
								d.y-=(d.y>size.y/2.0)?size.y:0;
								d.z-=(d.z>size.z/2.0)?size.z:0;
								d.x+=(d.x<-size.x/2.0)?size.x:0;
								d.y+=(d.y<-size.y/2.0)?size.y:0;
								d.z+=(d.z<-size.z/2.0)?size.z:0;
								if(d.x*d.x+d.y*d.y+d.z*d.z<patchSize*patchSize*4.0 && dotProduct(nearbyNormals[i],*myNormal)>0)
								{
									normalSum+=nearbyNormals[i];
									nNormals++;
								}
							}
							
							normalSum/=static_cast<double>(nNormals);
							
							//get the projection matrix, ROW_MAJOR
							threeVector<threeVector<double> > M;
							//M=I-(aP*aP^T)/aP^2
							//position<double> aPu=unitVector(aP);
							threeVector<double> aPu=unitVector(normalSum);//(*normal++);
							
							
							
							//aPu.x=-1.0/sqrt(3.0);
							//aPu.y=-1.0/sqrt(3.0);
							//aPu.z=-1.0/sqrt(3.0);
							
							//threeVector<double> aPuX=aPu,aPuZ=aPu;
							//aPuX.y=0;
							//aPuZ.x=0;
							//aPuX=unitVector(aPuX);
							//aPuZ=unitVector(aPuZ);
							//aPuX.x=1.0;
							//aPuX.y=0;
							//aPuX.z=0;
							//aPuY.x=0;
							//aPuY.y=1.0;
							//aPuY.z=0;
							//std::cout << "3 " << aP.x << ' ' << aP.y << ' ' << aP.z << std::endl;
							//std::cout << "Input: " << aPu.x << ' ' << aPu.y << ' ' << aPu.z << std::endl;
							//std::cout << "With magnitude: " << magnitude(aPu) << std::endl;
							//std::cout << "We want it here: " << oneZ.x << ' ' << oneZ.y << ' ' << oneZ.z << std::endl;
							//double aPuDotOneZ=dotProduct(aPu,oneZ);
							//double sinSqr=sin(acos(aPuDotOneZ)/2.0)*sin(acos(aPuDotOneZ)/2.0);//(1.0-aPuDotOneZ)/2.0;
							//double sinCos=sin(acos(aPuDotOneZ)/2.0)*cos(acos(aPuDotOneZ)/2.0);//sqrt(1.0+aPuDotOneZ*aPuDotOneZ)/2.0;
							//double qSinSqr=(1.0-aPuDotOneZ)/2.0;//q^2=sin(theta/2.0)^2
							//double qSin=sqrt(1.0+aPuDotOneZ*aPuDotOneZ)/2.0;//q=sin(theta/2.0)
							//double qCos=sqrt((1.0+aPuDotOneZ)/2.0);//q=cos(theta/2.0)
							
							
							
							//if(subpatch==0)
							//{
							//	std::cout << s*3 << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << std::endl;
							//	std::cout << s*3+1 << ' ' << aP.x+aPu.x*5.0 << ' ' << aP.y+aPu.y*5.0 << ' ' << aP.z+aPu.z*5.0 << std::endl;
							//}
							if(dotProduct(aPu,oneZ)<0)
							{
								oneZ.z=-oneZ.z;
							}
							threeVector<double> aPuZ=crossProduct(aPu,oneZ);
							
							double w=sqrt(1+dotProduct(aPu,oneZ));
							double b = aPuZ.x;
							double c = aPuZ.y;
							double d = aPuZ.z;
							double a = w;
							
							
							
							//stolen from boost, boosted from boost?
							double aa = a*a;
							double ab = a*b;
							double ac = a*c;
							double ad = a*d;
							double bb = b*b;
							double bc = b*c;
							double bd = b*d;
							double cc = c*c;
							double cd = c*d;
							double dd = d*d;
							
							double norme_carre = aa+bb+cc+dd;
							
							M.x.x=(aa + bb - cc - dd)/norme_carre;
							M.x.y=2 * (-ad + bc)/norme_carre;
							M.x.z=2 * (ac + bd)/norme_carre;
							M.y.x= 2 * (ad + bc)/norme_carre;
							M.y.y=(aa - bb + cc - dd)/norme_carre;
							M.y.z=2 * (-ab + cd)/norme_carre;
							M.z.x=2 * (-ac + bd)/norme_carre;
							M.z.y=2 * (ab + cd)/norme_carre;
							M.z.z=(aa - bb - cc + dd)/norme_carre;
							
							//four rotation axii (should work, don't need to do this)
							//M.x.x=1.0-qSqrJ-qSqrK;
							//M.x.y=qIJ-qKR;
							//M.x.z=qIK+qJR;
							//M.y.x=qIJ+qKR;
							//M.y.y=1.0-qSqrI-qSqrK;
							//M.y.z=qJK-qJR;
							//M.z.x=qIK-qKR;
							//M.z.y=qJK+qIR;
							//M.z.z=1.0-qSqrI-qSqrJ;
							
							//single axis (works)
							//M.x.x=1.0-2.0*sinSqr;
							//M.x.y=1.0*(sinSqr-sinCos);
							//M.x.z=1.0*(sinSqr+sinCos);
							//M.y.x=(sinSqr+sinCos);
							//M.y.y=1.0-2.0*sinSqr;
							//M.y.z=1.0*(sinSqr-sinCos);
							//M.z.x=1.0*(sinSqr-sinCos);
							//M.z.y=1.0*(sinSqr+sinCos);
							//M.z.z=1.0-2.0*sinSqr;
							
							//std::cout << "[ " << M.x.x << ' ' << M.x.y << ' ' << M.x.z << " ]" << std::endl;
							//std::cout << "[ " << M.y.x << ' ' << M.y.y << ' ' << M.y.z << " ]" << std::endl;
							//std::cout << "[ " << M.z.x << ' ' << M.z.y << ' ' << M.z.z << " ]" << std::endl;
							
							
							//threeVector<double> aPuT;
							//aPuT.x=aPu.x*M.x.x+aPu.y*M.x.y+aPu.z*M.x.z;
							//aPuT.y=aPu.x*M.y.x+aPu.y*M.y.y+aPu.z*M.y.z;
							//aPuT.z=aPu.x*M.z.x+aPu.y*M.z.y+aPu.z*M.z.z;
							//if(dotProduct(aPuT,oneZ)<sqrt(2.0)/2.0)
							//{
							//	std::cerr << "It came out here: " << aPuT.x << ' ' << aPuT.y << ' ' << aPuT.z << std::endl;
							//	std::cerr << "With magnitude: " << magnitude(aPuT) << std::endl;
							//	std::cerr << "from here: " << aPu.x << ' ' << aPu.y << ' ' << aPu.z << std::endl;
							//}
							//std::cout << "0 " << 0 << ' ' << 0 << ' ' << 0 << std::endl;
							
							//aPuZ=crossProduct(aPuT,oneZ);
							//w=sqrt(1+dotProduct(aPuT,oneZ));
							//b = aPuZ.x;
							//c = aPuZ.y;
							// d = aPuZ.z;
							// a = w;
							
							
							
							//stolen from boost, boosted from boost?
							// aa = a*a;
							// ab = a*b;
							// ac = a*c;
							// ad = a*d;
							// bb = b*b;
							// bc = b*c;
							// bd = b*d;
							// cc = c*c;
							// cd = c*d;
							// dd = d*d;
							
							// norme_carre = aa+bb+cc+dd;
							
							//M.x.x=(aa + bb - cc - dd)/norme_carre;
							//M.x.y=2 * (-ad + bc)/norme_carre;
							//M.x.z=2 * (ac + bd)/norme_carre;
							//M.y.x= 2 * (ad + bc)/norme_carre;
							//M.y.y=(aa - bb + cc - dd)/norme_carre;
							//M.y.z=2 * (-ab + cd)/norme_carre;
							//M.z.x=2 * (-ac + bd)/norme_carre;
							//M.z.y=2 * (ab + cd)/norme_carre;
							//M.z.z=(aa - bb - cc + dd)/norme_carre;
							
							//std::cout << "[ " << M.x.x << ' ' << M.x.y << ' ' << M.x.z << " ]" << std::endl;
							//std::cout << "[ " << M.y.x << ' ' << M.y.y << ' ' << M.y.z << " ]" << std::endl;
							//std::cout << "[ " << M.z.x << ' ' << M.z.y << ' ' << M.z.z << " ]" << std::endl;
							
							
							//threeVector<double> aPuTT;
							//aPuTT.x=aPuT.x*M.x.x+aPuT.y*M.x.y+aPuT.z*M.x.z;
							//aPuTT.y=aPuT.x*M.y.x+aPuT.y*M.y.y+aPuT.z*M.y.z;
							//aPuTT.z=aPuT.x*M.z.x+aPuT.y*M.z.y+aPuT.z*M.z.z;
							//std::cout << "It came out here: " << aPuTT.x << ' ' << aPuTT.y << ' ' << aPuTT.z << std::endl;
							//std::cout << "With magnitude: " << magnitude(aPuTT) << std::endl;
							//std::cout << "0 " << 0 << ' ' << 0 << ' ' << 0 << std::endl;
							
							//threeVector<double> buf;
							//buf.x=aPu.
							//threeVector<double> theta;
							//theta.x=acos(aPu
							
							//our z oriented patch (heights) I need to make this behavior work
							//threeVector<threeVector<double> > zPatch=0;
							//number of points used in this patch
							threeVector<threeVector<double> > nZPatch;
							//our z oriented patch (heights)
							threeVector<threeVector<double> > zPatch;
							
							for(int a=0;a<3;a++)
								for(int b=0;b<3;b++)
									zPatch.s[a].s[b]=0;
							for(int a=0;a<3;a++)
								for(int b=0;b<3;b++)
									nZPatch.s[a].s[b]=0;
								
							/*
							for(auto &fG:aGroup)
							{
								threeVector<double> dG;
								dG.x=fG.x-aP.x;
								dG.y=fG.y-aP.y;
								dG.z=fG.z-aP.z;
								
								dG.x-=(dG.x>size.x/2.0)?size.x:0;
								dG.y-=(dG.y>size.y/2.0)?size.y:0;
								dG.z-=(dG.z>size.z/2.0)?size.z:0;
								
								threeVector<double> jG;
								jG.x=dG.x*M.x.x+dG.y*M.x.y+dG.z*M.x.z;
								jG.y=dG.x*M.y.x+dG.y*M.y.y+dG.z*M.y.z;
								jG.z=dG.x*M.z.x+dG.y*M.z.y+dG.z*M.z.z;
								
								jG.z*=oneZ.z;
								
								for(auto &hG:aGroup)
								{
									
									threeVector<double> eG;
									
									eG.x=hG.x-aP.x;
									eG.y=hG.y-aP.y;
									eG.z=hG.z-aP.z;
									
									eG.x-=(eG.x>size.x/2.0)?size.x:0;
									eG.y-=(eG.y>size.y/2.0)?size.y:0;
									eG.z-=(eG.z>size.z/2.0)?size.z:0;
									
									threeVector<double> iG;
									
									iG.x=eG.x*M.x.x+eG.y*M.x.y+eG.z*M.x.z;
									iG.y=eG.x*M.y.x+eG.y*M.y.y+eG.z*M.y.z;
									iG.z=eG.x*M.z.x+eG.y*M.z.y+eG.z*M.z.z;
									
									iG.z*=oneZ.z;
									
									threeVector<double> dV;
									dV.x=iG.x-jG.x;
									dV.y=iG.y-jG.y;
									dV.z=iG.z-jG.z;
									//std::cout << "Hello! ";
									std::cout << sqrt(dV.x*dV.x+dV.y*dV.y+dV.z*dV.z) << ' ';
									
									dV.x=dG.x-eG.x;
									dV.y=dG.y-eG.y;
									dV.z=dG.z-eG.z;
									
									std::cout << sqrt(dV.x*dV.x+dV.y*dV.y+dV.z*dV.z) << std::endl;
								}
								
							}*/
								
							//auto mu=Mu.begin();
							//auto kS=K.begin();
							//transform the neighbors by the projection matrix and place on grid
							
							for(auto &nG:nearbyGroups)
							{
								threeVector<double> sG,qG;
								
								qG.x=nG.x-aP.x;
								qG.y=nG.y-aP.y;
								qG.z=nG.z-aP.z;
								
								qG.x-=(qG.x>size.x/2.0)?size.x:0;
								qG.y-=(qG.y>size.y/2.0)?size.y:0;
								qG.z-=(qG.z>size.z/2.0)?size.z:0;
								
								sG.x=qG.x*M.x.x+qG.y*M.x.y+qG.z*M.x.z;
								sG.y=qG.x*M.y.x+qG.y*M.y.y+qG.z*M.y.z;
								sG.z=qG.x*M.z.x+qG.y*M.z.y+qG.z*M.z.z;
								
								
								
								sG.z*=oneZ.z;
								
								//std::cout << sqrt(sG.x*sG.x+sG.y*sG.y+sG.z*sG.z) << std::endl;
								
								//sG=qG;
								
								//Rotate with respect to the x axis:
								//sG.y=nG.y*cos(theta.x)-nG.z*sin(theta.x);
								//sG.z=nG.y*sin(theta.x)+nG.z*cos(theta.x);
								
								//Rotate with respect to the y axis:
								//sG.x=nG.x*cos(theta.y)+nG.z*sin(theta.y);
								//sG.z=-nG.x*sin(theta.y)+nG.z*cos(theta.y);
								
								//Rotate with respect to the z axis:
								//sG.x=nG.x*cos(theta.z)-nG.y*sin(theta.z);
								//sG.y=nG.x*sin(theta.z)+nG.y*cos(theta.z);
								
								threeVector<int> patchCoor=-1.0;
								if(sG.x>-1.5*patchSize && sG.x<1.5*patchSize)
									patchCoor.x=static_cast<int>((sG.x+patchSize*1.5)/(patchSize));
								if(sG.y>-1.5*patchSize && sG.y<1.5*patchSize)
									patchCoor.y=static_cast<int>((sG.y+patchSize*1.5)/(patchSize));
								//if(sG.z>0 && sG.z<1.5*vertSpacing)
								//	patchCoor.z=static_cast<int>(sG.z/vertSpacing);
								
								if(qG.x*qG.x+qG.y*qG.y+qG.z*qG.z<4.0*patchSize*patchSize)
								if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
								{
									zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
									nZPatch.s[patchCoor.x].s[patchCoor.y]++;
									//if(patch==40 && subpatch==0)//patch==0 || patch==5 || patch==20)
										//std::cout << patch << ' ' << sG.x << ' ' << sG.y << ' ' << sG.z << std::endl;
										//std::cout << patchCoor.x+3*patchCoor.y+9 << ' ' << sG.x << ' ' << sG.y << ' ' << sG.z << std::endl;
										//std::cout << "2 " << nG.x << ' ' << nG.y << ' ' << nG.z << std::endl;
								}
								
								//threeVector<double> sN=unitVector(sG);
								//threeVector<double> sNdotaPu=dotProduct(sN,aPu);
								//double denominator=sqrt(sNdotaPu*sNdotaPu+1.0)*sqrt(sG.x*sG.x+sG.y*sG.y);
								//(*kS++)=sNdotaPu/denominator;
								
								//yxy
								//rrx
								//double xx=sG.x*sG.x;
								//double yy=sG.y*sG.y;
								//double rProjSqr=xx+yy;
								//cos^2
								//mu->x=xx/rProjSqr;
								//2*sin*cos
								//mu->y=2.0*sG.x*sG.y/rProjSqr;
								//sin^2
								//mu->z=yy/rProjSqr;
								
							}
							for(auto &nG:aGroup)
							{
								threeVector<double> sG,qG;
								
								qG.x=nG.x-aP.x;
								qG.y=nG.y-aP.y;
								qG.z=nG.z-aP.z;
								
								qG.x-=(qG.x>size.x/2.0)?size.x:0;
								qG.y-=(qG.y>size.y/2.0)?size.y:0;
								qG.z-=(qG.z>size.z/2.0)?size.z:0;
								
								sG.x=qG.x*M.x.x+qG.y*M.x.y+qG.z*M.x.z;
								sG.y=qG.x*M.y.x+qG.y*M.y.y+qG.z*M.y.z;
								sG.z=qG.x*M.z.x+qG.y*M.z.y+qG.z*M.z.z;
								
								sG.z*=oneZ.z;
								//sG=qG;
								
								threeVector<int> patchCoor=-1.0;
								if(sG.x>-1.5*patchSize && sG.x<1.5*patchSize)
									patchCoor.x=static_cast<int>((sG.x+patchSize*1.5)/(patchSize));
								if(sG.y>-1.5*patchSize && sG.y<1.5*patchSize)
									patchCoor.y=static_cast<int>((sG.y+patchSize*1.5)/(patchSize));
								//if(sG.z>0 && sG.z<1.5*vertSpacing)
								//	patchCoor.z=static_cast<int>(sG.z/vertSpacing);
								
								if(qG.x*qG.x+qG.y*qG.y+qG.z*qG.z<4.0*patchSize*patchSize)
								if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
								{
									zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
									nZPatch.s[patchCoor.x].s[patchCoor.y]++;
									//if(patch==40 && subpatch==0)// && (patch==0 || patch==5 || patch==20))
										//std::cout << s*3+2 << ' ' << nG.x << ' ' << nG.y << ' ' << nG.z << std::endl;
										//std::cout << patch << ' ' << sG.x << ' ' << sG.y << ' ' << sG.z << std::endl;
										//std::cout << patchCoor.x+3*patchCoor.y << ' ' << sG.x << ' ' << sG.y << ' ' << sG.z << std::endl;
								}
							}
							
							bool badPatch=false;
							
							for(int a=0;a<3;a++)
								for(int b=0;b<3;b++)
									if(nZPatch.s[a].s[b]==0)
										badPatch=true;
										//std::cerr << "Bad Patch! " << zPatch.s[a].s[b] << ' ' << a << ' ' << b << std::endl;
							
							//average of every patch
							
							
							if(nZPatch.s[0].s[0]!=0)
								zPatch.s[0].s[0]/=nZPatch.s[0].s[0];
							if(nZPatch.s[0].s[1]!=0)
								zPatch.s[0].s[1]/=nZPatch.s[0].s[1];
							if(nZPatch.s[0].s[2]!=0)
								zPatch.s[0].s[2]/=nZPatch.s[0].s[2];
							if(nZPatch.s[1].s[0]!=0)
								zPatch.s[1].s[0]/=nZPatch.s[1].s[0];
							if(nZPatch.s[1].s[1]!=0)
								zPatch.s[1].s[1]/=(nZPatch.s[1].s[1]+1.0);//aP is here, but contributes 0
							if(nZPatch.s[1].s[2]!=0)
								zPatch.s[1].s[2]/=nZPatch.s[1].s[2];
							if(nZPatch.s[2].s[0]!=0)
								zPatch.s[2].s[0]/=nZPatch.s[2].s[0];
							if(nZPatch.s[2].s[1]!=0)
								zPatch.s[2].s[1]/=nZPatch.s[2].s[1];
							if(nZPatch.s[2].s[2]!=0)
								zPatch.s[2].s[2]/=nZPatch.s[2].s[2];
							
							zPatch.s[0].s[1]-=zPatch.s[1].s[1];
							zPatch.s[1].s[0]-=zPatch.s[1].s[1];
							zPatch.s[1].s[2]-=zPatch.s[1].s[1];
							zPatch.s[2].s[1]-=zPatch.s[1].s[1];
							
							zPatch.s[0].s[0]-=zPatch.s[1].s[1];
							zPatch.s[0].s[2]-=zPatch.s[1].s[1];
							zPatch.s[2].s[0]-=zPatch.s[1].s[1];
							zPatch.s[2].s[2]-=zPatch.s[1].s[1];
							
							zPatch.s[1].s[1]=0.0;
							
							//double orthFactor=0.816633/(patchSize);
							//double crossFactor=0.816633/sqrt(patchSize*2.0);
							
							//zPatch.s[0].s[1]*=orthFactor;
							//zPatch.s[1].s[0]*=orthFactor;
							//zPatch.s[1].s[2]*=orthFactor;
							//zPatch.s[2].s[1]*=orthFactor;
							
							
							//zPatch.s[0].s[0]*=crossFactor;
							//zPatch.s[0].s[2]*=crossFactor;
							//zPatch.s[2].s[0]*=crossFactor;
							//zPatch.s[2].s[2]*=crossFactor;
							
							if(!badPatch)
							{
								//double myGradSquared=0.66689*pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
								//		zPatch.s[2].s[1]+zPatch.s[0].s[1])/2.0+
								//		(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
								//		zPatch.s[2].s[2]+zPatch.s[0].s[0])/4.0-
								//		3.0*zPatch.s[1].s[1])/(patchSize*patchSize),2.0);//average 0.66689 nm^2 per lipid
								
								double patchArea=patchSize*patchSize;
								double myGradSquared=pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
										zPatch.s[2].s[1]+zPatch.s[0].s[1])*2.0/3.0+
										(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
										zPatch.s[2].s[2]+zPatch.s[0].s[0])/6.0-
										(10.0/3.0)*zPatch.s[1].s[1])/(patchArea),2.0);//average 0.66689 nm^2 per lipid
								
								//gradSquared+=myGradSquared;
		//std::cout << patchSize << ' ' << myGradSquared << '\n';
								if(myGradSquared<minGradSquare)
								{
									minGradSquare=myGradSquared;
									minHistIndex=static_cast<int>(patchSize/patchInc);
								}
								
								//int histIndex=static_cast<int>(myGradSquared/dGrad);
								//int histIndex=static_cast<int>(0.66689*myGradSquared/dGrad);
								//int histIndex=static_cast<int>(patchSize/patchInc);
								/*while(histogram.size()<=histIndex)
								{
									threeVector<double> hist;
									//hist.x=static_cast<double>(histogram.size())*dGrad+0.5*dGrad;
									hist.x=static_cast<double>(histogram.size())*patchInc+0.5*patchInc;
									hist.y=0;
									hist.z=0;
									histogram.push_back(hist);
									gradSquareValues.push_back(std::vector<double>());
								}
								*/
								//histogram[histIndex].y+=0.66689*myGradSquared;
								//histogram[histIndex].y+=myGradSquared;
								//histogram[histIndex].z+=1.0;
								
							//	threeVector<double> d;
							//	d.x=aP.x-centerOfMass[iSurface].x;
							//	d.y=aP.y-centerOfMass[iSurface].y;
							//	d.z=aP.z-centerOfMass[iSurface].z;
							//	double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
							//	int radIndex=static_cast<double>(dr/dRad);
								
							//	while(radDist.size()<=radIndex)
							//	{
							//		twoVector<double> hist;
							//		hist.x=static_cast<double>(radDist.size())*dRad+0.5*dRad;
							//		hist.y=0;
							//		radDist.push_back(hist);
							//	}
								//radDist[radIndex].y+=1.0;
								//patchSize=(1/3)*actualPatchSize;
								//int histIndex=static_cast<int>(patchSize/patchInc);
								//gradSquareValues[histIndex].push_back(myGradSquared);
								//histogram[histIndex].y+=myGradSquared;
								//histogram[histIndex].z+=1.0;
								//nFound++;
								myNormal++;
							}
						}
		//std::cout << '\n';
						histogram[minHistIndex].y+=minGradSquare;
						histogram[minHistIndex].z+=1.0;
						gradSquared+=minGradSquare;
						//throw 0;
						//subpatch++;
						//covariant matrix values A,B,C
						// [A B]
						// [B C]
						//LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
						
						//get the eigenvectors, principal curvature directions (lapack) 
						
						//get the curvature
						
						//patch++;
						//std::cerr << time << ' ' << patch << ' '  << stack.size() << ' ' << filt.size() << ' ' << currentHash << std::endl;
						//if(patch==40)
						//throw 0;
					}
					//throw 0;
				}
				
			}
			iSurface++;
			//throw 0;
			if(filt.size()>5)
			{
				std::cout << ' ' << gradSquared << ' ' << surface.size();
				
				std::string buf;
				//buf="patchProfile_";
				buf="patchMin_";
				buf+=argv[1];
				buf+=".dat";
				std::fstream dataFile;
				
				dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
				//auto gradSquareValue=gradSquareValues.begin();
				
				//for(auto &value:gradSquareValues)
				//		dataFile << value << std::endl;
				
				for(auto &P:histogram)
				{
					//P.x=(P.x+patchInc/2.0);//scaled to 3 patch size
					//P.y/=P.z;//average value
					//double rms=0;
					//for(auto &value:*gradSquareValue)
					//	rms+=pow(value-P.y,2.0);
					//rms=sqrt(rms/static_cast<double>(gradSquareValue->size()));
					//if(P.z!=0)
					//	dataFile << P.x << ' ' << P.y*P.z << ' ' << rms << ' ' << gradSquareValue->size() << '\n';
					if(P.z!=0)
						dataFile << P.x << ' ' << P.y << ' ' << P.z << '\n';
					//gradSquareValue++;
				}
				
				dataFile << std::endl;
				
				dataFile.close();
				/*
				buf.clear();
				buf="patchProfile_";
				//buf="radDist_";
				buf+=argv[1];
				buf+=".dat";
				
				dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
				//auto gradSquareValue=gradSquareValues.begin();
				for(auto &P:radDist)
				{
						dataFile << P.x << ' ' << P.y << '\n';
				}
				dataFile << std::endl;
				
				dataFile.close();
				*/
			}
		}
		if(vFilt.size()>1)
		{
			std::cout << std::endl;
			throw 0;
		}
		//throw 0;
		/*
		if(triangles.size()>0 && dry==0 && frame==0)
		{
			
			//std::cout << gradSquared << '\t';// << hSum.x << '\t' << hSum.y << std::endl;
			//std::cerr << "Surface " << s << " gradSquared=" << gradSquared << " vectSums=" << aSum.x << ',' << aSum.y << ',' << bSum.x << ',' << bSum.y;
			//for(int neigh=0;neigh<nNeigh.size();neigh++)
			//{
			//	if(nNeigh[neigh].size()>0)
			//	{
			surfaceNumber=0;
			//for(auto& triangles:sTriangles)
			{
				std::string buf;
				buf="triangleMesh_";
				std::string whatever;
				std::stringstream parseNumber;
				parseNumber << surfaceNumber;
				parseNumber >> whatever;
				
				buf+=whatever;
				buf+="_";
				buf+=argv[1];
				buf+=".stl";
				std::fstream dataFile;
				dataFile.open(buf.c_str(), std::ios::out);// | std::ios::app);
				
				dataFile << "solid " << argv[1] << '\n';
				//dump the mesh to triangle file
				//for(int i=0;i<triangles.size();i++)
				
				//if(!delTriangle[i])
				for(auto& triangle:triangles)
				{
					//const auto& triangle=triangles[i];
					dataFile << "facet normal " << 1.0 << ' ' << 1.0 << ' ' << 1.0 << '\n';
					dataFile << "\touter loop " << '\n';
					for(int j=0;j<3;j++)
					{
						//auto k=triangles[i].s[j];
						//if(k!=-1)
						{
							dataFile << "\t\tvertex " <<
							triangle.s[j]->second.x << ' ' <<
							triangle.s[j]->second.y << ' ' <<
							triangle.s[j]->second.z << '\n';
						}
						//else
						//{
						//	dataFile << "\t\tnotvertex " <<
						//		k << ' ' << '\n';
						//k << ' ' <<
						//k << std::endl;
						//}
					}
					dataFile << "\tendloop" << '\n';
					dataFile << "endfacet" << '\n';
				}
				
				dataFile << "endsolid " << argv[1] << '\n';
				dataFile.close();
				surfaceNumber++;
				//}
			}
			surfaceNumber=0;
			for(auto &filt:vFilt)
			{
				std::string buf;
				buf="vertex_";
				std::string whatever;
				std::stringstream parseNumber;
				parseNumber << surfaceNumber;
				parseNumber >> whatever;
				
				buf+=whatever;
				buf+="_";
				buf+=argv[1];
				buf+=".xyz";
				std::fstream dataFile;
				dataFile.open(buf.c_str(), std::ios::out);// | std::ios::app);
				
				dataFile << filt.size() << "\nvertex\n";
				
				for(auto &vertex:filt)
					dataFile << 1 << ' ' << vertex.second.x << ' ' << vertex.second.y << ' ' << vertex.second.z << '\n';
				dataFile.close();
				surfaceNumber++;
			}
		}
		*/
		frame++;
	}
	return 0;
}
