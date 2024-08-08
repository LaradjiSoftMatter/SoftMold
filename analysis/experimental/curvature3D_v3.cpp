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
		std::cerr << argv[0] << " name surfSpacing minAngle vertShellRatio vertSpacing dry type1 cType1 ...\n";
		std::cerr << "\tdry is for no output, if the surface spacing needs to be optimized. 0=false 1=true\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double surfSpacing, minAngle,vertSpacing,vertShellRatio;
	int dry;
	cmdArg >> surfSpacing >> minAngle >> vertShellRatio >> vertSpacing >> dry;
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
		if(cTypes[i]<System.readNMolecules())
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
	bool dumpedMesh=false;
	while(xyzFile.load())
	{
		//subsample pIndex to pFilt on a grid, weighted by center of mass
		//std::map<int, position<double> >::iterator->first is an integer for the hashed grid position
		//std::map<int, position<double> >::iterator->second is the weighted central position
		//position pos=pFilt[gridId]; would get you the actual position
		std::map<int, position<double> > pFilt;
		
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
			{
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				//std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
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
		
		
		
		std::map<int,cellType > cells,etoe;
		
		threeVector<double> zero=0;
		
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
			}
			else
			{
				//initialize operation
				cellType buf;
				buf.push_back(p[i]);
				
				cells.insert(it, std::map<int,cellType >::value_type(hash,buf));
				
				buf[0]=p[pIndex[j].y];
				etoe[hash]=buf;
			}
		}
		
		std::cerr << "using " << pIndex.size() << " points for surfaces" << std::endl;
		
		//subsample our pIndex
		//c++11 syntax reduces this:
		//std::map<int,cellType >::iterator it
		for(auto& it: cells)
		{
			int gridPos=it.first;
			
			position<double> pAvg;
			pAvg.x=0;
			pAvg.y=0;
			pAvg.z=0;
			
			//std::vector<position>::iterator pI
			for(auto& pI: it.second)
			{
				//int i=pIndex[pI];
				pAvg.x+=pI.x;
				pAvg.y+=pI.y;
				pAvg.z+=pI.z;
			}
			pAvg.x/=static_cast<double>(it.second.size());
			pAvg.y/=static_cast<double>(it.second.size());
			pAvg.z/=static_cast<double>(it.second.size());
			pAvg.type=-1;
			pFilt[gridPos]=pAvg;
		}
		
		std::cerr << "using " << pFilt.size() << " points for surface discrimination" << std::endl;
		
		//std::vector< vector<position> >
		std::vector< cellType > surfaces, surfacesETOE;
		
		//obtaining surfaces
		for(auto pFIT=pFilt.begin();pFIT!=pFilt.end();pFIT++)
		{
			//unaccounted surface
			//if(sFlag[pFIT.first]==-1)
			if((*pFIT).second.type==-1)
			{
				(*pFIT).second.type=surfaces.size();
				//std::deque<int> stack;
				//std::vector<int> stack;
				std::deque<std::map<int, position<double> >::iterator > stack;
				
				stack.push_back(pFIT);
				
				//std::vector<std::map<int, position<double> >::iterator> surface;
				cellType surface,surfaceETOE;
				
				while(stack.size()>0)
				{
					auto currentPoint=stack.back();
					
					int currentHash=currentPoint->first;
					
					position<double> aP=currentPoint->second;
					
					threeVector<int> iCoor=unhash(currentHash,surfN);
					
					//threeVector<double> normalSum=0;
					//int nNormals=0;
					//auto otherEnd=etoe[currentHash].begin();
					
					for(auto& cellPoint:cells[currentHash])
					{
					//	threeVector<double> normal;
					//	normal.x=otherEnd->x-cellPoint.x;
					//	normal.y=otherEnd->y-cellPoint.y;
					//	normal.z=otherEnd->z-cellPoint.z;
					//	normal.x-=(normal.x>size.x/2.0)?size.x:0;
					//	normal.y-=(normal.y>size.y/2.0)?size.y:0;
					//	normal.z-=(normal.z>size.z/2.0)?size.z:0;
					//	normal.x+=(normal.x<-size.x/2.0)?size.x:0;
					//	normal.y+=(normal.y<-size.y/2.0)?size.y:0;
					//	normal.z+=(normal.z<-size.z/2.0)?size.z:0;
					//	
					//	normalSum+=normal;
					//	nNormals++;
						
					//	otherEnd++;
						
						surface.push_back(cellPoint);
						
					}
					
					//normalSum/=static_cast<double>(nNormals);
					//normalSum=unitVector(normalSum);
					
					
					for(auto& cellPoint:etoe[currentHash])
						surfaceETOE.push_back(cellPoint);
					
					stack.pop_back();
					
					//surface search
					for(int k=0;k<27;k++)
					{
						threeVector<int> nCoor=iCoor;
						nCoor.x+=k%3-1;
						nCoor.y+=int(k/3)%3-1;
						nCoor.z+=int(k/9)-1;
						int neighborHash=hash(nCoor,surfN);
						auto it=pFilt.find(neighborHash);
						if(it!=pFilt.end() && it->first==neighborHash && currentHash!=neighborHash)
						{
							if(it->second.type==-1)
							{
								/*
								threeVector<double> otherSum=0;
								nNormals=0;
								otherEnd=etoe[neighborHash].begin();
								
								for(auto& cellPoint:cells[neighborHash])
								{
									threeVector<double> normal;
									normal.x=otherEnd->x-cellPoint.x;
									normal.y=otherEnd->y-cellPoint.y;
									normal.z=otherEnd->z-cellPoint.z;
									normal.x-=(normal.x>size.x/2.0)?size.x:0;
									normal.y-=(normal.y>size.y/2.0)?size.y:0;
									normal.z-=(normal.z>size.z/2.0)?size.z:0;
									normal.x+=(normal.x<-size.x/2.0)?size.x:0;
									normal.y+=(normal.y<-size.y/2.0)?size.y:0;
									normal.z+=(normal.z<-size.z/2.0)?size.z:0;
									
									otherSum+=normal;
									nNormals++;
									otherEnd++;
								}
								otherSum/=static_cast<double>(nNormals);
								otherSum=unitVector(otherSum);
								*/
								
								position<double> neP=it->second;
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
								
								//The multiplying factor maintains surface identity
								//sqrt(1.0)=1     (linear)
								//sqrt(2.0)=1.414 (surface)
								//sqrt(3.0)=1.73  (volume)
								//if(d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
								d.x=k%3-1;
								d.y=int(k/3)%3-1;
								d.z=int(k/9)-1;
								
								if(d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)// && dotProduct(normalSum,otherSum)>sqrt(2.0)/2.0)
								{
									it->second.type=surfaces.size();
									//stack.push_back(neighborHash);
									stack.push_back(it);
								}
							}
						}
					}
				}
				surfaces.push_back(surface);
				surfacesETOE.push_back(surfaceETOE);
				//nSurfaces++;
			}
		}
		std::cerr << "Found " << surfaces.size() << " surfaces." << std::endl;
		
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
		
		auto sETOE=surfacesETOE.begin();
		
		for(auto& surface:surfaces)
		{
			if(surface.size()>10)
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
				itETOE++;
			}
			//nSurface++;
			surfCells.push_back(sCells);
			surfNormals.push_back(normals);
			sETOE++;
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
		
		std::cout << time;
		
		//try every surface to get surface normals
		for(int s=0;s<vFilt.size();s++)
		{
			double gradSquared=0;
			#define vertexPoint 0
			#define edgePoint 1
			#define ignorePoint 2
			auto& filt=vFilt[s];//for nearest neighbor traversal
			auto& sCells=surfCells[s];//for vertex candidates
			auto& normals=surfNormals[s];
			//auto& filtNormal=vFiltNormals[s];
			
			//unaccounted surface
			for(auto element=filt.begin();element!=filt.end() && filt.size()>5;element++)
			{
				//unaccounted surface element?
				if(element->second.type==unusedPoint)
				{
					/*
					//this will use filt to align the normals
					//we will assume the first normal on filt is aligned if it is positive Z
					threeVector<double> normalSum=0;
					double nNormals=0;
					std::vector<position<double> > firstBcP;
					
					for(int k=0;k<27;k++)
					{
						threeVector<int> nCoor=unhash(element->first,vertN);
						nCoor.x+=k%3-1;
						nCoor.y+=int(k/3)%3-1;
						nCoor.z+=int(k/9)-1;
						int neighborHash=hash(nCoor,vertN);
						
						position<double> aP=element->second;
						
						
						auto it=filt.find(neighborHash);
						if(it!=filt.end() && it->first==neighborHash && neighborHash!=element->first)
							firstBcP.push_back(it->second);
					}*/
					
					threeVector<double> oneZ=0;
					oneZ.z=1.0;
					/*
					for(auto bP=firstBcP.begin();bP!=firstBcP.end();bP++)
					{
						position<double> aP=element->second;
						threeVector<double> dab;
						dab.x=aP.x-bP->x;
						dab.y=aP.y-bP->y;
						dab.z=aP.z-bP->z;
						
						while(dab.x>size.x/2.0)dab.x-=size.x;
						while(dab.y>size.y/2.0)dab.y-=size.y;
						while(dab.z>size.z/2.0)dab.z-=size.z;
						while(dab.x<-size.x/2.0)dab.x+=size.x;
						while(dab.y<-size.y/2.0)dab.y+=size.y;
						while(dab.z<-size.z/2.0)dab.z+=size.z;
						
						for(auto cP=bP+1;cP!=firstBcP.end();cP++)
						{
							threeVector<double> dbc;
							dbc.x=bP->x-cP->x;
							dbc.y=bP->y-cP->y;
							dbc.z=bP->z-cP->z;
							
							while(dbc.x>size.x/2.0)dbc.x-=size.x;
							while(dbc.y>size.y/2.0)dbc.y-=size.y;
							while(dbc.z>size.z/2.0)dbc.z-=size.z;
							while(dbc.x<-size.x/2.0)dbc.x+=size.x;
							while(dbc.y<-size.y/2.0)dbc.y+=size.y;
							while(dbc.z<-size.z/2.0)dbc.z+=size.z;
							
							threeVector<double> thisNormal=unitVector(crossProduct(dab,dbc));
							if(dotProduct(thisNormal,oneZ)<0)
							{
								thisNormal.x=-thisNormal.x;
								thisNormal.y=-thisNormal.y;
								thisNormal.z=-thisNormal.z;
							}
							normalSum+=thisNormal;
							nNormals++;
						}
					}
					
					normalSum/=nNormals;
					
					filtNormal[element->first]=normalSum;
					*/
					element->second.type=usedPoint;
					std::deque<filtType::iterator > stack;
					
					stack.push_back(element);
					int patch=0;	
					while(stack.size()>0)
					{
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
						
						auto currentPoint=stack.back();
						
						int currentHash=currentPoint->first;
						
						threeVector<int> iCoor=unhash(currentHash,vertN);
						//std::cout << iCoor.x << ' ' << iCoor.y << ' ' << iCoor.z << std::endl;
						//std::cin.get();
						
						currentPoint->second.type=usedPoint;
							
						stack.pop_back();
						/*
						//get the local alignment
						std::vector<threeVector<double> > nList;
						for(int k=0;k<27;k++)
						{
							threeVector<int> nCoor=iCoor;
							nCoor.x+=k%3-1;
							nCoor.y+=int(k/3)%3-1;
							nCoor.z+=int(k/9)-1;
							int neighborHash=hash(nCoor,vertN);
							
							auto it=filtNormal.lower_bound(neighborHash);
							if(it!=filtNormal.end() && it->first==neighborHash)
								if(it->second.x!=0.0 && it->second.y!=0 && it->second.z!=0)
									nList.push_back(it->second);
						}
						
						threeVector<double> localNormal=0;
						for(auto &nL: nList)
							localNormal+=nL;
						localNormal/=static_cast<double>(nList.size());
						
						//quick checks, these can be turned off
						if(nList.size()<1)
						{
							std::cerr << "No local alignment found!" << std::endl;
							//throw 0;
						}
						else
						{
							for(auto aL=nList.begin();aL!=nList.end();aL++)
							{
								for(auto bL=nList.begin();bL!=nList.end();bL++)
								{
									if(dotProduct(*aL,*bL)<0)
									{
										std::cerr << "Local alignment doesn't agree!\n Try smaller normalSpacing!" << std::endl;
										std::cerr << "a: " << (*aL).x << ' ' << (*aL).y << ' ' << (*aL).z << std::endl;
										std::cerr << "b: " << (*bL).x << ' ' << (*bL).y << ' ' << (*bL).z << std::endl;
										threeVector<double> aCrossB=crossProduct(*aL,*bL);
										std::cerr << "axb: " << aCrossB.x << ' ' << aCrossB.y << ' ' << aCrossB.z << std::endl;
										throw 0;
									}
								}
							}
						}
						*/
						//keep up with the local alignment
						//normalSum=oneZ;
						//nNormals=0;
						std::vector<position<double> > bcP;
						
						for(int k=0;k<27;k++)
						{
							threeVector<int> nCoor=iCoor;
							nCoor.x+=k%3-1;
							nCoor.y+=int(k/3)%3-1;
							nCoor.z+=int(k/9)-1;
							int neighborHash=hash(nCoor,vertN);
							
							auto it=filt.lower_bound(neighborHash);
							if(it!=filt.end() && it->first==neighborHash && neighborHash!=currentHash)
							{
								bcP.push_back(it->second);
								if(it->second.type==unusedPoint)
								{
									stack.push_back(it);
									it->second.type=usedPoint;
								}
							}
						}
						auto& aGroup=sCells[currentHash];
						/*
						//auto& aGroup=sCells[currentHash];
						auto test=sCells.begin();
						if(test->second.size()<10)
							test++;
						auto& aGroup=test->second;
						threeVector<threeVector<double> > F;
						F.x.x=0;
						F.x.y=0;
						F.x.z=0;
						F.y.x=0;
						F.y.y=0;
						F.y.z=0;
						F.z.x=0;
						F.z.y=0;
						F.z.z=0;
						//      [x^2 xy  x]          [xz]
						//F=sum([xy  y^2 y]), zV=sum([yz])
						//      [x   y   n]          [z ]
						//        [a]
						//F^-1*zV=[b], where z-ax-by-c=0, normal=[[a/c][b/c][1/c]]
						//        [c]
						threeVector<double> zV;
						zV.x=0;
						zV.y=0;
						zV.z=0;
						for(auto &aP:aGroup)
						{
							F.x.x+=aP.x*aP.x;//f
							F.x.y+=aP.x*aP.y;//g
							F.x.z+=aP.x;//h
							F.y.y+=aP.y*aP.y;//r
							F.y.z+=aP.y;//s
							F.z.z+=1.0;//v
							zV.x+=aP.x*aP.z;//l
							zV.y+=aP.y*aP.z;//m
							zV.z+=aP.z;//n
							std::cout << "1 " << aP.x << ' ' << aP.y << ' ' << aP.z << std::endl;
						}
						
						//localNormal.x=a/c,localNormal.y=b/c, localNormal.z=1/c
						//Something is missing, these should all be something like c1*c2-c3*c4
						threeVector<double> localNormal;
						double msnr=zV.y*F.y.z-zV.z*F.y.y;
						double fshg=F.x.x*F.y.z-F.x.z*F.x.y;
						double hgs=F.x.z-F.y.z*F.y.z;
						double lsng=zV.x*F.y.z-zV.z*F.x.y;
						double ssvr=F.y.z*F.y.z-F.z.z*F.y.y;
						double vhs=F.z.z-F.x.z*F.y.z;
						double fnlh=F.x.x*zV.z-zV.x*F.x.z;
						double hfv=F.x.z-F.x.x*F.z.z;
						double fsgh=F.x.x*F.y.z-F.x.y*F.x.z;
						
						localNormal.z=(fsh*ssv-vhs*hgs)/(msn*fsh+hgs*lsn);//1/c
						localNormal.x=(lsn+vhs/localNormal.z)/fsgh;
						localNormal.y=(fnl+hfv/localNormal.z)/fsgh;
						
						auto aPtemp=aGroup[0];
						std::cout << "3 " << aPtemp.x+localNormal.x << ' ' << aPtemp.y+localNormal.y << ' ' << aPtemp.z+localNormal.z << std::endl;
						
						//testing the fit
						for(auto &aP:aGroup)
						{
							std::cout << "2 " << aP.x << ' ' << aP.y << ' ' <<
								localNormal.x*aP.x+localNormal.y*aP.y+1.0/localNormal.z << std::endl;
						}
						
						throw 0;
						/*
						/*
						for(auto aP=aGroup.begin();aP!=aGroup.end();aP++)
						for(auto bP=aP+1;bP!=aGroup.end();bP++)
						{
							threeVector<double> dab;
							dab.x=aP->x-bP->x;
							dab.y=aP->y-bP->y;
							dab.z=aP->z-bP->z;
							
							while(dab.x>size.x/2.0)dab.x-=size.x;
							while(dab.y>size.y/2.0)dab.y-=size.y;
							while(dab.z>size.z/2.0)dab.z-=size.z;
							while(dab.x<-size.x/2.0)dab.x+=size.x;
							while(dab.y<-size.y/2.0)dab.y+=size.y;
							while(dab.z<-size.z/2.0)dab.z+=size.z;
							
							for(auto cP=bP+1;cP!=aGroup.end();cP++)
							{
								threeVector<double> dbc;
								dbc.x=bP->x-cP->x;
								dbc.y=bP->y-cP->y;
								dbc.z=bP->z-cP->z;
								
								while(dbc.x>size.x/2.0)dbc.x-=size.x;
								while(dbc.y>size.y/2.0)dbc.y-=size.y;
								while(dbc.z>size.z/2.0)dbc.z-=size.z;
								while(dbc.x<-size.x/2.0)dbc.x+=size.x;
								while(dbc.y<-size.y/2.0)dbc.y+=size.y;
								while(dbc.z<-size.z/2.0)dbc.z+=size.z;
								
								threeVector<double> thisNormal=unitVector(crossProduct(dab,dbc));
								if(dotProduct(thisNormal,normalSum)<0)
								{
									thisNormal.x=-thisNormal.x;
									thisNormal.y=-thisNormal.y;
									thisNormal.z=-thisNormal.z;
								}
								
								normalSum+=thisNormal;
								normalSum=unitVector(normalSum);
								nNormals++;
							}
						}
						*/
						
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
						//we don't need to do everybodies normal
						/*
						for(auto aP=aGroup.begin();aP!=aGroup.end();aP++)
						{
							
							threeVector<double> aNormal=0;
							nNormals=0;
							
							for(auto bP=aP+1;bP!=aGroup.end();bP++)
							{
								threeVector<double> dab;
								dab.x=aP->x-bP->x;
								dab.y=aP->y-bP->y;
								dab.z=aP->z-bP->z;
								
								while(dab.x>size.x/2.0)dab.x-=size.x;
								while(dab.y>size.y/2.0)dab.y-=size.y;
								while(dab.z>size.z/2.0)dab.z-=size.z;
								while(dab.x<-size.x/2.0)dab.x+=size.x;
								while(dab.y<-size.y/2.0)dab.y+=size.y;
								while(dab.z<-size.z/2.0)dab.z+=size.z;
								
								for(auto cP=bP+1;cP!=aGroup.end();cP++)
								{
									threeVector<double> dbc;
									dbc.x=bP->x-cP->x;
									dbc.y=bP->y-cP->y;
									dbc.z=bP->z-cP->z;
									
									while(dbc.x>size.x/2.0)dbc.x-=size.x;
									while(dbc.y>size.y/2.0)dbc.y-=size.y;
									while(dbc.z>size.z/2.0)dbc.z-=size.z;
									while(dbc.x<-size.x/2.0)dbc.x+=size.x;
									while(dbc.y<-size.y/2.0)dbc.y+=size.y;
									while(dbc.z<-size.z/2.0)dbc.z+=size.z;
									
									threeVector<double> thisNormal=unitVector(crossProduct(dab,dbc));
									if(dotProduct(thisNormal,normalSum)<0)
									{
										thisNormal.x=-thisNormal.x;
										thisNormal.y=-thisNormal.y;
										thisNormal.z=-thisNormal.z;
									}
									aNormal+=thisNormal;
									nNormals++;
								}
							}
							aNormal/=nNormals;
							
							//std::cout << "Here? " << magnitude(aNormal) << ' ' << aNormal.x << ' ' << aNormal.y << ' ' << aNormal.z << std::endl;
							//std::cin.get();
							
							(*normal++)=aNormal;
						}
						*/
						//vector projections
						//std::vector<threeVector<double> > S(nearbyGroups.size(),zero);
						
						//together these are min(|Mu-K|)
						//curvatures between aP and S
						//std::vector<double> K(nearbyGroups.size(),0);
						//our minimization matrix, M * mu
						//std::vector<threeVector<double> > Mu(nearbyGroups.size(),zero);
						
						//reset this
						//normal=normals[currentHash].begin();//iterator
						
						//auto normal=normals[currentHash].begin();
						
						//identity matrix, ROW_MAJOR
						//threeVector<threeVector<double> > I;
						//I.x=0;
						//I.y=0;
						//I.z=0;
						//I.x.x=1.0;
						//I.y.y=1.0;
						//I.z.z=1.0;
						
						auto normalGroup=normals[currentHash];
						
						//std::cerr << nearbyGroups.size() << std::endl;
						int subpatch=0;
						//do this for every particle in aGroup
						for(auto& aP:aGroup)
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
								if(d.x*d.x+d.y*d.y+d.z*d.z<vertSpacingSqr*0.75*0.75)
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
								if(d.x*d.x+d.y*d.y+d.z*d.z<vertSpacingSqr*0.75*0.75)
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
							
							
							threeVector<double> aPuT;
							aPuT.x=aPu.x*M.x.x+aPu.y*M.x.y+aPu.z*M.x.z;
							aPuT.y=aPu.x*M.y.x+aPu.y*M.y.y+aPu.z*M.y.z;
							aPuT.z=aPu.x*M.z.x+aPu.y*M.z.y+aPu.z*M.z.z;
							if(dotProduct(aPuT,oneZ)<sqrt(2.0)/2.0)
							{
								std::cerr << "It came out here: " << aPuT.x << ' ' << aPuT.y << ' ' << aPuT.z << std::endl;
								std::cerr << "With magnitude: " << magnitude(aPuT) << std::endl;
								std::cerr << "from here: " << aPu.x << ' ' << aPu.y << ' ' << aPu.z << std::endl;
							}
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
								if(sG.x>-0.75*vertSpacing && sG.x<0.75*vertSpacing)
									patchCoor.x=static_cast<int>((sG.x+vertSpacing*0.75)/(vertSpacing*0.5));
								if(sG.y>-0.75*vertSpacing && sG.y<0.75*vertSpacing)
									patchCoor.y=static_cast<int>((sG.y+vertSpacing*0.75)/(vertSpacing*0.5));
								//if(sG.z>0 && sG.z<1.5*vertSpacing)
								//	patchCoor.z=static_cast<int>(sG.z/vertSpacing);
								
								if(qG.x*qG.x+qG.y*qG.y+qG.z*qG.z<1.5*1.5*vertSpacingSqr)
								if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
								{
									zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
									nZPatch.s[patchCoor.x].s[patchCoor.y]++;
									//if(subpatch==0)// && (patch==0 || patch==5 || patch==20))
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
								
								//sG=qG;
								
								threeVector<int> patchCoor=-1.0;
								if(sG.x>-0.75*vertSpacing && sG.x<0.75*vertSpacing)
									patchCoor.x=static_cast<int>((sG.x+vertSpacing*0.75)/(vertSpacing*0.5));
								if(sG.y>-0.75*vertSpacing && sG.y<0.75*vertSpacing)
									patchCoor.y=static_cast<int>((sG.y+vertSpacing*0.75)/(vertSpacing*0.5));
								//if(sG.z>0 && sG.z<1.5*vertSpacing)
								//	patchCoor.z=static_cast<int>(sG.z/vertSpacing);
								
								if(qG.x*qG.x+qG.y*qG.y+qG.z*qG.z<1.5*1.5*vertSpacingSqr)
								if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
								{
									zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
									nZPatch.s[patchCoor.x].s[patchCoor.y]++;
									//if(subpatch==0)// && (patch==0 || patch==5 || patch==20))
									//	std::cout << s*3+2 << ' ' << nG.x << ' ' << nG.y << ' ' << nG.z << std::endl;
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
							
							zPatch.s[0].s[1]*=orthFactor;
							zPatch.s[1].s[0]*=orthFactor;
							zPatch.s[1].s[2]*=orthFactor;
							zPatch.s[2].s[1]*=orthFactor;
							
							
							zPatch.s[0].s[0]*=crossFactor;
							zPatch.s[0].s[2]*=crossFactor;
							zPatch.s[2].s[0]*=crossFactor;
							zPatch.s[2].s[2]*=crossFactor;
							
							if(!badPatch)
							gradSquared+=pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
										zPatch.s[2].s[1]+zPatch.s[0].s[1])/2.0+
										(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
										zPatch.s[2].s[2]+zPatch.s[0].s[0])/4.0-
										3.0*zPatch.s[1].s[1]),2.0)/(0.6431);//temporary 0.6431 nm^2 per lipid
							
							//throw 0;
							//subpatch++;
							//covariant matrix values A,B,C
							// [A B]
							// [B C]
							//LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
							
							//get the eigenvectors, principal curvature directions (lapack) 
							
							//get the curvature
						}
						//patch++;
						//std::cerr << time << ' ' << patch << ' '  << stack.size() << ' ' << filt.size() << ' ' << currentHash << std::endl;
						//if(patch==40)
						//	throw 0;
					}
					//throw 0;
				}
			}
			if(filt.size()>5)
				std::cout << ' ' << gradSquared;
		}
		
		std::cout << std::endl;
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
