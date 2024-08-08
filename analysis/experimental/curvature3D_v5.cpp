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
bool compByT(fourVector<double> a, fourVector<double> b)
{
	//if(a.x<b.x)
	//	return true;
	//return false;
	return a.t<b.t;
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
//	while(c.x>=nFactor.x)c.x-=nFactor.x;
//	while(c.y>=nFactor.y)c.y-=nFactor.y;
//	while(c.z>=nFactor.z)c.z-=nFactor.z;
//	while(c.x<0)c.x+=nFactor.x;
//	while(c.y<0)c.y+=nFactor.y;
//	while(c.z<0)c.z+=nFactor.z;
	c.x=(c.x+nFactor.x)%nFactor.x;
	c.y=(c.y+nFactor.y)%nFactor.y;
	c.z=(c.z+nFactor.z)%nFactor.z;

	
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

int main(int argc, char **argv)
{
	if(argc<8)
	{
		std::cerr << argv[0] << " name surfSpacing noTraverseRatio frame vertSpacing nVertSpacing type1 cType1 ..." << std::endl;
		std::cerr << "\tframe performs the measure at a frame (-1 for all frames)." << std::endl;
		std::cerr << "\tnVertSpacing indicates how many measures from [1,vertSpacing] will be performed." << std::endl;
		std::cerr << "\tif nVertSpacing is 1 or 0, the default is 1 measure at vertSpacing." << std::endl;
		std::cerr << "\tif nVertSpacing is -1, it will dump a map in (N/Nt, N/V ,C^2) space in n^2*log(N)^2 time!" << std::endl;
		std::cerr << "\tPossible files:" << std::endl;
		std::cerr << "\t\tpatchProfile_S#_name.dat" << std::endl;
		std::cerr << "\t\tgradSquare_S#_name.dat" << std::endl;
		std::cerr << "\t\tcMap_S#_name.xyz" << std::endl;
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double surfSpacing, noTraverseRatio,vertSpacing;
	int frame, nVertSpacing;
	cmdArg >> surfSpacing >> noTraverseRatio >> frame >> vertSpacing >> nVertSpacing;
	double surfSpacingSqr=surfSpacing*surfSpacing;
	double vertSpacingSqr=vertSpacing*vertSpacing;
	
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
						if(pIndex[k].x>start && pIndex[k].x<start+nChains*length)
						{
							int chain=(pIndex[k].x-start)/length;
							int seg=(pIndex[k].x-start)%length;
							if(seg==0)//its at the 0 end, use the length-1 end
								pIndex[k].y=(length-1)+chain*length;
							else if(seg==length-1)//its at the length-1 end, use the 0 end
								pIndex[k].y=chain*length;
							else//its in the middle, just use the 0 end
								pIndex[k].y=chain*length;
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
	int currentFrame=0;
	double avgNeighbors=0;
	
	while(xyzFile.load())
	{
		time+=System.readStoreInterval();
		std::cerr << time << '\t' << currentFrame << std::endl;
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
		
		std::map<int,cellType > cells, etoe;
		
		threeVector<double> zero=0;
		
		#define unusedFlag -1
		#define noTraverseFlag -2
		#define traverseFlag 0
		
		for(int j=0;j<pIndex.size() && (frame==currentFrame || frame==-1 || currentFrame%20==0);j++)
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
		
		if(currentFrame%20==0)
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
			std::cerr << "Density: " << avgNeighbors/(M_PI*surfSpacingSqr) << std::endl;
			std::cerr << "Using " << avgNeighbors*noTraverseRatio << " neighbors for no traverse condition!" << std::endl;
		}
		
		std::cerr << "using " << pIndex.size() << " points for surface discrimination" << std::endl;
		
		std::vector< std::vector< position<double> > > surfaces, surfacesETOE;
		
		
		//obtaining surfaces
		for(auto pFIT=pFilt.begin();pFIT!=pFilt.end() && (frame==currentFrame || frame==-1);pFIT++)//a cell
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
		
		double vertInc=2.0*vertSpacing;
		if(nVertSpacing>0)
			vertInc=(vertSpacing-1.0)/(static_cast<double>(nVertSpacing));
		
		for(double cVert=vertSpacing;cVert>1.0 && (frame==currentFrame || frame==-1);cVert-=vertInc)
		{
			//set up a quick search using the minimum surface spacing
			threeVector<int> vertN;
			vertN.x=static_cast<int>(size.x/(cVert));
			vertN.y=static_cast<int>(size.y/(cVert));
			vertN.z=static_cast<int>(size.z/(cVert));
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
					normal.x+=(normal.x<-size.x/2.0)?size.x:0;
					normal.y+=(normal.y<-size.y/2.0)?size.y:0;
					normal.z+=(normal.z<-size.z/2.0)?size.z:0;
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
					int myHash=hash(c,vertN);//c.x+c.y*vertN.x+c.z*vertN.x*vertN.y;
					
					std::map<int,cellType >::iterator it=sCells.lower_bound(myHash);
					if(it!=sCells.end() && it->first==myHash)//!(cells.key_comp()(it->first, hash)))
					{
						//push operation
						it->second.push_back(pI);
						
						normals[myHash].push_back(normal);
					}
					else
					{
						//initialize operation
						cellType buf;
						buf.push_back(pI);
						//cells[hash]=buf;
						sCells.insert(it, std::map<int,cellType >::value_type(myHash,buf));
						
						
						normalType normalBuf;
						normalBuf.push_back(normal);
						normals[myHash]=normalBuf;
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
			typedef std::map<int,int>  filtType;
			std::vector<filtType> vFilt;
			
			//typedef std::map<int,threeVector<double> > filtNormalType;
			//std::vector<filtNormalType> vFiltNormals;
			
			#define unusedPoint -1
			#define usedPoint 0;
			for(int s=0;s<surfCells.size();s++)
			{
				auto& sCells=surfCells[s];
				
				filtType filt;
				//filtNormalType filtNormal;
				//subsample our pIndex
				//c++11 syntax reduces this:
				//std::map<int,cellType >::iterator it
				for(auto& it: sCells)
				{
					int gridPos=it.first;
					
					filt[gridPos]=unusedPoint;
				}
				vFilt.push_back(filt);
			}
			int nBadPatches=0;
			
			double patchSize=cVert*sqrt(2.0)/3.0;
			
			//try every surface to get surface normals
			for(int s=0;s<vFilt.size() && vFilt.size()>1;s++)
			{
				std::vector<double> gradSquareValues;
				std::vector< fourVector<double> > gradSquarePos;
				double gradSquared=0;
				double laplacianInt=0;
				#define vertexPoint 0
				#define edgePoint 1
				#define ignorePoint 2
				auto& filt=vFilt[s];//for nearest neighbor traversal
				auto& sCells=surfCells[s];//for vertex candidates
				auto& normals=surfNormals[s];
				auto& surface=surfaces[s];
				
				double localRms=0;
				
				int patch=0;
				int nElements=0;
				int percentCompletion=0;
				//unaccounted surface
				for(auto element=filt.begin();element!=filt.end() && filt.size()>10;element++)
				{
					if(percentCompletion<=100.0*static_cast<double>(nElements++)/static_cast<double>(filt.size()))
					{
						std::cerr << percentCompletion << "%" << std::endl;
						percentCompletion+=10;
					}
					//unaccounted surface element?
					if(element->second==unusedPoint)
					{
						threeVector<double> oneZ=0;
						oneZ.z=1.0;
						
						element->second=usedPoint;
						int currentHash=element->first;
						
						threeVector<int> iCoor=unhash(currentHash,vertN);
						
						auto& aGroup=sCells[currentHash];
						
						//now align the clustered group
						
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
						
						int gradSquareValuesStart=gradSquareValues.size();
						double gradSquaredLow=gradSquared;
						
						
						
						//do this for every particle in aGroup
						for(auto& aP:aGroup)
						{
							double deltaHdeltaX=1.1;
							//double patchSize=initPatchSize;
							//while(deltaHdeltaX>1.0)
							
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
								//closest edge
								if(d.x*d.x+d.y*d.y+d.z*d.z<8.0 && dotProduct(normalGroup[i],*myNormal)>0)
									normalSum+=normalGroup[i];
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
								//closest edge
								if(d.x*d.x+d.y*d.y+d.z*d.z<8.0 && dotProduct(nearbyNormals[i],*myNormal)>0)
									normalSum+=nearbyNormals[i];
							}
							
							//normalSum/=static_cast<double>(nNormals);
							
							//get the projection matrix, ROW_MAJOR
							threeVector<threeVector<double> > M;
							//M=I-(aP*aP^T)/aP^2
							//position<double> aPu=unitVector(aP);
							threeVector<double> aPu=unitVector(normalSum);
							
							if(dotProduct(aPu,oneZ)<0)
								oneZ.z=-oneZ.z;
							
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
							
							//our z oriented patch (heights)
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
								
							//transform the neighbors by the projection matrix and place on grid
							auto nearbyNormal=nearbyNormals.begin();
							
							for(auto &nG:nearbyGroups)
							{
								threeVector<double> sG,qG;
								
								qG.x=nG.x-aP.x;
								qG.y=nG.y-aP.y;
								qG.z=nG.z-aP.z;
								
								qG.x-=(qG.x>size.x/2.0)?size.x:0;
								qG.y-=(qG.y>size.y/2.0)?size.y:0;
								qG.z-=(qG.z>size.z/2.0)?size.z:0;
								qG.x+=(qG.x<-size.x/2.0)?size.x:0;
								qG.y+=(qG.y<-size.y/2.0)?size.y:0;
								qG.z+=(qG.z<-size.z/2.0)?size.z:0;
								
								sG.x=qG.x*M.x.x+qG.y*M.x.y+qG.z*M.x.z;
								sG.y=qG.x*M.y.x+qG.y*M.y.y+qG.z*M.y.z;
								sG.z=qG.x*M.z.x+qG.y*M.z.y+qG.z*M.z.z;
								
								
								
								sG.z*=oneZ.z;
								//we are here
								
								threeVector<int> patchCoor=-1.0;
								if(sG.x>-1.5*patchSize && sG.x<1.5*patchSize)
									patchCoor.x=static_cast<int>((sG.x+patchSize*1.5)/(patchSize));
								if(sG.y>-1.5*patchSize && sG.y<1.5*patchSize)
									patchCoor.y=static_cast<int>((sG.y+patchSize*1.5)/(patchSize));
								
								//farthest corner
								//if(sG.x*sG.x+sG.y*sG.y+sG.z*sG.z<4.5*patchSize*patchSize)
								if(sG.z<patchSize && sG.z>-patchSize)
								{
									
									if(dotProduct(*nearbyNormal,aPu)>0)
									{
										
										if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
										{
											zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
											nZPatch.s[patchCoor.x].s[patchCoor.y]++;
										}
									}
								}
								nearbyNormal++;
							}
							
							nearbyNormal=normalGroup.begin();
							for(auto &nG:aGroup)
							{
								threeVector<double> sG,qG;
								
								qG.x=nG.x-aP.x;
								qG.y=nG.y-aP.y;
								qG.z=nG.z-aP.z;
								
								qG.x-=(qG.x>size.x/2.0)?size.x:0;
								qG.y-=(qG.y>size.y/2.0)?size.y:0;
								qG.z-=(qG.z>size.z/2.0)?size.z:0;
								qG.x+=(qG.x<-size.x/2.0)?size.x:0;
								qG.y+=(qG.y<-size.y/2.0)?size.y:0;
								qG.z+=(qG.z<-size.z/2.0)?size.z:0;
								
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
								//farthest corner
								//if(sG.x*sG.x+sG.y*sG.y+sG.z*sG.z<4.5*patchSize*patchSize)// && dotProduct(*nearbyNormal,aPu)>0)
								if(sG.z<patchSize && sG.z>-patchSize)
								{
									
									if(dotProduct(*nearbyNormal,aPu)>0)
									{
										
										if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
										{
											zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
											nZPatch.s[patchCoor.x].s[patchCoor.y]++;
										}
									}
								}
								nearbyNormal++;
							}
							//throw 0;
							bool badPatch=false;
							int nL=0;
							for(int a=0;a<3;a++)
							{
								for(int b=0;b<3;b++)
								{
									nL+=nZPatch.s[a].s[b];
									//if(nZPatch.s[a].s[b]==0 && a!=b && !(a==0 && b==2) && !(a==2 && b==0))
									//{
									//	nBadPatches++;
										//badPatch=true;
									//}
								}
							}
							nL+=1;
							
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
							
				/*for(int a=0;a<3;a++)
				{
				for(int b=0;b<3;b++)
				{
					std::cout << zPatch.s[a].s[b] << std::endl;
				//if(zPatch.s[a].s[b]>patchSize || zPatch.s[a].s[b]<-patchSize)
				{
					//std::cout << a+b*3 << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << '\n';
				}
				}
				}
				*/
							
							zPatch.s[1].s[1]=0.0;
							
							
							
							if(!badPatch)
							{
								//double myGradSquared=0.66689*pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
								//		zPatch.s[2].s[1]+zPatch.s[0].s[1])/2.0+
								//		(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
								//		zPatch.s[2].s[2]+zPatch.s[0].s[0])/4.0-
								//		3.0*zPatch.s[1].s[1])/(patchSize*patchSize),2.0);//average 0.66689 nm^2 per lipid
								
								double patchArea=patchSize*patchSize;
								
								double area=sqrt(1.0+(zPatch.s[0].s[1]*zPatch.s[0].s[1]+zPatch.s[1].s[0]*zPatch.s[1].s[0])/patchArea);
								
								double myLaplacian=((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
								zPatch.s[2].s[1]+zPatch.s[0].s[1])*2.0/3.0+
								(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
								zPatch.s[2].s[2]+zPatch.s[0].s[0])/6.0-
								(10.0/3.0)*zPatch.s[1].s[1])/(patchArea);
								
								
								
								double myGradSquared=pow(myLaplacian,2.0)*area;
								
								//double myGradSquared=pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
								//zPatch.s[2].s[1]+zPatch.s[0].s[1])*2.0/3.0+
								//(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
								//zPatch.s[2].s[2]+zPatch.s[0].s[0])/6.0-
								//(10.0/3.0)*zPatch.s[1].s[1])/(patchArea),2.0);//average 0.66689 nm^2 per lipid
								//double myGradSquared=pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
								//zPatch.s[2].s[1]+zPatch.s[0].s[1])-
								//4.0*zPatch.s[1].s[1])/(patchArea),2.0);//average 0.66689 nm^2 per lipid
								
								gradSquared+=myGradSquared;
								laplacianInt+=myLaplacian;
								gradSquareValues.push_back(myGradSquared);
								if(nVertSpacing<0)
								{
									fourVector<double> aV;
									aV.x=aP.x;
									aV.y=aP.y;
									aV.z=aP.z;
									aV.t=myGradSquared;
									gradSquarePos.push_back(aV);
									
								}
							}
							myNormal++;
						}
						//int gradSquareValuesEnd=gradSquareValues.size();
						
						double avgGradSquare=(gradSquared-gradSquaredLow)/static_cast<double>(gradSquareValues.size()-gradSquareValuesStart);
						double rms=0;
						for(int i=gradSquareValuesStart;i<gradSquareValues.size();i++)
							rms+=pow(gradSquareValues[i]-avgGradSquare,2.0);
						if(gradSquareValues.size()-gradSquareValuesStart>0)
							localRms+=sqrt(rms/static_cast<double>(gradSquareValues.size()-gradSquareValuesStart));
						
					}
				}
				//throw 0;
				std::cerr << percentCompletion << "%\n" << std::endl;
				if(filt.size()>10)
				{
					if(nVertSpacing<0)
					{
						std::string buf;
						buf="heatMap_SN";
						std::string whatever;
						std::stringstream parseNumber;
						parseNumber << s << "_" << patchSize;
						parseNumber >> whatever;
						
						buf+=whatever;
						buf+="_";
						//buf="patchMin_";
						buf+=argv[1];
						buf+=".xyz";
						std::fstream dataFile;
						
						dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
						//auto gradSquareValue=gradSquareValues.begin();
						
						double avgGradSquare=gradSquared/static_cast<double>(gradSquareValues.size());
						double minGrad=avgGradSquare;
						double maxGrad=avgGradSquare;
						for(auto &value:gradSquareValues)
						{
							if(value<minGrad)
								minGrad=value;
							if(value>maxGrad)
								maxGrad=value;
						}
						double dGrad=(maxGrad-minGrad)/20.0;
						
						
						
						dataFile << gradSquareValues.size() << "\nalksdfj\n";
						
						auto cPos=gradSquarePos.begin();
						for(auto &value:gradSquareValues)
						{
							dataFile << floor((value-minGrad)/dGrad) << ' ' << cPos->x << ' ' << cPos->y << ' ' << cPos->z << '\n';
							cPos++;
						}
						dataFile << std::endl;
						dataFile.close();
						
						
						buf.clear();
						buf="gradSorted_SN";
						
						buf+=whatever;
						buf+="_";
						//buf="patchMin_";
						buf+=argv[1];
						buf+=".dat";
						
						dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
						
						std::sort(gradSquareValues.begin(), gradSquareValues.end());
						for(auto &value:gradSquareValues)
						{
							dataFile << value << '\n';
						}
						dataFile << std::endl;
						dataFile.close();
						
					}
					if(nVertSpacing>1)
					{
						std::string buf;
						buf="patchProfile_SN";
						std::string whatever;
						std::stringstream parseNumber;
						parseNumber << s;
						parseNumber >> whatever;
						
						buf+=whatever;
						buf+="_";
						//buf="patchMin_";
						buf+=argv[1];
						buf+=".dat";
						std::fstream dataFile;
						
						dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
						//auto gradSquareValue=gradSquareValues.begin();
						
						
						double avgGradSquare=gradSquared/static_cast<double>(gradSquareValues.size());
						double rms=0;
						for(auto &value:gradSquareValues)
							rms+=pow(value-avgGradSquare,2.0);
						rms=sqrt(rms/static_cast<double>(gradSquareValues.size()));
						localRms/=static_cast<double>(gradSquareValues.size());
						if(gradSquareValues.size()!=0)
						{
							dataFile << patchSize << ' ' << gradSquared << ' ' << rms << ' '; 
							dataFile << localRms << ' ' << gradSquareValues.size() << ' ' << laplacianInt << std::endl;
						}
						
						dataFile.close();
					}
					else if (nVertSpacing==1 || nVertSpacing==0)
					{
						std::string buf;
						buf="gradSquare_SN";
						std::string whatever;
						std::stringstream parseNumber;
						parseNumber << s << "_" << patchSize;
						parseNumber >> whatever;
						
						buf+=whatever;
						buf+="_";
						//buf="patchMin_";
						buf+=argv[1];
						buf+=".dat";
						std::fstream dataFile;
						
						dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
						//auto gradSquareValue=gradSquareValues.begin();
						
						double avgGradSquare=gradSquared/static_cast<double>(gradSquareValues.size());
						double rms=0;
						for(auto &value:gradSquareValues)
							rms+=pow(value-avgGradSquare,2.0);
						rms=sqrt(rms/static_cast<double>(gradSquareValues.size()));
						localRms/=static_cast<double>(gradSquareValues.size());
						if(gradSquareValues.size()!=0)
							dataFile << time << ' ' << gradSquared << ' ' << rms << ' ' << ' ';
							dataFile << localRms << ' ' << gradSquareValues.size() << ' ' << laplacianInt << std::endl;
						
						dataFile.close();
					}
				}
			}
		}
		currentFrame++;
		if(currentFrame>frame && frame>0)
			break;
	}
	return 0;
}
