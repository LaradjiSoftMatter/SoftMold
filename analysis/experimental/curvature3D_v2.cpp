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
		std::cerr << argv[0] << " name surfSpacing minAngle vertShellRatio vertSpacing dry type1 type2 ...\n";
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
	//double maxSphere=sqrt(4.0*2.0*surfSpacingSqr);//sqrt((2.0*surfSpacing)^2*dimension)
	//double maxSphereSqr=maxRad*maxRad;
	//double minSphereSqr=minRad*minRad;
	
	//surfSkip=1;//floor(2.0*maxRad/surfSpacing);
	
	//std::cout << surfSkip << std::endl;
	
	std::vector<int> types;
	
	for(int i=7;i<argc;i++)
	{
		int type;
		cmdArg >> type;
		types.push_back(type);
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
	std::vector<int> pIndex;
	
	
	
	for(int i=0;i<System.readNParticles();i++)
	{
		bool match=false;
		for(int tIndex=0;!match && tIndex<types.size();tIndex++)
			match=(p[i].type==types[tIndex]);
		
		if(match)
			pIndex.push_back(i);
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
		
		
		std::vector< std::vector< position<double> > > cellList;
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
		
		std::map<int,cellType > cells;
		
		for(int j=0;j<pIndex.size();j++)
		{
			int i=pIndex[j];
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
				it->second.push_back(j);
			}
			else
			{
				//initialize operation
				cellType buf;
				buf.push_back(p[i]);
				//cells[hash]=buf;
				cells.insert(it, std::map<int,cellType >::value_type(hash,buf));
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
		std::vector< cellType > surfaces;
		
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
				cellType surface;
				
				while(stack.size()>0)
				{
					auto currentPoint=stack.back();
					
					int currentHash=currentPoint->first;
					
					position<double> aP=currentPoint->second;
					
					threeVector<int> iCoor=unhash(currentHash,surfN);
					
					for(auto& cellPoint:cells[currentHash])
						surface.push_back(cellPoint);
						
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
						if(it!=pFilt.end() && it->first==neighborHash && neighborHash!=currentHash)
						{
							if(it->second.type==-1)
							{
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
								if(d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr*2.75)
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
		
		std::vector<std::map<int,cellType > > surfCells;
		for(auto& surface:surfaces)
		{
			std::map<int,cellType > sCells;
			//std::sort(surface.begin(),surface.end());
			for(auto& j:surface)
			{
				int i=pIndex[j];
				threeVector<int> c;
				c.x=static_cast<int>(p[i].x/vertR.x)%vertN.x;
				c.y=static_cast<int>(p[i].y/vertR.y)%vertN.y;
				c.z=static_cast<int>(p[i].z/vertR.z)%vertN.z;
				int hash=c.x+c.y*vertN.x+c.z*vertN.x*vertN.y;
				
				std::map<int,cellType >::iterator it=sCells.lower_bound(hash);
				if(it!=sCells.end() && it->first==hash)//!(cells.key_comp()(it->first, hash)))
				{
					//push operation
					it->second.push_back(j);
				}
				else
				{
					//initialize operation
					cellType buf;
					buf.push_back(j);
					//cells[hash]=buf;
					sCells.insert(it, std::map<int,cellType >::value_type(hash,buf));
				}
			}
			//nSurface++;
			surfCells.push_back(sCells);
		}
		
		if(surfaces.size()<4)
		{
			for(auto &sCells:surfCells)
				std::cerr << "using " << sCells.size() << " points for surface triangulation" << std::endl;
		}
		
		//examined surfaces for triangulation
		typedef std::map<int, position<double> >  filtType;
		std::vector<filtType> vFilt;
		
		for(auto& sCells:surfCells)
		{
			filtType filt;
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
				}
				pAvg.x/=static_cast<double>(it.second.size());
				pAvg.y/=static_cast<double>(it.second.size());
				pAvg.z/=static_cast<double>(it.second.size());
				pAvg.type=-1;
				
				filt[gridPos]=pAvg;
			}
			vFilt.push_back(filt);
		}
		throw 0;
		
		if(surfaces.size()<4)
		{
			for(auto &filt:vFilt)
				std::cerr << "using " << filt.size() << " points for surface triangulation" << std::endl;
			
		}
		
		//triangles, each is composed of integral s[3], with each element a hashed value from sCell
		typedef threeVector<cellType::iterator >  triangleType;
		typedef std::vector<triangleType> triangleGroupType;
		std::vector< std::map<int,triangleGroupType> > sTriangles;
		
		//an empty iterator to check if a triangle vertex is blank
		//apparantly, there is no null or default value I can use, so constructed this
		cellType emptyCellContainer;
		cellType::iterator emptyCell=emptyCellConter.end();
		
		//try every surface
		for(int s=0;s<vFilt.size();s++)
		{
			#define unusedPoint -1
			#define vertexPoint 0
			#define edgePoint 1
			#define ignorePoint 2
			const& auto filt=vFilt[s];//for nearest neighbor traversal
			const& auto sCells=surfCells[s];//for vertex candidates
			
			//this will have duplicates in it if vertices cross our hash
			std::map<int,triangleGroupType> triangles;
			
			//unaccounted surface
			for(auto element=filt.begin();element!=filt.end();element++)
			{
				//unaccounted surface element?
				if(element->second.type==unusedPoint)
				{
					//grab a seed triangle from the first element we found for surface passing through cell
					//the other vertices are incomplete
					triangleType seedTriangle;
					seedTriangle.s[0]=sCells[element->first].begin();
					seedTriangle.s[1]=emptyCell;
					seedTriangle.s[2]=emptyCell;
					
					triangleGroup seedGroup;
					seedGroup.push_back(seedTriangle);
					
					//triangles[element->first]=seedGroup;
					std::map<int,triangleGroupType>::iterator triangleSeedIt=triangles.lower_bound(element->first);
					triangleSeedIt=triangles.insert(triangleSeedIt,element->first);
					
					element->second.type=vertexPoint;//traversed
					std::deque<std::map<int,triangleGroupType>::iterator> stack;//our boundary
					stack.push_back(triangleSeedIt);
					
					while(stack.size()>0)
					{
						//grab the current group from stack, we will eliminate it if it can't be completed
						std::map<int,triangleGroupType>::iterator currentGroup=stack.back();
						stack.pop_back();
						
						int currentHash=currentGroup->first;
						
						threeVector<int> iCoor=unhash(currentHash,vertN);
						
						//mark it as a vertex
						filt[neighborHash].second.type=vertexPoint;
						
						//candidate triangles
						//std::vector<triangleType> tCandidates;
						
						//candidate points for triangles, these are our edges with the current triangles vertices
						//std::vector<cellType::iterator> candidates;
						
						//we will start making candidate triangles and edges with this 
						triangleGroupType candidateEdges;
						
						//start with the first one
						//candidates.push_back(cell[element->first].second[0]);
						
						bool edgeMatch=false;
						//find candidates around current group
						for(int k=0;k<27;k++)
						{
							threeVector<int> nCoor=iCoor;
							nCoor.x+=(k%3-1);
							nCoor.y+=(int(k/3)%3-1);
							nCoor.z+=(int(k/9)-1);
							
							int neighborHash=hash(nCoor,vertN);
							
							auto it=sCells.lower_bound(neighborHash);
							if(it!=sCells.end() && it->first==neighborHash)// && neighborHash!=currentHash)
							{
								//test a point in the neighbor group
								for(auto pC=it->second.begin();pC!=it->second.end();pC++)
								{
									//against every vertex in current group
									for(triangleGroupType::iterator triangle=currentGroup->second.begin();
										triangle!=currentGroup->second.end();triangle++)
									for(int vertex=0;vertex<3;vertex++)
									{
										if(triangle.s[vertex]!=emptyCell)
										{
											position<double> neP=(*triangle->s[vertex]);
											threeVector<double> d;
											d.x=neP.x-pC->second.x;
											d.y=neP.y-pC->second.y;
											d.z=neP.z-pC->second.z;
											
											while(d.x>size.x/2.0)d.x-=size.x;
											while(d.y>size.y/2.0)d.y-=size.y;
											while(d.z>size.z/2.0)d.z-=size.z;
											while(d.x<-size.x/2.0)d.x+=size.x;
											while(d.y<-size.y/2.0)d.y+=size.y;
											while(d.z<-size.z/2.0)d.z+=size.z;
										
											double dr=d.x*d.x+d.y*d.y+d.z*d.z;
											
											//this is a thickish shell
											if(dr>vertSpacingSqrMin && dr<vertSpacingSqrMax)
											{
												triangleType edge;
												edge.x=pC;//from candidate group
												edge.y=triangle->s[vertex];//from this triangle
												edge.z=emptyCell;//We don't know yet
												edges.push_back(edge);
												candidateEdges.push_back(triangle);
												edgeMatch=true;
											}
										}
									}
								}
							}
						}
						
						bool vertexMatch=false;
						//check candidates against other vertices
						//find candidates around current group
						for(int k=0;k<27 && candidates.size()>0;k++)
						{
							threeVector<int> nCoor=iCoor;
							nCoor.x+=(k%3-1);
							nCoor.y+=(int(k/3)%3-1);
							nCoor.z+=(int(k/9)-1);
							
							int neighborHash=hash(nCoor,vertN);
							
							auto neighGroup=triangles.lower_bound(neighborHash);
							if(neighGroup!=triangles.end() && neighGroup->first==neighborHash)// && neighborHash!=currentHash)
							{
								//test a point in the neighbor group
								for(auto aEdge=candidateEdges.begin();aEdge!=candidateEdges.end();aEdge++)
								{
									position<double> aP=(*aEdge.x)->second;
									position<double> bP=(*aEdge.y)->second;
									threeVector<double> da;
									da.x=aP.x-bP.x;
									da.y=aP.y-bP.y;
									da.z=aP.z-bP.z;
												
									while(da.x>size.x/2.0)da.x-=size.x;
									while(da.y>size.y/2.0)da.y-=size.y;
									while(da.z>size.z/2.0)da.z-=size.z;
									while(da.x<-size.x/2.0)da.x+=size.x;
									while(da.y<-size.y/2.0)da.y+=size.y;
									while(da.z<-size.z/2.0)da.z+=size.z;
									
									//against every vertex in neighbor group
									for(triangleGroupType::iterator triangle=neighGroup->second.begin();
										triangle!=neighGroup->second.end();triangle++)
									{
										for(int vertex=0;vertex<3;vertex++)
										{
											if(triangle->s[vertex]!=emptyCell)
											{
												position<double> cP=(*triangle->s[vertex]);
												threeVector<double> db;
												db.x=cP.x-aP.x;
												db.y=cP.y-aP.y;
												db.z=cP.z-aP.z;
															
												while(db.x>size.x/2.0)db.x-=size.x;
												while(db.y>size.y/2.0)db.y-=size.y;
												while(db.z>size.z/2.0)db.z-=size.z;
												while(db.x<-size.x/2.0)db.x+=size.x;
												while(db.y<-size.y/2.0)db.y+=size.y;
												while(db.z<-size.z/2.0)db.z+=size.z;
												
												threeVector<double> dc;
												dc.x=bP.x-cP.x;
												dc.y=bP.y-cP.y;
												dc.z=bP.z-cP.z;
														
												while(dc.x>size.x/2.0)dc.x-=size.x;
												while(dc.y>size.y/2.0)dc.y-=size.y;
												while(dc.z>size.z/2.0)dc.z-=size.z;
												while(dc.x<-size.x/2.0)dc.x+=size.x;
												while(dc.y<-size.y/2.0)dc.y+=size.y;
												while(dc.z<-size.z/2.0)dc.z+=size.z;
												
												threeVector<double> scal=circumscribe(da,db,dc);
											
												position<double> sphere;
												
												sphere.x=scal.x*aP.x+scal.y*bP.x+scal.z*cP.x;
												sphere.y=scal.x*aP.y+scal.y*bP.y+scal.z*cP.y;
												sphere.z=scal.x*aP.z+scal.y*bP.z+scal.z*cP.z;
												
												threeVector<double> d;
												
												d.x=sphere.x-aP.x;
												d.y=sphere.y-aP.y;
												d.z=sphere.z-aP.z;
												
												double sRadiusSqr=d.x*d.x+d.y*d.y+d.z*d.z-0.0001;
												
												for(triangleGroupType::iterator dEdge=neighGroup->second.begin();
													dEdge!=neighGroup->second.end();dEdge++)
												{
													for(int vertex=0;vertex<3;vertex++)
													{
														if(triangle->s[vertex]!=emptyCell)
														{
															
														}
													}
												}
												threeVector<double> d;
												d.x=neP.x-currentEdge.y->second.x;
												d.y=neP.y-currentEdge.y->second.y;
												d.z=neP.z-currentEdge.y->second.z;
												
												while(d.x>size.x/2.0)d.x-=size.x;
												while(d.y>size.y/2.0)d.y-=size.y;
												while(d.z>size.z/2.0)d.z-=size.z;
												while(d.x<-size.x/2.0)d.x+=size.x;
												while(d.y<-size.y/2.0)d.y+=size.y;
												while(d.z<-size.z/2.0)d.z+=size.z;
												
												double dr=d.x*d.x+d.y*d.y+d.z*d.z;
												if(dr>vertSpacingSqrMin && dr<vertSpacingSqrMax)
												{
													aEdge.z=triangle->s[vertex];
													vertexMatch=true;
												}
											}
										}
									}
								}
							}
						}
						//Conditions:
						//0. No match, overlap eliminates candidates from triangulation, (edgeMatch || vertexMatch)==false
						//	a. Go to next on stack
						//1. Matches one current triangle groups (edgeMatch==true && vertexMatch=false)
						//	a. 2nd Vertex only {Failed, delete triangle}, candidatesEdges.size()==1
						//	b. 2nd and 3rd Vertex {completes all edges, optimal spacing ~equilateral}, candidateEdges.size()>1
						//2. Matches one or more neighbor triangle groups, but no current group (never going to happen)
						//	a. no additions left, try next on stack
						//3. Matches one or more groups including at least one current group
						//	a. Forms its own group and completes a current group, check delauny
						//	b. Doesn't complete a current group, delete 
						//4. Check for circularity in current group, if true, don't add to stack
						//	else add to stack
						//5. Go to next on stack
						
						bool newGroup=false;
						bool circular=false;
						//get equilateral spacing
						if(edgeMatch && !vertexMatch && candidateEdges.size()>1)
						{
							std::vector<triangleGroupType> nearestX,nearestY,nearestZ;
							double nearestAngle;//we want this to be ~acos(M_PI/3)
							//remember,this needs to be true: (*aEdge.x)==(*bEdge.x) !
							for(auto aEdge=candidateEdges.begin();aEdge!=candidateEdges.end();aEdge++)
							{
								threeVector<double> da;
								da.x=(*aEdge.x)->second.x-(*aEdge.y)->second.x;
								da.y=(*aEdge.x)->second.y-(*aEdge.y)->second.y;
								da.z=(*aEdge.x)->second.z-(*aEdge.y)->second.z;
											
								while(da.x>size.x/2.0)da.x-=size.x;
								while(da.y>size.y/2.0)da.y-=size.y;
								while(da.z>size.z/2.0)da.z-=size.z;
								while(da.x<-size.x/2.0)da.x+=size.x;
								while(da.y<-size.y/2.0)da.y+=size.y;
								while(da.z<-size.z/2.0)da.z+=size.z;
								
								for(auto bEdge=aEdge+1;bEdge!=candidateEdges.end();bEdge++)
								{
									if((*aEdge.x)==(*bEdge.x))
									{
										threeVector<double> db;
										db.x=(*bEdge.x)->second.x-(*bEdge.y)->second.x;
										db.y=(*bEdge.x)->second.y-(*bEdge.y)->second.y;
										db.z=(*bEdge.x)->second.z-(*bEdge.y)->second.z;
													
										while(db.x>size.x/2.0)db.x-=size.x;
										while(db.y>size.y/2.0)db.y-=size.y;
										while(db.z>size.z/2.0)db.z-=size.z;
										while(db.x<-size.x/2.0)db.x+=size.x;
										while(db.y<-size.y/2.0)db.y+=size.y;
										while(db.z<-size.z/2.0)db.z+=size.z;
										
										threeVector<double> dc;
										dc.x=(*aEdge.y)->second.x-(*bEdge.y)->second.x;
										dc.y=(*aEdge.y)->second.y-(*bEdge.y)->second.y;
										dc.z=(*aEdge.y)->second.z-(*bEdge.y)->second.z;
													
										while(dc.x>size.x/2.0)dc.x-=size.x;
										while(dc.y>size.y/2.0)dc.y-=size.y;
										while(dc.z>size.z/2.0)dc.z-=size.z;
										while(dc.x<-size.x/2.0)dc.x+=size.x;
										while(dc.y<-size.y/2.0)dc.y+=size.y;
										while(dc.z<-size.z/2.0)dc.z+=size.z;
										
										threeVector<double> scal=circumscribe(da,db,dc);
										
										position<double> sphere;
										
										sphere.x=scal.x*aP.x+scal.y*bP.x+scal.z*cP.x;
										sphere.y=scal.x*aP.y+scal.y*bP.y+scal.z*cP.y;
										sphere.z=scal.x*aP.z+scal.y*bP.z+scal.z*cP.z;
										
										threeVector<double> d;
										
										d.x=sphere.x-aP.x;
										d.y=sphere.y-aP.y;
										d.z=sphere.z-aP.z;
										
										double sRadiusSqr=d.x*d.x+d.y*d.y+d.z*d.z-0.0001;
										
										for(
									}
								}
							}
							
						}
						if(edgeMatch && vertexMatch)
						{
							
						}
						
						//check if a point already has neighboring triangles
						auto triNeighIt=triangles.lower_bound(neighborHash);
								
						//if the point hasn't been tested, insert new neighbor list
						if(triNeighIt==triangles.end() || triNeighIt->first!=neighborHash)
							triNeighIt=triangles.insert(triNeighIt, std::map<int,triangleGroupType >::value_type(neighborHash,triangleGroupType()));
						
						//try deluanay triangulation with candidate edges
						//This just finds all potential surface types of the form:
						//   a
						//     |> c
						//   b 
						//This is only the first triangulation pass assuming a is a corePoint
						//for(auto& ab:boundary)
						for(auto b=candidates.begin();b!=candidates.end();b++)
							{
								//try ab from boundary and c from stack
								//position<double> aP=ab.x;
								//position<double> bP=ab.y;
								position<double> bP=(*b)->second;
								
								for(auto c=candidates.begin();c!=candidates.end();c++)
								{
									position<double> cP=(*c)->second;
									
									if((*b)->first != (*c)->first)
									{
										bool accept=true;
										//int accept=0;
										
										threeVector<double> dab,dbc,dca;
										//ab
										dab.x=aP.x-bP.x;
										dab.y=aP.y-bP.y;
										dab.z=aP.z-bP.z;
										while(dab.x>size.x/2.0)dab.x-=size.x;
										while(dab.y>size.y/2.0)dab.y-=size.y;
										while(dab.z>size.z/2.0)dab.z-=size.z;
										while(dab.x<-size.x/2.0)dab.x+=size.x;
										while(dab.y<-size.y/2.0)dab.y+=size.y;
										while(dab.z<-size.z/2.0)dab.z+=size.z;
										//bc
										dbc.x=bP.x-cP.x;
										dbc.y=bP.y-cP.y;
										dbc.z=bP.z-cP.z;
										while(dbc.x>size.x/2.0)dbc.x-=size.x;
										while(dbc.y>size.y/2.0)dbc.y-=size.y;
										while(dbc.z>size.z/2.0)dbc.z-=size.z;
										while(dbc.x<-size.x/2.0)dbc.x+=size.x;
										while(dbc.y<-size.y/2.0)dbc.y+=size.y;
										while(dbc.z<-size.z/2.0)dbc.z+=size.z;
										//ca
										dca.x=cP.x-aP.x;
										dca.y=cP.y-aP.y;
										dca.z=cP.z-aP.z;
										while(dca.x>size.x/2.0)dca.x-=size.x;
										while(dca.y>size.y/2.0)dca.y-=size.y;
										while(dca.z>size.z/2.0)dca.z-=size.z;
										while(dca.x<-size.x/2.0)dca.x+=size.x;
										while(dca.y<-size.y/2.0)dca.y+=size.y;
										while(dca.z<-size.z/2.0)dca.z+=size.z;
										
										threeVector<double> scal=circumscribe(dab,dbc,dca);
										
										position<double> sphere;
										
										sphere.x=scal.x*aP.x+scal.y*bP.x+scal.z*cP.x;
										sphere.y=scal.x*aP.y+scal.y*bP.y+scal.z*cP.y;
										sphere.z=scal.x*aP.z+scal.y*bP.z+scal.z*cP.z;
										
										threeVector<double> d;
										
										d.x=sphere.x-aP.x;
										d.y=sphere.y-aP.y;
										d.z=sphere.z-aP.z;
										
										double sRadiusSqr=d.x*d.x+d.y*d.y+d.z*d.z-0.0001;
										
										//make sure triangles have angles of at least minAngle
										double abcAngle=acos(dotProduct(dab,dbc)/(magnitude(dab)*magnitude(dbc)));
										double bcaAngle=acos(dotProduct(dbc,dca)/(magnitude(dbc)*magnitude(dca)));
										double cabAngle=acos(dotProduct(dca,dab)/(magnitude(dca)*magnitude(dab)));
										
										//std::cout << abcAngle << ' ' << bcaAngle << ' ' << cabAngle << std::endl;
										
										if((abcAngle<minAngle || bcaAngle<minAngle || cabAngle<minAngle) ||
											dab.x*dab.x+dab.y*dab.y+dab.z*dab.z>vertSpacingSqr*3.0 ||
											dbc.x*dbc.x+dbc.y*dbc.y+dbc.z*dbc.z>vertSpacingSqr*3.0 ||
											dca.x*dca.x+dca.y*dca.y+dca.z*dca.z>vertSpacingSqr*3.0 ||
											sphere.x!=sphere.x || sphere.y!=sphere.y || sphere.z!=sphere.z)// ||
											//2.25*vertSpacingSqr*static_cast<double>(1.0)<sRadiusSqr)
										{
											accept=false;
										}
										else for(auto& test:excludeRange)
										{
std::cout << 2 << ' ' << (*test).second.x << ' ' << (*test).second.y << ' ' << (*test).second.z << std::endl;
											d.x=sphere.x-(*test).second.x;
											d.y=sphere.y-(*test).second.y;
											d.z=sphere.z-(*test).second.z;
											
											while(d.x>size.x/2.0)d.x-=size.x;
											while(d.y>size.y/2.0)d.y-=size.y;
											while(d.z>size.z/2.0)d.z-=size.z;
											while(d.x<-size.x/2.0)d.x+=size.x;
											while(d.y<-size.y/2.0)d.y+=size.y;
											while(d.z<-size.z/2.0)d.z+=size.z;
											
											if(d.x*d.x+d.y*d.y+d.z*d.z<sRadiusSqr)
											{
												accept=false;
												//accept++;
												break;
											}
										}
										if(accept)//surfSkip)
										{
											threeVector<std::map<int, position<double> >::iterator> triangle;
											triangle.x=currentPoint;
											triangle.y=*b;
											triangle.z=*c;
											
//std::cout << 1 << ' ' << sphere.x << ' ' << sphere.y << ' ' << sphere.z << std::endl;
											
											//threeVector<int> triangle;
											//triangle.x=currentHash;
											//triangle.y=(*b)->first;
											//triangle.z=(*c)->first;
											
											//check if we have duplicates
											bool duplicate=false;
											for(auto& aTriangle:tCandidates)
											{
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
												
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
											}
											//duplicate=false;
											/*
											for(auto& aTriangles:neighbors[currentHash]) {
												const auto& aTriangle=triangles[aTriangles];
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
												
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
											}
											
											for(auto& bTriangles:neighbors[(*b)->first]) {
												const auto& aTriangle=triangles[bTriangles];
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
												
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
												
											}
											
											for(auto& cTriangles:neighbors[(*c)->first]) {
												const auto& aTriangle=triangles[cTriangles];
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
												
												if(aTriangle.x->first==triangle.z->first && aTriangle.y->first==triangle.y->first && aTriangle.z->first==triangle.x->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.y->first && aTriangle.y->first==triangle.x->first && aTriangle.z->first==triangle.z->first)
													duplicate=true;
												if(aTriangle.x->first==triangle.x->first && aTriangle.y->first==triangle.z->first && aTriangle.z->first==triangle.y->first)
													duplicate=true;
												
											}
											*/
											
											if(!duplicate)
											{
												tCandidates.push_back(triangle);
												/*
												//check normal vector orientation
												threeVector<double> myNormal=unitVector(crossProduct(dab,dbc));
												int flip=0;
												//int nNormals=0;
												bool isGreaterThanMinAngle=true;//false;
												for(auto& nTriangles:neighbors[currentHash]) 
												{
													const auto& aTriangle=triangles[nTriangles];
													position<double> aPb=aTriangle.x->second;
													bP=aTriangle.y->second;
													cP=aTriangle.z->second;
													
													dab.x=aPb.x-bP.x;
													dab.y=aPb.y-bP.y;
													dab.z=aPb.z-bP.z;
													while(dab.x>size.x/2.0)dab.x-=size.x;
													while(dab.y>size.y/2.0)dab.y-=size.y;
													while(dab.z>size.z/2.0)dab.z-=size.z;
													while(dab.x<-size.x/2.0)dab.x+=size.x;
													while(dab.y<-size.y/2.0)dab.y+=size.y;
													while(dab.z<-size.z/2.0)dab.z+=size.z;
													//bc
													dbc.x=bP.x-cP.x;
													dbc.y=bP.y-cP.y;
													dbc.z=bP.z-cP.z;
													while(dbc.x>size.x/2.0)dbc.x-=size.x;
													while(dbc.y>size.y/2.0)dbc.y-=size.y;
													while(dbc.z>size.z/2.0)dbc.z-=size.z;
													while(dbc.x<-size.x/2.0)dbc.x+=size.x;
													while(dbc.y<-size.y/2.0)dbc.y+=size.y;
													while(dbc.z<-size.z/2.0)dbc.z+=size.z;
													
													threeVector<double> normal=unitVector(crossProduct(dab,dbc));
													double orient=dotProduct(myNormal,normal);
													if(orient<0)
													{
														normal.x=-normal.x;
														normal.y=-normal.y;
														normal.z=-normal.z;
													}
													//if(acos(dotProduct(normal,myNormal))>vertShellRatio)
													//	isGreaterThanMinAngle=true;
													if(orient<0)
													{
														flip++;
														break;
													}
													//nNormals++;
												}
												if(flip>0)
												{
													triangle.z=currentPoint;
													//triangle.y=*b;
													triangle.x=*c;
													//myNormal.x=-myNormal.x;
													//myNormal.y=-myNormal.y;
													//myNormal.z=-myNormal.z;
												}
												*/
												
												//if(isGreaterThanMinAngle)
												//{
												//	neighbors[currentHash].push_back(triangles.size());
												//	neighbors[(*b)->first].push_back(triangles.size());
												//	neighbors[(*c)->first].push_back(triangles.size());
													//triangles.push_back(triangle);
std::cout << 1 << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << std::endl;
std::cout << 1 << ' ' << bP.x << ' ' << bP.y << ' ' << bP.z << std::endl;
std::cout << 1 << ' ' << cP.x << ' ' << cP.y << ' ' << cP.z << std::endl;
												//}
												//normals.push_back(myNormal);
												//boundary.push_back(currentPoint);
												//boundary.push_back(*b);
												//boundary.push_back(*c);
											}
										}
										else
										{
std::cout << 0 << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << std::endl;
std::cout << 0 << ' ' << bP.x << ' ' << bP.y << ' ' << bP.z << std::endl;
std::cout << 0 << ' ' << cP.x << ' ' << cP.y << ' ' << cP.z << std::endl;
										}
										
									}
								}
							}
							//throw 0;
						/*
						int checkLayer=1;
						bool nonAdded;
						do
						{
							nonAdded=true;
							int kLength=checkLayer*2+1;
							//int kBlock=kLength*kLength*kLength;
							threeVector<int> kInd;
							for(int dim=0;dim<3;dim++)
							for(kInd.s[dim]=-kLength;kInd.s[dim]<=kLength;kInd.s[dim]+=kLength)
							for(kInd.s[(dim+1)%3]=-kLength+(dim>0);kInd.s[(dim+1)%3]<=kLength-(dim>0);kInd.s[(dim+1)%3]++)
							for(kInd.s[(dim+2)%3]=-kLength+(dim>0);kInd.s[(dim+2)%3]<=kLength-(dim>0);kInd.s[(dim+2)%3]++)
							{
								threeVector<int> nCoor=iCoor;
								nCoor.x+=(kInd.x%kLength-checkLayer);
								nCoor.y+=(int(kInd.y/kLength)%kLength-checkLayer);
								nCoor.z+=(int(kInd.z/(kLength*kLength))-checkLayer);
								
								int neighborHash=hash(nCoor,vertN);
								
								auto it=neighbors.find(neighborHash);
								if(it!=neighbors.end() && it->first==neighborHash && neighborHash!=currentHash)// && it->second.type!=corePoint)
								{
									for(int triangleIndex=it->second.size()-1;triangleIndex>=0;--triangleIndex)
									{
										if((triangles[it->second[triangleIndex]].x->first==currentHash ||
											triangles[it->second[triangleIndex]].y->first==currentHash ||
											triangles[it->second[triangleIndex]].z->first==currentHash))// && 
											//!delTriangle[it->second[triangleIndex]])
										{
											tCandidates.push_back(triangles[it->second[triangleIndex]]);
											//delTriangle[it->second[triangleIndex]]=true;
											//it->second[vertex]=tCandidates.back();
											//it->second.pop_back();
											nonAdded=false;
										}
									}
								}
							}
							checkLayer++;
						} while (!nonAdded);
						*/
						
						
						std::vector<int> validRotations;
						//it is inside?
						if(0<iCoor.x<vertN.x-1 && 0<iCoor.y<vertN.y-1 && 0<iCoor.z<vertN.z-1)
						{
							//make sure any path is complete
							//
							//handle incomplete edges by just filling them in
							std::vector< std::vector<int> > vertexPath(tCandidates.size(),std::vector<int>());
							//all corners, looking for shared edges
							
							for(int aT=0;aT<tCandidates.size();aT++)
							{
								const auto& aTriangle=tCandidates[aT];
								for(int bT=aT+1;bT<tCandidates.size();bT++)
								{
									const auto& bTriangle=tCandidates[bT];
									int aIOrient,bIOrient;
									//Find the vertex matching the current index
									for(int iVert=0;iVert<3;iVert++)
									{
										if(aTriangle.s[iVert]->first==currentHash)
											aIOrient=iVert;
										if(bTriangle.s[iVert]->first==currentHash)
											bIOrient=iVert;
									}
									//Found an orientation index, all triangles are oriented as ijk and ik'j for neighbors
									//looking for j
									int aJOrient=-1;
									int bJOrient=-1;
									for(int jVert=1;jVert<3;jVert++)
									{
										aJOrient=(jVert+aIOrient)%3;
										bJOrient=(bIOrient-jVert+3)%3;
										
										//std::cout << aJOrient << ' ' << bJOrient << std::endl;
										//std::cin.get();
										//they don't match, default to -1
										if(aTriangle.s[aJOrient]->first!=bTriangle.s[bJOrient]->first)
										{
											aJOrient=-1;
											bJOrient=-1;
										}
										else
										{
											break;
										}
									}
									//they match
									if(aJOrient!=-1 && bJOrient!=-1)
									{
										int aKOrient, bKOrient;
										//nEdgeMatches++;
										//aKOrient=(aIOrient+1)%3;
										//if(aKOrient==aJOrient)
										//	aKOrient=(aKOrient+1)%3;
										//bKOrient=(bIOrient+1)%3;
										//if(bKOrient==bJOrient)
										//	bKOrient=(bKOrient+1)%3;
										//check:
										//    i
										//  k<|>k'
										//    j
										// for consistency
										
										//position<double> aIPos=aTriangle.s[aIOrient]->second;
										//position<double> aJPos=aTriangle.s[aJOrient]->second;
										//position<double> aKPos=aTriangle.s[aKOrient]->second;
										//position<double> bKPos=bTriangle.s[bKOrient]->second;
										/*
										dab.x=aP.x-bP.x;
										dab.y=aP.y-bP.y;
										dab.z=aP.z-bP.z;
										while(dab.x>size.x/2.0)dab.x-=size.x;
										while(dab.y>size.y/2.0)dab.y-=size.y;
										while(dab.z>size.z/2.0)dab.z-=size.z;
										while(dab.x<-size.x/2.0)dab.x+=size.x;
										while(dab.y<-size.y/2.0)dab.y+=size.y;
										while(dab.z<-size.z/2.0)dab.z+=size.z;
										//bc
										dbc.x=bP.x-cP.x;
										dbc.y=bP.y-cP.y;
										dbc.z=bP.z-cP.z;
										while(dbc.x>size.x/2.0)dbc.x-=size.x;
										while(dbc.y>size.y/2.0)dbc.y-=size.y;
										while(dbc.z>size.z/2.0)dbc.z-=size.z;
										while(dbc.x<-size.x/2.0)dbc.x+=size.x;
										while(dbc.y<-size.y/2.0)dbc.y+=size.y;
										while(dbc.z<-size.z/2.0)dbc.z+=size.z;
										//ca
										dca.x=cP.x-aP.x;
										dca.y=cP.y-aP.y;
										dca.z=cP.z-aP.z;
										while(dca.x>size.x/2.0)dca.x-=size.x;
										while(dca.y>size.y/2.0)dca.y-=size.y;
										while(dca.z>size.z/2.0)dca.z-=size.z;
										while(dca.x<-size.x/2.0)dca.x+=size.x;
										while(dca.y<-size.y/2.0)dca.y+=size.y;
										while(dca.z<-size.z/2.0)dca.z+=size.z;
										*/
										vertexPath[aT].push_back(bT);
										vertexPath[bT].push_back(aT);
//std::cout << aT << ' ' << bT << std::endl;
									}
								}
							}
							
							for(int vertex=vertexPath.size()-1;vertex>=0;--vertex)
							{
								//incomplete connection
								if(vertexPath[vertex].size()==0 || vertexPath[vertex].size()==1)
								{
									tCandidates[vertex]=tCandidates.back();
									tCandidates.pop_back();
								}
								else
								{
									std::vector<int> stackTraverse;
									stackTraverse.push_back(vertex);
									
									while(stackTraverse.size()>0)
									{
										int vertexNumber=stack.back();
										
									}
								}
								std::cout << vertex;
								for(auto& vertVal:vertexPath[vertex])
									std::cout << ' ' << vertVal;
								std::cout << std::endl;
							}
						}
						
						auto &nVertex=neighbors[currentHash];
						for(auto& tCand:tCandidates)
						{
							stack.push_back(tCand.y);
							stack.push_back(tCand.z);
							tCand.x->second.type=edgePoint;
							tCand.y->second.type=edgePoint;
							tCand.z->second.type=edgePoint;
							nVertex.push_back(triangles.size());
							//delTriangle.push_back(false);
							triangles.push_back(tCand);
						}
							
						throw 0;
					}
				}
			}
//throw 0;
			//sTriangles.push_back(triangles);
			vNeighbors.push_back(neighbors);
			//sNormals.push_back(normals);
		}
		
		//int nTriangles=0;
		//for(auto& triangles:sTriangles)
		//	nTriangles+=triangles.size();	
		std::cerr << "Found " << surfaces.size() << " surfaces with " << triangles.size() << " potential triangles." << std::endl;
		
		std::vector<double> normalsDist;
		int nNormalsIndex=100;
		int surfaceNumber=0;
		
		std::vector<bool > traversed(triangles.size(),false);
		
		//std::vector<double> totalCurvatureE;
		int s=0;
		//all surfaces
		for(auto& neighbors:vNeighbors)
		{
			double curvatureE=0;
			double tCurvature=0;
			
			//all vertices
			for(auto& nVertex:neighbors)
			{
				//handle incomplete edges by just filling them in
				std::vector<int> incompleteEdges;
				//all corners, looking for shared edges
				for(int i=0;i<nVertex.second.size();i++)
				{
					const auto& aTriangle=triangles[nVertex.second[i]];
					int nEdgeMatches=0;//this should be nVertex.second.size()+1 on a flat surface
					threeVector<double> normalSqrSum=0;
					double areaSum=0;
					for(int j=i+1;j<nVertex.second.size();j++)
					{
						const auto& bTriangle=triangles[nVertex.second[j]];
						int aIOrient,bIOrient;
						//Find the vertex matching the current index
						for(int iVert=0;iVert<3;iVert++)
						{
							if(aTriangle.s[iVert]->first==nVertex.first)
								aIOrient=iVert;
							if(bTriangle.s[iVert]->first==nVertex.first)
								bIOrient=iVert;
						}
						//Found an orientation index, all triangles are oriented as ijk and ik'j for neighbors
						//looking for j
						int aJOrient=-1;
						int bJOrient=-1;
						for(int jVert=1;jVert<3;jVert++)
						{
							aJOrient=(jVert+aIOrient)%3;
							bJOrient=(bIOrient-jVert+3)%3;
							
							//std::cout << aJOrient << ' ' << bJOrient << std::endl;
							//std::cin.get();
							//they don't match, default to -1
							if(aTriangle.s[aJOrient]->first!=bTriangle.s[bJOrient]->first)
							{
								aJOrient=-1;
								bJOrient=-1;
							}
							else
							{
								break;
							}
						}
						//they match, compute the curvature energy
						if(aJOrient!=-1 && bJOrient!=-1)
						{
							int aKOrient, bKOrient;
							nEdgeMatches++;
							aKOrient=(aIOrient+1)%3;
							if(aKOrient==aJOrient)
								aKOrient=(aKOrient+1)%3;
							bKOrient=(bIOrient+1)%3;
							if(bKOrient==bJOrient)
								bKOrient=(bKOrient+1)%3;
							
							//     a
							//   b<|>d
							//     c
							//
							
							threeVector<double> dab,dcb,dac, dad,dcd;
							position<double> aP=aTriangle.s[aIOrient]->second;
							position<double> bP=aTriangle.s[aKOrient]->second;
							position<double> cP=aTriangle.s[aJOrient]->second;
							
							threeVector<double> aFaceCenter=0, bFaceCenter=0;
							aFaceCenter.x=aP.x;
							aFaceCenter.y=aP.y;
							aFaceCenter.z=aP.z;
							aFaceCenter.x+=bP.x;
							aFaceCenter.y+=bP.y;
							aFaceCenter.z+=bP.z;
							aFaceCenter.x+=cP.x;
							aFaceCenter.y+=cP.y;
							aFaceCenter.z+=cP.z;
							aFaceCenter.x/=3.0;
							aFaceCenter.y/=3.0;
							aFaceCenter.z/=3.0;
							
							dab.x=aP.x-bP.x;
							dab.y=aP.y-bP.y;
							dab.z=aP.z-bP.z;
							while(dab.x>size.x/2.0)dab.x-=size.x;
							while(dab.y>size.y/2.0)dab.y-=size.y;
							while(dab.z>size.z/2.0)dab.z-=size.z;
							while(dab.x<-size.x/2.0)dab.x+=size.x;
							while(dab.y<-size.y/2.0)dab.y+=size.y;
							while(dab.z<-size.z/2.0)dab.z+=size.z;
							
							dcb.x=cP.x-bP.x;
							dcb.y=cP.y-bP.y;
							dcb.z=cP.z-bP.z;
							while(dcb.x>size.x/2.0)dcb.x-=size.x;
							while(dcb.y>size.y/2.0)dcb.y-=size.y;
							while(dcb.z>size.z/2.0)dcb.z-=size.z;
							while(dcb.x<-size.x/2.0)dcb.x+=size.x;
							while(dcb.y<-size.y/2.0)dcb.y+=size.y;
							while(dcb.z<-size.z/2.0)dcb.z+=size.z;
							
							//dab=unitVector(dab);
							//dcb=unitVector(dcb);
							
							//double acos=dotProduct(dab,dcb);
							//double cotSum=sqrt(acos/(1-acos));
							//double cotSum
							
							//flip this vector, dab->dba
//std::cout << acos << ' ' << cotSum << std::endl;
//std::cin.get();
							dac.x=aP.x-cP.x;
							dac.y=aP.y-cP.y;
							dac.z=aP.z-cP.z;
							while(dac.x>size.x/2.0)dac.x-=size.x;
							while(dac.y>size.y/2.0)dac.y-=size.y;
							while(dac.z>size.z/2.0)dac.z-=size.z;
							while(dac.x<-size.x/2.0)dac.x+=size.x;
							while(dac.y<-size.y/2.0)dac.y+=size.y;
							while(dac.z<-size.z/2.0)dac.z+=size.z;
							double lac=sqrt(dotProduct(dac,dac));
							
							
							double aArea=sqrt(dotProduct(dab,dab)*dotProduct(dcb,dcb)-pow(dotProduct(dab,dcb),2.0))/2.0;
							dab=unitVector(dab);
							dcb=unitVector(dcb);
							
							threeVector<double> aNormal=unitVector(crossProduct(dab,dcb));
							
							//dab.x=-dab.x;
							//dab.y=-dab.y;
							//dab.z=-dab.z;
							//double aArea=sqrt(dotProduct(dab,dab)*dotProduct(dbc,dbc)-pow(dotProduct(dab,dbc),2.0))/2.0;
							
							//flip it back
							//dab.x=-dab.x;
							//dab.y=-dab.y;
							//dab.z=-dab.z;
							//areaSum+=sqrt(dotProduct(dab,dab)*dotProduct(dbc,dbc)-pow(dotProduct(dab,dbc),2.0))/2.0;
							//aP=bTriangle.x.second;
							//bP=bTriangle.y.second;
							bP=bTriangle.s[bKOrient]->second;
							
							bFaceCenter.x=aP.x;
							bFaceCenter.y=aP.y;
							bFaceCenter.z=aP.z;
							bFaceCenter.x+=bP.x;
							bFaceCenter.y+=bP.y;
							bFaceCenter.z+=bP.z;
							bFaceCenter.x+=cP.x;
							bFaceCenter.y+=cP.y;
							bFaceCenter.z+=cP.z;
							bFaceCenter.x/=3.0;
							bFaceCenter.y/=3.0;
							bFaceCenter.z/=3.0;
							
							//dab.x=aP.x-bP.x;
							//dab.y=aP.y-bP.y;
							//dab.z=aP.z-bP.z;
							//while(dab.x>size.x/2.0)dab.x-=size.x;
							//while(dab.y>size.y/2.0)dab.y-=size.y;
							//while(dab.z>size.z/2.0)dab.z-=size.z;
							//while(dab.x<-size.x/2.0)dab.x+=size.x;
							//while(dab.y<-size.y/2.0)dab.y+=size.y;
							//while(dab.z<-size.z/2.0)dab.z+=size.z;
							//bc
							dad.x=aP.x-bP.x;
							dad.y=aP.y-bP.y;
							dad.z=aP.z-bP.z;
							while(dad.x>size.x/2.0)dad.x-=size.x;
							while(dad.y>size.y/2.0)dad.y-=size.y;
							while(dad.z>size.z/2.0)dad.z-=size.z;
							while(dad.x<-size.x/2.0)dad.x+=size.x;
							while(dad.y<-size.y/2.0)dad.y+=size.y;
							while(dad.z<-size.z/2.0)dad.z+=size.z;
							
							dcd.x=cP.x-bP.x;
							dcd.y=cP.y-bP.y;
							dcd.z=cP.z-bP.z;
							while(dcd.x>size.x/2.0)dcd.x-=size.x;
							while(dcd.y>size.y/2.0)dcd.y-=size.y;
							while(dcd.z>size.z/2.0)dcd.z-=size.z;
							while(dcd.x<-size.x/2.0)dcd.x+=size.x;
							while(dcd.y<-size.y/2.0)dcd.y+=size.y;
							while(dcd.z<-size.z/2.0)dcd.z+=size.z;
							
							
							
//std::cout << aArea << ' ' << bArea << ' ' << static_cast<double>(magnitude(dad)) << ' ' << static_cast<double>(magnitude(dcd)) << std::endl;
//std::cin.get();
							
							dad=unitVector(dad);
							dcd=unitVector(dcd);
							
							double bArea=sqrt(dotProduct(dad,dad)*dotProduct(dcd,dcd)-pow(dotProduct(dad,dcd),2.0))/2.0;
							
							//acos=dotProduct(dad,dcd);
							//cotSum+=sqrt(acos/(1-acos));
							//cotSum*=lac/2.0;
							
							
							threeVector<double> bNormal=unitVector(crossProduct(dad,dcd));
							if(dotProduct(aNormal,bNormal)<0)
							{
								bNormal.x=-bNormal.x;
								bNormal.y=-bNormal.y;
								bNormal.z=-bNormal.z;
							}
if(nVertex.second[i]==1000)//!traversed[nVertex.second[i]])
{
int base=1;
//for(double xyz=0;xyz<1.0;xyz+=0.1)
//	std::cout << base++ << ' ' << (2.0*aFaceCenter.x+aNormal.x)*xyz/10.0+aFaceCenter.x << 
//	' ' << (2.0*aFaceCenter.y+aNormal.y)*xyz/10.0+aFaceCenter.y << 
//	' ' << (2.0*aFaceCenter.z+aNormal.z)*xyz/10.0+aFaceCenter.z << '\n';
	
	std::cout << base++ << ' ' << aFaceCenter.x << 
	' ' << aFaceCenter.y << 
	' ' << aFaceCenter.z << '\n';
	std::cout << base++ << ' ' << aNormal.x+aFaceCenter.x << 
	' ' << aNormal.y+aFaceCenter.y << 
	' ' << aNormal.z+aFaceCenter.z << '\n';
traversed[nVertex.second[i]]=true;
}
							
							//threeVector<double> dNormal=aNormal-bNormal;
							
							//double bArea=sqrt(dotProduct(dab,dab)*dotProduct(dbc,dbc)-pow(dotProduct(dab,dbc),2.0))/2.0;
							
							//using |ni-nj|^2
							//curvatureE+=dotProduct(dNormal,dNormal)*sqrt(dotProduct(dab,dab))/(4.0*(aArea+bArea));
							//using (1-ni.nj)
							//curvatureE+=(1.0-dotProduct(aNormal,bNormal))*sqrt(dotProduct(dab,dab))/(4.0*(aArea+bArea));
							//using ni.nj
							//double buffe=dotProduct(aNormal,bNormal);///(4.0*(aArea+bArea));
							tCurvature+=dotProduct(aNormal,bNormal);//dotProduct(aNormal,bNormal)*sqrt(dotProduct(dab,dab))/(4.0*(aArea+bArea));
							
							//Itzykson 1986
							normalSqrSum.x+=dac.x;
							normalSqrSum.y+=dac.y;
							normalSqrSum.z+=dac.z;
							areaSum+=aArea+bArea;
							
							int histInd=dotProduct(aNormal,bNormal)*nNormalsIndex;
							if(surfaces.size()<4 && surfaces[surfaceNumber].size()>10)
							{
								while(histInd>=normalsDist.size())
									normalsDist.push_back(0.0);
								normalsDist[histInd]++;
							}
						}
					}
					//first triangle is special, 2 matches can be found
					//if(i==0 && nEdgeMatches<2)
					//	incompleteEdges.push_back(i);
					//any other triangle should only have 1
					//if(i>0 && nEdgeMatches<1)
					//	incompleteEdges.push_back(i);
					//unusual cases
					//if((i==0 && nEdges>2) || (i>0 && nEdges>1))
					{
						
					}
					//curvatureE+=normalSum/area
					
					if(nEdgeMatches>0)
					{
						areaSum/=static_cast<double>(nEdgeMatches);
						curvatureE+=dotProduct(normalSqrSum,normalSqrSum)/areaSum;
					}
				}
				//some edges are broken, try to fix them
				//for(auto& iFix:incompleteEdges)
				{
					//walk around vertex
					
				}
				
			}
			if(surfaces.size()<4)
			{
				std::cerr << "Curvature Energy: " << curvatureE << std::endl;
			std::string buf;
			buf="curvatureE_";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << surfaceNumber;
			parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".dat";
			std::fstream dataFile;
			dataFile.open(buf.c_str(), std::ios::out | std::ios::app);
			
			dataFile << time << ' ' << curvatureE << ' ' << tCurvature << std::endl;
			dataFile.close();
			
			}
			surfaceNumber++;
		}
		
		/*
		for(int i=0;i<normalsDist.size() && surfaces.size()<4;i++)
		{
			std::cout << static_cast<double>(i)/static_cast<double>(nNormalsIndex) << ' ' << normalsDist[i] << '\n';
		}
		if(surfaces.size()<4)
		{
			std::cout << std::endl;
		}
		*/
		//std::cout << time << '\t';
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
		frame++;
		throw 0;
	}
	return 0;
}
