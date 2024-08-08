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

template <typename T>
struct threeDHash {
		int operator()(T input)
		{
			threeVector<int> iCoor;
			
			iCoor.x=(static_cast<int>(input.x/cellSize.x)+n.x)%n.x;
			iCoor.y=(static_cast<int>(input.y/cellSize.y)+n.y)%n.y;
			iCoor.z=(static_cast<int>(input.z/cellSize.z)+n.z)%n.z;
			
			int hash=iCoor.x+iCoor.y*n.x+iCoor.z*n.x*n.y;
			if(hash<0)
			{
				std::cerr << "Error: " << hash << std::endl;
			}
			
			return iCoor.x+iCoor.y*n.x+iCoor.z*n.x*n.y;
		};
		
		int operator()(int a, int b, int c)
		{
			threeVector<int> iCoor;
			
			iCoor.x=(a+10*n.x)%n.x;
			iCoor.y=(b+10*n.y)%n.y;
			iCoor.z=(c+10*n.z)%n.z;
			
			int hash=iCoor.x+iCoor.y*n.x+iCoor.z*n.x*n.y;
			if(hash<0)
			{
				std::cerr << "Error: " << hash << std::endl;
			}
			
			return iCoor.x+iCoor.y*n.x+iCoor.z*n.x*n.y;
		};
		
		threeVector<int> cellCoord(T input)
		{
			threeVector<int> iCoor;
			
			iCoor.x=(static_cast<int>(input.x/cellSize.x)+n.x)%n.x;
			iCoor.y=(static_cast<int>(input.y/cellSize.y)+n.y)%n.y;
			iCoor.z=(static_cast<int>(input.z/cellSize.z)+n.z)%n.z;
			
			return iCoor;
		};
		
		threeVector<int> cellCoord(int a, int b, int c)
		{
			threeVector<int> iCoor;
			
			iCoor.x=(a+n.x)%n.x;
			iCoor.y=(b+n.y)%n.y;
			iCoor.z=(c+n.z)%n.z;
			
			return iCoor;
		};
		
		threeVector<int> cellCoord(int input)
		{
			threeVector<int> iCoor;
			
			iCoor.x=input%n.x;
			iCoor.y=(input/n.x)%n.y;
			iCoor.z=(input/(n.x*n.y))%n.z;
			
			return iCoor;
		}
		
		T cellSize;
		threeVector<int> n;
};

int main(int argc, char **argv)
{
	if(argc<7)
	{
		std::cerr << argv[0] << " name nX nY surfSpacing dry type1 type2 ...\n";
		std::cerr << "\tdry is for no output, if the surface spacing needs to be optimized. 0=false 1=true\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	threeVector<int> n;
	double surfSpacing;
	int dry;
	cmdArg >> n.x >> n.y >> surfSpacing >> dry;
	double surfSpacingSqr=surfSpacing*surfSpacing;
	
	std::vector<int> types;
	
	for(int i=4;i<argc;i++)
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
	
	while(xyzFile.load())
	{
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
		
		//threeDHash< position<double> > cellHash;
		
		//cellHash.n.x=static_cast<int>(size.x/surfSpacing);
		//cellHash.n.y=static_cast<int>(size.y/surfSpacing);
		//cellHash.n.z=static_cast<int>(size.z/surfSpacing);
		
		//cellHash.cellSize.x=size.x/static_cast<double>(cellHash.n.x);
		//cellHash.cellSize.y=size.y/static_cast<double>(cellHash.n.y);
		//cellHash.cellSize.z=size.z/static_cast<double>(cellHash.n.z);
		
		typedef std::vector<int> cellType;//fastest is vector, deque is a little slower
		
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
				buf.push_back(j);
				//cells[hash]=buf;
				cells.insert(it, std::map<int,cellType >::value_type(hash,buf));
			}
		}
		/*
		int color=0;
		for(std::map<int,cellType >::iterator current=cells.begin();current!=cells.end();++current)
		{
			color++;
			for(int k=0;current!=cells.end() && k<current->second.size();k++)
			{
				int l=pIndex[current->second[k]];
				std::cout << color%5 << '\t' << p[l].x << '\t' << p[l].y << '\t' << p[l].z << std::endl;
			}
		}
		
		throw 0;
		*/
		std::cerr << pIndex.size() << '\t';
		
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		
		
		//get surfaces
		std::vector< std::vector<int> > surfaces;//this just grows with indices per surface
		std::vector<int> sFlag(pIndex.size(),-1);//cooresponds to the surfaces found i.e. 0=surface 1, 1=surface 2...
		for(int i=0;i<pIndex.size();++i)
		{
			//unaccounted surface
			if(sFlag[i]==-1)
			{
				std::vector<int> surface;//new surface
				//std::vector<int> stack;
				std::deque<int> stack;//almost the same as vector
				//std::list<int> stack;
				stack.push_back(i);
				sFlag[i]=surfaces.size();
				surface.push_back(pIndex[i]);
				while(stack.size()>0)
				{
					//std::cerr << stack.size() << std::endl;
					
					int element=stack.back();
					stack.pop_back();
					
					int j=pIndex[element];
					
					//std::cerr << stack.size() << '\t' << surface.size() << '\t' << surfaces.size() << std::endl;
					
					//calculate hash
					threeVector<int> iCoor;
					threeVector<double> corP;
					for(corP.x=p[j].x;corP.x>size.x;corP.x-=size.x);
					for(corP.y=p[j].y;corP.y>size.y;corP.y-=size.y);
					for(corP.z=p[j].z;corP.z>size.z;corP.z-=size.z);
					for(corP.x=corP.x;corP.x<0;corP.x+=size.x);
					for(corP.y=corP.y;corP.y<0;corP.y+=size.y);
					for(corP.z=corP.z;corP.z<0;corP.z+=size.z);
					iCoor.x=static_cast<int>(corP.x/surfR.x)%surfN.x;
					iCoor.y=static_cast<int>(corP.y/surfR.y)%surfN.y;
					iCoor.z=static_cast<int>(corP.z/surfR.z)%surfN.z;
					//int hash=iCoor.x+iCoor.y*surfN.x;
					//threeVector<int> iCoor=cellHash.cellCoord(p[j]);
					//find connected elements and push onto stack
					//check the nearby cells for traversal
					for(int k=0;k<27;k++)
					{
						int a=k%3-1;
						int b=int(k/3)%3-1;
						int c=int(k/9)-1;
						int neighborHash=
						((iCoor.x+a+surfN.x)%surfN.x)+
						(((iCoor.y+b+surfN.y)%surfN.y)*surfN.x)+
						(((iCoor.z+c+surfN.z)%surfN.z)*surfN.x*surfN.y);
						//cellHash(iCoor.x+a,iCoor.y+b,iCoor.z+c);
						//std::cerr << hash << '\t' << neighborHash << std::endl;
						
						std::map<int,cellType >::iterator it=cells.find(neighborHash);
						for(int m=0;it!=cells.end() && m<it->second.size();m++)
							//for(int m=0;m<pIndex.size();m++)
						{
							threeVector<double> d;
							
							int l=pIndex[it->second[m]];
							//int l=pIndex[m];
							d.x=p[l].x-p[j].x;
							d.y=p[l].y-p[j].y;
							d.z=p[l].z-p[j].z;
							
							while(d.x>size.x/2.0)d.x-=size.x;
							while(d.y>size.y/2.0)d.y-=size.y;
							while(d.z>size.z/2.0)d.z-=size.z;
							while(d.x<-size.x/2.0)d.x+=size.x;
							while(d.y<-size.y/2.0)d.y+=size.y;
							while(d.z<-size.z/2.0)d.z+=size.z;
							
							if(l!=j && sFlag[it->second[m]]==-1 && d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
								//if(l!=j && sFlag[m]==-1 && d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
							{
								sFlag[it->second[m]]=surfaces.size();
								surface.push_back(l);
								
								stack.push_back(it->second[m]);
								
								do
								{
									d.z=p[l].z-p[j].z;
									if(d.z>size.z/2.0)
										p[l].z-=size.z;
									d.z=p[l].z-p[j].z;
								} while(d.z>size.z/2.0);
								
								do
								{
									d.z=p[l].z-p[j].z;
									if(d.z<-size.z/2.0)
										p[l].z+=size.z;
									d.z=p[l].z-p[j].z;
								} while(d.z<-size.z/2.0);
								
								
								//std::cout << d.x*d.x+d.y*d.y+d.z*d.z << '\t' << hash << '\t' << neighborHash << std::endl;
								
								//sFlag[m]=surfaces.size();
								//stack.push_back(m);
							}
						}
					}
				}
				
				//for(int l=0;l<surface.size();l++)
				//	std::cout << 1 << '\t' << p[surface[l]].x << '\t' << p[surface[l]].y << '\t' << p[surface[l]].z << std::endl;
				
				//throw 0;
				//std::cerr << surfaces.size() << std::endl;
				surfaces.push_back(surface);
			}
		}
		
		//our first surface is assumed to be correct
		if(sFlagOrigin.size()==0)
		{
			/*
			std::cout << surfaces.size() << std::endl;
			for(int i=0;i<surfaces.size();i++)
			{
				std::cout << surfaces[i].size() << '\t' << pSig*pIndex.size() << std::endl;
				if(surfaces[i].size()<pSig*pIndex.size())
				{
					int candidate=0;
					double candidateDistance=size.x*size.x+size.y*size.y+size.z*size.z;
					for(int j=0;j<surfaces[i].size();j++)
					{
						int m=pIndex[surfaces[i][j]];
						for(int k=0;k<pIndex.size();k++)
						{
							int l=pIndex[k];
							if(sFlag[k]!=sFlag[j])
							{
								
								threeVector<double> d;
								
								d.x=p[l].x-p[m].x;
								d.y=p[l].y-p[m].y;
								d.z=p[l].z-p[m].z;
								
								while(d.x>size.x/2.0)d.x-=size.x;
								while(d.y>size.y/2.0)d.y-=size.y;
								while(d.z>size.z/2.0)d.z-=size.z;
								while(d.x<-size.x/2.0)d.x+=size.x;
								while(d.y<-size.y/2.0)d.y+=size.y;
								while(d.z<-size.z/2.0)d.z+=size.z;
								
								double dr=d.x*d.x+d.y*d.y+d.z*d.z;
								
								if(dr<candidateDistance)
								{
									std::cout << "lalala " << m << '\t' << l << std::endl;
									candidateDistance=dr;
									candidate=k;
								}
							}
						}
					}
					std::cout << "asdfasdf " << sFlag[candidate] << '\t' <<  candidate << '\t' << i << std::endl;
					while(surfaces[i].size()>0)
					{
						std::cout << "removing " << surfaces[i].size() << " where " << sFlag[pIndex[surfaces[i].back()]] << " and " << sFlag[candidate] << std::endl;
						surfaces[sFlag[candidate]].push_back(surfaces[i].back());
						sFlag[surfaces[i].back()]=sFlag[candidate];
						surfaces[i].pop_back();
						std::cout << "now " << surfaces[i].size() << std::endl;
					}
				}
			}
			
			for(int i=0;i<surfaces.size();i++)
			{
				if(surfaces[i].size()==0)
				{
					surfaces.erase(surfaces.begin()+i);
					i--;
				}
			}
			*/
			sFlagOrigin=sFlag;
			nSurfaces=surfaces.size();
			//throw 0;
		}
		else//make new surfaces match
		{
			std::cout << std::endl;
			std::vector< std::vector<int> > cMap(nSurfaces, std::vector<int>(surfaces.size(),0) );
			for(int i=0;i<sFlag.size();i++)
			{
				cMap[sFlagOrigin[i]][sFlag[i]]++;
			}
			std::vector<int> sMap(nSurfaces, 0);
			for(int i=0;i<cMap.size();i++)
				for(int j=0;j<cMap[i].size();j++)
					sMap[i]+=cMap[i][j];
			
			std::cout << "nP\t\t";
			for(int i=0;i<surfaces.size();i++)
				std::cout << surfaces[i].size() << '\t';
			std::cout << "\n\to/n\t";
			for(int i=0;i<surfaces.size();i++)
				std::cout << i << '\t';
			std::cout << std::endl;
			for(int i=0;i<cMap.size();i++)
			{
				std::cout << sMap[i] << '\t' << i << "\t[";
				for(int j=0;j<cMap[i].size();j++)
				{
					std::cout << std::fixed << std::setprecision(4) << static_cast<double>(cMap[i][j])/static_cast<double>(sMap[i]) << '\t';
				}
				std::cout << "]" << std::endl;
			}
			//no match, swap with other
			if(cMap[0][0]<pSig*pIndex.size())
			{
				std::cerr << "\nSwapping: " << cMap[0][0] << '\t' << pSig*pIndex.size() << std::endl;
				for(int i=0;i<cMap[0].size();i++)
				{
					if(cMap[0][i]>pSig*pIndex.size())
					{
						surfaces[i].swap(surfaces[0]);
					}
				}
			}
		}
		
		//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		//std::cerr << time_span.count() << std::endl;
		
		std::cerr << surfaces.size() << '\t';
		
		//make the height map
		threeVector<double> r;
		r.x=size.x/static_cast<double>(n.x);
		r.y=size.y/static_cast<double>(n.y);
		
		//std::cout << pIndex.size() << std::endl << "test" << std::endl;
		
		//for(int s=0;s<surfaces.size();s++)
		//	for(int j=0;j<surfaces[s].size();j++)
		//		std::cout << s << '\t' << p[surfaces[s][j]].x << '\t' << p[surfaces[s][j]].y << '\t' << p[surfaces[s][j]].z << std::endl;
		
		for(int s=0;s<surfaces.size() && dry==0 && surfaces.size()==2;s++)
		{
			std::cout << time << '\t';
			//std::vector<double> heightMap(n.x*n.y, 0.0);
			std::vector< std::vector<double> > heightMap(n.x,std::vector<double>(n.y,0.0));
			std::vector< std::vector<double> > gradSquaredMap(n.x,std::vector<double>(n.y,0.0));
			std::vector< std::vector<int> > countMap(n.x,std::vector<int>(n.y,0));
			//std::vector<int> countMap(n.x*n.y, 0.0);
			
			for(int j=0;j<surfaces[s].size();j++)
			{
				int i=surfaces[s][j];
				threeVector<int> c;
				c.x=static_cast<int>(p[i].x/r.x);
				c.y=static_cast<int>(p[i].y/r.y);
				
				c.x=(c.x>=n.x)?n.x-1:c.x;
				c.y=(c.y>=n.y)?n.y-1:c.y;
				/*
				int hash=c.x+c.y*n.x;
				if(hash>=heightMap.size())
				{
					std::cerr << "Out of bounds particle!\n";
					std::cerr << "Size: " << size.x << '\t' << size.y << std::endl;
					std::cerr << "Position: " << p[i].x << '\t' << p[i].y << std::endl;
				}
				
				heightMap[hash]+=p[i].z;
				countMap[hash]++;
				*/
				heightMap[c.x][c.y]+=p[i].z;
				
				countMap[c.x][c.y]++;
			}
			
			for(int j=0;j<heightMap.size();j++)
				for(int k=0;k<heightMap[j].size();k++)
					if(countMap[j][k]!=0)
						heightMap[j][k]/=static_cast<double>(countMap[j][k]);
			for(int j=0;j<heightMap.size();j++)
				for(int k=0;k<heightMap[j].size();k++)
					if(heightMap[j][k]==0)
						heightMap[j][k]=(heightMap[(j-1+n.x)%n.x][k]+
										heightMap[(j+1+n.x)%n.x][k]+
										heightMap[j][(k-1+n.y)%n.y]+
										heightMap[j][(k+1+n.y)%n.y])/4.0;
			/*
			for(int j=0;j<heightMap.size();j++)
			{
				for(int k=0;k<heightMap[j].size();k++)
				{
					heightMap[j][k]+=size.z/3.0;
					heightMap[j][k]-=(heightMap[j][k]>size.z)?size.z:0.0;
				}
			}
			*/		
			double gradSquared=0;
			twoVector<double> aSum,bSum,hSum;
			aSum.x=0;
			aSum.y=0;
			bSum.x=0;
			bSum.y=0;
			hSum.x=0;
			hSum.y=0;
			
			//integrate gradSquared height over map
			for(int j=0;j<heightMap.size();j++)
			{
				for(int k=0;k<heightMap[j].size();k++)
				{
					double area=r.x*r.y;//crude
					
					
					twoVector<double> da,db,dc;
					//(1,0)
					da.x=heightMap[(j+1+n.x)%n.x][k]-heightMap[(j-1+n.x)%n.x][k];
					//(0,1)
					da.y=heightMap[j][(k+1+n.y)%n.y]-heightMap[j][(k-1+n.y)%n.y];
					
					da.x-=(da.x>size.z/2.0)?size.z:0;
					da.y-=(da.y>size.z/2.0)?size.z:0;
					da.x+=(da.x<-size.z/2.0)?size.z:0;
					da.y+=(da.y<-size.z/2.0)?size.z:0;
					
					area=sqrt(da.x*da.x+r.x*r.x)*sqrt(da.y*da.y+r.y*r.y);//corrected
					
					//(-1,0)
					//db.x=-(heightMap[j][k]-heightMap[(j-1+n.x)%n.x][k]);
					//(0,-1)
					//db.y=-(heightMap[j][k]-heightMap[j][(k-1+n.y)%n.y]);
					
					//db.x-=(db.x>size.z/2.0)?size.z:0;
					//db.y-=(db.y>size.z/2.0)?size.z:0;
					//db.x+=(db.x<-size.z/2.0)?size.z:0;
					//db.y+=(db.y<-size.z/2.0)?size.z:0;
					
					//(-1,1)
					//dc.x=heightMap[(j+1+n.x)%n.x][k]-heightMap[(j-1+n.x)%n.x][k];
					//(1,-1)
					//dc.y=heightMap[j][(k+1+n.y)%n.y]-heightMap[j][(k-1+n.y)%n.y];
					
					//dc.x-=(dc.x>size.z/2.0)?size.z:0;
					//dc.y-=(dc.y>size.z/2.0)?size.z:0;
					//dc.x+=(dc.x<-size.z/2.0)?size.z:0;
					//dc.y+=(dc.y<-size.z/2.0)?size.z:0;
					
					gradSquaredMap[j][k]=pow(((heightMap[(j+1+n.x)%n.x][k]+heightMap[(j-1+n.x)%n.x][k]+
										heightMap[j][(k-1+n.y)%n.y]+heightMap[j][(k+1+n.y)%n.y])/2.0+
										(heightMap[(j-1+n.x)%n.x][(k-1+n.y)%n.y]+heightMap[(j-1+n.x)%n.x][(k+1+n.y)%n.y]+
										heightMap[(j+1+n.x)%n.x][(k-1+n.y)%n.y]+heightMap[(j+1+n.x)%n.x][(k+1+n.y)%n.y])/4.0-
										3.0*heightMap[j][k]),2.0)/(area);
					
					gradSquared+=pow(((heightMap[(j+1+n.x)%n.x][k]+heightMap[(j-1+n.x)%n.x][k]+
										heightMap[j][(k-1+n.y)%n.y]+heightMap[j][(k+1+n.y)%n.y])/2.0+
										(heightMap[(j-1+n.x)%n.x][(k-1+n.y)%n.y]+heightMap[(j-1+n.x)%n.x][(k+1+n.y)%n.y]+
										heightMap[(j+1+n.x)%n.x][(k-1+n.y)%n.y]+heightMap[(j+1+n.x)%n.x][(k+1+n.y)%n.y])/4.0-
										3.0*heightMap[j][k]),2.0)/(area);
					//gradSquared+=pow(da.x+da.y+db.x+db.y,2.0)/(r.x*r.y);
					//gradSquaredMap[j][k]=pow(da.x+da.y+db.x+db.y,2.0)/(r.x*r.y);
					
					//gradSquared-=(da.x*db.x/(r.x*r.x)+da.y*db.y/(r.y*r.y));
					//gradSquared+=(h.x+h.y)/2.0;//mean curvature
				}
			}
			//std::cout << gradSquared << '\t';// << hSum.x << '\t' << hSum.y << std::endl;
			//std::cerr << "Surface " << s << " gradSquared=" << gradSquared << " vectSums=" << aSum.x << ',' << aSum.y << ',' << bSum.x << ',' << bSum.y;
			std::cout << gradSquared << std::endl;
			std::string buf;
			
			buf="gradSquaredMap_N";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << s;
			parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".xyz";
			std::fstream dataFile;
			//dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			//dataFile << n.x*n.y << "\ntest\n";
			//dump the heightmap to output
			//for(int i=0;i<gradSquaredMap.size();i++)
			//{
			//	for(int j=0;j<gradSquaredMap[i].size();j++)
			//	{
			//		dataFile << "1\t" << static_cast<double>(i)*r.x+r.x/2.0 << '\t' << static_cast<double>(j)*r.y+r.y/2.0 << '\t' << gradSquaredMap[i][j] << std::endl;
			//	}
			//}
			//dataFile << std::endl;
			//dataFile.close();
			
			buf.clear();
			buf="gradSquared_N";
			////std::string whatever;
			//std::stringstream parseNumber;
			//parseNumber << s;
			//parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".dat";
			//std::fstream dataFile;
			//dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			//dataFile << n.x*n.y << "\ntest\n";
			//dump the heightmap to output
			//for(int i=0;i<heightMap.size();i++)
			//{
			//	int i=40;
			//	for(int j=0;j<heightMap[i].size();j++)
			//	{
			//		
			//		//dataFile << "1\t" << static_cast<double>(i)*r.x+r.x/2.0 << '\t' << static_cast<double>(j)*r.y+r.y/2.0 << '\t' << heightMap[i][j] << std::endl;
			//		dataFile << static_cast<double>(j)*r.y+r.y/2.0 << '\t' << heightMap[i][j] << '\t' << gradSquaredMap[i][j] << std::endl;
			//	}
			//	dataFile << std::endl;
			//}
			//dataFile << std::endl;
			//dataFile.close();
			
			buf.clear();
			buf="size_N";
			//std::string whatever;
			//std::stringstream parseNumber;
			//parseNumber << s;
			//parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".xyz";
			//std::fstream dataFile;
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile << time << '\t' << size.x << '\t' << size.y << std::endl;
			dataFile.close();
			
			buf.clear();
			buf="gradSquared_N";
			//std::string whatever;
			//std::stringstream parseNumber;
			//parseNumber << s;
			//parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".dat";
			//std::fstream dataFile;
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			dataFile << time << '\t' << gradSquared << std::endl;
			dataFile.close();
			
			std::cerr << std::endl;
		}
		if(surfaces.size()==2)
			std::cout << std::endl;
		
	}	
	return 0;
}
