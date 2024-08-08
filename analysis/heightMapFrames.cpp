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
		
		//set up a quick search using the minimum surface spacing
		twoVector<int> surfN;
		surfN.x=static_cast<int>(size.x/surfSpacing);
		surfN.y=static_cast<int>(size.y/surfSpacing);
		twoVector<double> surfR;
		surfR.x=size.x/static_cast<double>(surfN.x);
		surfR.y=size.y/static_cast<double>(surfN.y);
		
		typedef std::vector<int> cellType;//fastest is vector, deque is a little slower
		
		std::map<int,cellType > cells;
		
		for(int j=0;j<pIndex.size();j++)
		{
			int i=pIndex[j];
			twoVector<int> c;
			c.x=static_cast<int>(p[i].x/surfR.x)%surfN.x;
			c.y=static_cast<int>(p[i].y/surfR.y)%surfN.y;
			int hash=c.x+c.y*surfN.x;
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
					twoVector<int> c;
					c.x=static_cast<int>(p[j].x/surfR.x)%surfN.x;
					c.y=static_cast<int>(p[j].y/surfR.y)%surfN.y;
					//int hash=c.x+c.y*surfN.x;
					
					//find connected elements and push onto stack
					for(int a=-1;a<2;a++)
					{
						for(int b=-1;b<2;b++)
						{
							twoVector<int> neighbor;
							neighbor.x=(c.x+a+surfN.x)%surfN.x;
							neighbor.y=(c.y+b+surfN.y)%surfN.y;
							int neighborHash=neighbor.x+neighbor.y*surfN.x;
							
							//std::cerr << hash << '\t' << neighborHash << std::endl;
							
							std::map<int,cellType >::iterator it=cells.find(neighborHash);
							for(int k=0;it!=cells.end() && k<it->second.size();k++)
							//for(int k=0;k<pIndex.size();k++)
							{
								threeVector<double> d;
								
								int l=pIndex[it->second[k]];
								//int l=pIndex[k];
								
								d.x=p[l].x-p[j].x;
								d.y=p[l].y-p[j].y;
								d.z=p[l].z-p[j].z;
								d.x-=(d.x>size.x/2.0)?size.x:0.0;
								d.y-=(d.y>size.y/2.0)?size.y:0.0;
								d.z-=(d.z>size.z/2.0)?size.z:0.0;
								d.x+=(d.x<-size.x/2.0)?size.x:0.0;
								d.y+=(d.y<-size.y/2.0)?size.y:0.0;
								d.z+=(d.z<-size.z/2.0)?size.z:0.0;
								
								if(l!=j && sFlag[it->second[k]]==-1 && d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
								//if(l!=j && sFlag[k]==-1 && d.x*d.x+d.y*d.y+d.z*d.z<surfSpacingSqr)
								{
									sFlag[it->second[k]]=surfaces.size();
									surface.push_back(l);
									
									stack.push_back(it->second[k]);
									
									//std::cout << d.x*d.x+d.y*d.y+d.z*d.z << '\t' << hash << '\t' << neighborHash << std::endl;
									
									//sFlag[k]=surfaces.size();
									//stack.push_back(k);
								}
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
		
		//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		//std::cerr << time_span.count() << std::endl;
		
		std::cerr << surfaces.size() << '\t';
		
		//make the height map
		threeVector<double> r;
		r.x=size.x/static_cast<double>(n.x);
		r.y=size.y/static_cast<double>(n.y);
		
		
		for(int s=0;s<surfaces.size() && dry==0;s++)
		{
			std::vector<double> heightMap(n.x*n.y, 0.0);
			std::vector<int> countMap(n.x*n.y, 0.0);
			
			for(int j=0;j<surfaces[s].size();j++)
			{
				int i=surfaces[s][j];
				threeVector<int> c;
				c.x=static_cast<int>(p[i].x/r.x);
				c.y=static_cast<int>(p[i].y/r.y);
				c.x=(c.x>=n.x)?n.x-1:c.x;
				c.y=(c.y>=n.y)?n.y-1:c.y;
				int hash=c.x+c.y*n.x;
				if(hash>=heightMap.size())
				{
					std::cerr << "Out of bounds particle!\n";
					std::cerr << "Size: " << size.x << '\t' << size.y << std::endl;
					std::cerr << "Position: " << p[i].x << '\t' << p[i].y << std::endl;
				}
				heightMap[hash]+=p[i].z;
				countMap[hash]++;
			}
			
			std::string buf;
			buf="surface_N";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << s;
			parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			buf+=argv[1];
			buf+=".dat";
			std::fstream dataFile;
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			
			//dump the heightmap to output
			for(int i=0;i<heightMap.size();i++)
			{
				if(countMap[i]!=0)
					heightMap[i]/=static_cast<double>(countMap[i]);
				//dataFile << static_cast<double>(i%n.x)*r.x+r.x/2.0 << '\t' << floor(i/n.x)*r.y+r.y/2.0 << '\t' << heightMap[i] << std::endl;
				dataFile << i << '\t' << heightMap[i] << std::endl;
				//if(i%n.y==n.y-1)
				//	std::cout << std::endl;
			}
			//dataFile << std::endl;
			dataFile.close();
		}
		
		std::cerr << std::endl;
	}	
	return 0;
}
