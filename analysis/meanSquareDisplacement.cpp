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

//For potential energy
#include "cellInclude/cell.h"

int main(int argc, char **argv)
{
	if(argc!=4 && argc!=5)
	{
		std::cerr << "Usage: " << argv[0] << " name fromTime type\n";
		std::cerr << "Usage: " << argv[0] << " name fromTime type cutoff\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime;
	int type;
	cmdArg >> fromTime >> type;
	
	double cutoff=0;
	if(argc==5)
		cmdArg >> cutoff;
	
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
		
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<int> ind;
	for(int i=0;i<nParticles;i++)
		if(p[i].type==type)
			ind.push_back(i);
	std::vector<threeVector<double>> initial,last;
	std::vector<std::vector<threeVector<double>>> paths;//paths per each particle
	//color flag per each particle, for example:
	//0,1,2,3,4,4,4,4,4,5,6,7,8,9,...,n
	//where at times 0 to 4 and 9 to n, the particle is in contact, and
	// at times 5 to 8, the particle isn't in contact.
	//Subtracting a time, t, from an initial time, tf, will tell you if it is a
	// continuous change. 
	std::vector<std::vector<int>> colors;
	std::vector<double> times;
		
	int frames=0;
	while(xyzFile.load())
	{
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
		if(time>=fromTime)
		{
			position<double> s;
			s.x=size.x;
			s.y=size.y;
			s.z=size.z;
			double cutoffSqr=cutoff*cutoff;
			threeVector<int> nCells=0;
			threeVector<double> cSize=0;
			umapfwd<int> iMap;
			if(cutoff!=0)
			{
				nCells.x=size.x/cutoff;
				nCells.y=size.y/cutoff;
				nCells.z=size.z/cutoff;
				cSize.x=size.x/nCells.x;
				cSize.y=size.y/nCells.y;
				cSize.z=size.z/nCells.z;
				iMap=createIMap(p,nParticles,nCells,cSize);
			}
			
			std::vector<int> color(ind.size());
			std::vector<threeVector<double>> current(ind.size()),path;
			for(int i=0;i<ind.size();i++)
			{
				if(cutoff==0 || initial.size()==0)
				{
					current[i].x=p[ind[i]].x;
					current[i].y=p[ind[i]].y;
					current[i].z=p[ind[i]].z;
					if(initial.size()==0)
						color[i]=0;
					else if(cutoff==0)
						color[i]=colors.back()[i]+1;
						
				}
				else
				{
					bool inRange=false;
					auto iCell=getCell(p[ind[i]],cSize);
					auto neigh=neighIndices(hash(iCell,nCells),nCells);
					for(auto nHash:neigh)
					{
						auto iMapIt=iMap.find(nHash);
						if(iMapIt!=iMap.end())
						for(auto &j:iMapIt->second)
						{
							if(ind[i]!=j)
							{
								threeVector<double> d;
								d.x=p[ind[i]].x-p[j].x;
								d.x-=(d.x>size.x/2.0)?size.x:0;
								d.x+=(d.x<-size.x/2.0)?size.x:0;
								d.y=p[ind[i]].y-p[j].y;
								d.y-=(d.y>size.y/2.0)?size.y:0;
								d.y+=(d.y<-size.y/2.0)?size.y:0;
								d.z=p[ind[i]].z-p[j].z;
								d.z-=(d.z>size.z/2.0)?size.z:0;
								d.z+=(d.z<-size.z/2.0)?size.z:0;
								if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
									inRange=true;
							}
						}
					}
					current[i].x=p[ind[i]].x;
					current[i].y=p[ind[i]].y;
					current[i].z=p[ind[i]].z;
					if(inRange)
						color[i]=colors.back()[i]+1;
					else//pretend it didn't move
						color[i]=colors.back()[i];
				}
			}
			if(initial.size()==0)
			{
				initial=current;
				path=current;
			}
			else
			{
				path=paths.back();
				for(int i=0;i<ind.size();i++)
				{
					//use last configuration to remove boundary
					threeVector<double> d;
					d.x=current[i].x-last[i].x;
					d.x-=(d.x>size.x/2.0)?size.x:0;
					d.x+=(d.x<-size.x/2.0)?size.x:0;
					d.y=current[i].y-last[i].y;
					d.y-=(d.y>size.y/2.0)?size.y:0;
					d.y+=(d.y<-size.y/2.0)?size.y:0;
					d.z=current[i].z-last[i].z;
					d.z-=(d.z>size.z/2.0)?size.z:0;
					d.z+=(d.z<-size.z/2.0)?size.z:0;
					path[i].x+=d.x;
					path[i].y+=d.y;
					path[i].z+=d.z;
				}
			}
			last=current;
			paths.push_back(path);
			colors.push_back(color);
				
			frames++;
			std::cerr << time << std::endl;
			times.push_back(time);
		}
		time+=System.readStoreInterval();
	}
	
	std::vector<double> msDs(paths.size(),0), counts(paths.size(),0);
	for(int tf=0;tf<paths.size()-1;tf++)
	for(int t=tf+1;t<paths.size();t++)
	{
		std::vector<threeVector<double>> &first=paths[tf];
		double msD=0;
		double count=0;
		for(int i=0;i<paths[t].size();i++)
		{
			if(colors[t][i]-colors[tf][i]==t-tf)
			{
				threeVector<double> d;
				d.x=paths[t][i].x-first[i].x;
				d.y=paths[t][i].y-first[i].y;
				d.z=paths[t][i].z-first[i].z;
				msD+=d.x*d.x+d.y*d.y+d.z*d.z;
				count++;
			}
		}
		if(count!=0)
		{
			msDs[t-tf]+=msD/count;
			counts[t-tf]++;
		}
	}
	for(int t=0;t<msDs.size();t++)
		if(counts[t]!=0)
			msDs[t]/=counts[t];
	
	for(int t=0;t<msDs.size();t++)
		std::cout << times[t] << ' ' << msDs[t] << std::endl;
	
	return 0;
}

