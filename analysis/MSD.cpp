#define ANCHOR_DATA

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name molIndex\n";
		std::cout << "Extracts data in place of MD.cpp for already executed systems.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	std::stringstream cmdArg;
	cmdArg << argv[2];
	//int skip;
	int molIndex;
	cmdArg >> molIndex;
	//cmdArg >> skip;
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	//frames[frame][particle].coordinate
	std::vector< std::vector< position<double> > > frames;
	
	//load and store every frame
	for(int i=0;xyzFile.load();i++)
	{
		double time=static_cast<double>(i)*System.readStoreInterval();
		std::cerr << static_cast<double>(i)*System.readStoreInterval() << std::endl;
		
		std::vector<position<double> > frame;
		
		for(int j=0;j<nParticles;j++)
			frame.push_back(p[j]);
		frames.push_back(frame);
	}
	
	xyzFile.close();
	
	//MSD[particle][frame]
	std::vector< std::vector<double> > MSD;
	
	//calculate MSD per particle, while averaging every skipped frame
//	for(int i=1;i<frames[0].size();i++)//particle
//	{
		for(int skip=1;skip<frames.size()/2;skip++)
		//int skip=frames.size()/2-1;
		{
			for(int i=1;i<frames[0].size();i++)//particle
			{
			double displacement=0;
			int nValues=0;
			for(int j=skip;j<frames.size();j+=skip)//skip displacement
			{
				threeVector<double> d;
				d.x=frames[j][i].x-frames[j-skip][i].x;
				d.y=frames[j][i].y-frames[j-skip][i].y;
				d.z=frames[j][i].z-frames[j-skip][i].z;
				
				displacement+=(d.x*d.x+d.y*d.y+d.z*d.z);
				nValues++;
			}
			if(nValues!=0)
				displacement/=static_cast<double>(nValues);
//			std::cout << static_cast<double>(skip)*System.readStoreInterval() << '\t' << displacement << std::endl;
			std::cout << i << '\t' << displacement << std::endl;
		}
		std::cout << std::endl;
		//if(i%10000==0)
		//	std::cerr << i << std::endl;
		std::cerr << skip << std::endl;
		
	}
	
	return 0;
}
