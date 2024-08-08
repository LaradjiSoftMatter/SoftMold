#define ANCHOR_DATA

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name\n";
		std::cout << "Lays out trajectories as system states.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	//std::stringstream cmdArg;
	//cmdArg << argv[2];
	//int skip;
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
	
	for(int i=0;i<frames[0].size();i++)
	{
		std::cout << frames.size() << "\ntest\n";
		for(int j=0;j<frames.size();j++)
		{
			std::cout << (j!=0) << '\t' << frames[j][i].x << '\t' << frames[j][i].y << '\t' << frames[j][i].z << std::endl;
		}
	}
	
	return 0;
}
