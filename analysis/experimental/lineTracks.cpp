#define ANCHOR_DATA

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name nTrials frame depth\n";
		std::cout << "Calculates possible hops along the membrane.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	std::stringstream cmdArg;
	cmdArg << argv[2] << '\t' << argv[3] << '\t' << argv[4];
	int nTrials, frame, depth;
	cmdArg >> nTrials >> frame >> depth;
	
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
	double cutoff=System.readCutoff();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	xyzFile.load(frame);
	
	xyzFile.close();
	
	//trial paths for every particle
	for(int i=0;i<nParticles;i++)
	{
		double avgPathLength=0;
		for(int trial=0;trial<nTrials;trial++)
		{
			double pathLength=0;
			for(int hop=0;hop<depth;hop++)
			{
				
			}
			avgPathLength+=pathLength;
		}
		avgPathLength/=nTrials;
		std::cout << i << '\t' << avgPathLength << std::endl;
	}
	
	return 0;
}
