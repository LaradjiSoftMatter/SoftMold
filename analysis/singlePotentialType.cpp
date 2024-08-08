//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name type\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	int type=atoi(argv[2]);
	
	///Configuration variables
	Blob<double> System;
	
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	CellOpt<double, Potential<double>, Force <double> > pairInteractions(System.getPositions(), System.getAccelerations(), 
		System.getTwoBodyFconst(), System.getTwoBodyUconst(), System.readNParticles(), System.readNTypes(), System.readSize(),
		System.readPeriodic(), System.readCutoff());
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	int nType=0;
	for(int i=0;i<System.readNParticles();i++)
		if(p[i].type==type)
			nType++;
	
	int *index=new int[nType];
	
	nType=0;
	
	for(int i=0;i<System.readNParticles();i++)
		if(p[i].type==type)
		{
			index[nType]=i;
			nType++;
		}
	
	for(int i=0;xyzFile.load();i++)
	{
		pairInteractions.build();
		std::cout << System.readStoreInterval()*(double)i;
		int a=0;
		for(int j=0;j<nType;j++)
		{
			double potential=pairInteractions.computePotential(index[j]);
			if(potential>-10.0)
				a++;
			//std::cout << '\t' << potential;
		}
		std::cout << '\t' << a <<'\n';
	}
	
	delete index;
	
	return 0;
}

