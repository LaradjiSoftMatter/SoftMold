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
	if(argc!=2 || argc!=3 || argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "To see available chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name\n";
		
		std::cout << "To get the average length of chain type molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex\n";
		
		std::cout << "To get the length of a particular chain:\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex chainIndex\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
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
	molecule< fourVector<int> > *m=System.getMolecule();
	
	if(argc==2)
	{
		
	}
	if(argc==3)
	{
		
	}
	if(argc==4)
	{
		
	}
	
	/*
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	while(xyzFile.load())
	{
		//Do data measurements in here
		int a=0;
		for(int i=0;i<nParticles;i++)
			if(p[i].x>System.readSize().x/2.0)
				a++;
		std::cout << a << std::endl;
	}
	*/
	return 0;
}

