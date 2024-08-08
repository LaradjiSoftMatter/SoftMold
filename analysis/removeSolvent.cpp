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
		std::cout << "usage: " << argv[0] << " name outname.xyz\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	char *outname=argv[2];
	
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
	
	//allocate some new positions for the non solvent variety
	int nOut=0;
	for(int i=0;i<nParticles;i++)
		#ifdef SOLVENT_FLAG
			if(p[i].type!=SOLVENT_FLAG)
		#endif
		#ifdef SOLVENT
			if(p[i].type!=SOLVENT)
		#endif
			nOut++;
	position<double> *pOut=new position<double>[nOut];
	
	xyzFormat<double> output(pOut, nOut);
	output.open(outname, std::ios::out);
	
	for(int frame=0;xyzFile.load();frame++)
	{
		std::cout << "Storing frame " << frame << ".\n";
		
		//Do data measurements in here
		int a=0;
		for(int i=0, j=0;i<nParticles;i++)
			#ifdef SOLVENT_FLAG
				if(p[i].type!=SOLVENT_FLAG)
			#endif
			#ifdef SOLVENT
				if(p[i].type!=SOLVENT)
			#endif
				pOut[j++]=p[i];
			
		output.store();
	}
	
	delete pOut;
	
	return 0;
}

