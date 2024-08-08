/**
 */

#include <iostream>
#include <cstdlib>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#define CUTOFF 2.0
#define RMIN 1.0

#define HEAD 2
#define TAIL 3
#define HEAD2 6
#define TAIL2 7

#define NANOTYPE 4
//#define SINGLEBEAD 5

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		std::cout << "nArgs: " << argc << std::endl;
		std::cout << "Usage: " << argv[0] << " name mIndex constant(aka k/2) length" << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int mIndex;
	double constant, length;
	cmdArg >> mIndex >> constant >> length;
	
	//the variables for the simulation
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	if(System.readNMolecules()<=mIndex)
	{
		std::cerr << "mIndex out of range!" << std::endl;
		std::cerr << "\tSystem: " << System.readNMolecules() << '\t' << "mIndex: " << mIndex << std::endl;
		return 0;
	}
	
	if(System.getMolecule()[mIndex].readType()!=BOND)
	{
		std::cerr << "Wrong type!" << std::endl;
		std::cerr << "Should be of BOND type (" << BOND << ")!" << std::endl;
		return 0;
	}
	
	molecule< double, fourVector<int> > *m=System.getMolecule();
	double *mConstants=m[mIndex].getConstants();
	
	mConstants[0]=length;
	mConstants[1]=constant;
	
	std::cerr << "Storing configuration...";
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


