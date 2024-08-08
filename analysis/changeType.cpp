#include <iostream>
#include <cstdlib>
#include <random>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char **argv)
{
	if(argc!=4 && argc!=5)
	{
		std::cout << "Usage: " << argv[0] << " name fromType toType" << std::endl;
		std::cout << "Usage: " << argv[0] << " name fromType toType fraction" << std::endl;
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fraction=1.0;
	int fromType,toType;
	cmdArg >> fromType >> toType;
	if(argc==5)
		cmdArg >> fraction;
	
	if(fraction>1.0 || fraction<0.0)
	{
		std::cerr << "Fraction must be in range (0,1)!" << std::endl;
		return -1;
	}
	
	//the variables for the simulation
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	if(System.readNTypes()<=fromType)
	{
		std::cerr << "fromType out of range!" << std::endl;
		std::cerr << "\tnTypes: " << System.readNTypes() << '\t' << "fromType: " << fromType << std::endl;
		return -1;
	}
	
	if(System.readNTypes()<=toType)
	{
		std::cerr << "toType out of range!" << std::endl;
		std::cerr << "\tnTypes: " << System.readNTypes() << '\t' << "toType: " << toType << std::endl;
		return -1;
	}
	
	auto p=System.getPositions();
	std::vector<int> pInd;
	for(int i=0;i<System.readNParticles();i++)
		if(p[i].type==fromType)
			pInd.push_back(i);
	std::shuffle(pInd.begin(),pInd.end(),std::mt19937(System.readSeed()));
	int nChanges=0;
	for(const auto &i:pInd)
		if(nChanges<fraction*pInd.size())
		{
			p[i].type=toType;
			nChanges++;
		}
	std::cerr << "Changed " << nChanges << " particle types." << std::endl;
	
	std::cerr << "Storing configuration...";
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


