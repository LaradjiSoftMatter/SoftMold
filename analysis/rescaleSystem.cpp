#include <iostream>
#include <cstdlib>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		std::cout << "Usage: " << argv[0] << " name xratio yratio zratio" << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double xratio,yratio,zratio;
	cmdArg >> xratio >> yratio >> zratio;
	
	//the variables for the simulation
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	if(xratio<=0 || yratio<=0 || zratio<=0)
	{
		std::cerr << "Ratios cannot be zero or negative!" << std::endl;
		std::cerr << "\tnTypes: (" << xratio << ", " << yratio << ", " << zratio << ")" << std::endl;
		return -1;
	}
	
	auto p=System.getPositions();
	for(int i=0;i<System.readNParticles();i++)
	{
		p[i].x*=xratio;
		p[i].y*=yratio;
		p[i].z*=zratio;
	}
	auto size=System.readSize();
	size.x*=xratio;
	size.y*=yratio;
	size.z*=zratio;
	System.setSize(size);
	
	std::cerr << "Storing configuration...";
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


