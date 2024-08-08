/**
 * \brief Simple implicit molecular dynamic system. Just a bilayer in various geometries.
 * The steps to make a system are pretty simple, especially given the lipidModels.h header.
 * Just put together the parameters for initial conditions, generate the geometry, and
 * then run the system. That's it! If you want more, less, or something different, just
 * modify this program. There is a section that takes the command line arguments (just
 * after the main function definition) and checks that there are enough. After that,
 * atoi() and atof() (see some cstdlib documentation) convert the text to values. Those
 * are basically your initial conditions. The Blob class template is also useful for this.
 * Just take a parameter you want to modify and either use setParameter, addParameter, or
 * delParameter members in Blob to modify it. I've used the variable definition
 * 'Blob\<double\> System' for convienience. Functions in lipidModels.h accept System as
 * a reference, and can modify parameters as one would in main. For convienience, this
 * program outputs a name.xyz file to quickly view the initial geometry; this is useful
 * for debugging.
 */

#include <iostream>
#include <cstdlib>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include "../models/lipidModels.h"
#include "../models/nanoModels.h"
using namespace std;

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
		std::cout << "Usage: " << argv[0] << " name size.x size.y size.z" << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	threeVector<double> size;
	cmdArg >> size.x >> size.y >> size.z;
	
	//the variables for the simulation
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	threeVector<double> oldSize=System.readSize();
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	for(int i=0;i<nParticles;i++)
	{
		p[i].x*=size.x/oldSize.x;
		p[i].y*=size.y/oldSize.y;
		p[i].z*=size.z/oldSize.z;
	}
	System.setSize(size);
	
	std::cerr << "Storing configuration...";
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


