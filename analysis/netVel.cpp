//determine if there is a net velocity in the system

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <ctime>

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name\n";
		return 0;
	}
	
	char *name=argv[1];//Aaaaaaaargh, a V matey!
	
	//the variables for the simulation
	Blob<double> System;
	
	//load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	threeVector<double> netVel;
	
	netVel=pairwiseSum< threeVector<double> >(System.getVelocities(),0,System.readNParticles());
	netVel/=(double)System.readNParticles();
	std::cout << "Net velocity (x,y,z): (" << netVel.x << ", " << netVel.y << ", " << netVel.z << ")\n";
	
	//fileIO.open(name,std::ios::out);
	//fileIO.write();
	//fileIO.close();
		
	return 0;
}
