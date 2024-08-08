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
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name offsetX offsetY offsetZ\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	threeVector<double> offset;
	cmdArg >> offset.x >> offset.y >> offset.z;
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	position<double> *p=System.getPositions();
	threeVector<double> *v=System.getVelocities();
	threeVector<double> *a=System.getAccelerations();
	
	int nParticles=System.readNParticles();
	
	for(int i=0;i<nParticles;i++)
	{
		position<double> buf=p[i];
		buf.x+=offset.x;
		buf.y+=offset.y;
		buf.z+=offset.z;
		//this automatically puts particles within system bounds
		System.setParticle(i, buf, v[i], a[i]);
	}
	
	//Save it
	fileIO.open(name,std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

