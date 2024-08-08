//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#include "../include/systemMD.h"
#include "lipidModels.h"

int main(int argc, char* argv[])
{
	if(argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name density removeSolventFraction\n";
		return 0;
	}
	
	char *name=argv[1];
	double density=atof(argv[2]);
	double removeSolventFraction=atof(argv[3]);
	
	///Previous Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	//Steal underwear
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	///Check to make sure configuration is correct
	//?????????
	
	///This initializes the new variables
	//Profit!
	
	//Adjust system to minimize space for solvent
	threeVector<double> lower,upper;
	lower=System.readSize();
	upper=0;
		
	for(int i=0;i<System.readNParticles();i++)
	{
		position <double> p=System.getPositions()[i];
		lower.x=(lower.x>p.x)?p.x:lower.x;
		lower.y=(lower.y>p.y)?p.y:lower.y;
		lower.z=(lower.z>p.z)?p.z:lower.z;
		
		upper.x=(upper.x<p.x)?p.x:upper.x;
		upper.y=(upper.y<p.y)?p.y:upper.y;
		upper.z=(upper.z<p.z)?p.z:upper.z;
	}
	
	for(int i=0;i<System.readNParticles();i++)
	{
		System.getPositions()[i].x+=((upper.x-lower.x)*1.0/8.0-lower.x);
		System.getPositions()[i].y+=((upper.y-lower.y)*1.0/8.0-lower.y);
		System.getPositions()[i].z+=((upper.z-lower.z)*1.0/8.0-lower.z);
		if(System.getPositions()[i].x<0 || System.getPositions()[i].x>(upper.x-lower.x)*(10.0/8.0) ||
			System.getPositions()[i].y<0 || System.getPositions()[i].y>(upper.y-lower.y)*(10.0/8.0) ||
			System.getPositions()[i].z<0 || System.getPositions()[i].z>(upper.z-lower.z)*(10.0/8.0))
		{
			std::cout << "Particle out of box!\n";
			std::cout << i << '\n';
			std::cout << lower.x << ' ' << lower.y << ' ' << lower.z << '\n';
			std::cout << upper.x << ' ' << upper.y << ' ' << upper.z << '\n';
			std::cout << System.getPositions()[i].x << ' ' << System.getPositions()[i].y << ' ' << System.getPositions()[i].z << '\n';
			return 0;
		}
	}
	
	upper.x=(upper.x-lower.x)*(10.0/8.0);
	upper.y=(upper.y-lower.y)*(10.0/8.0);
	upper.z=(upper.z-lower.z)*(10.0/8.0);
	
	System.setSize(upper);
	
	//Add solvent
	solventFill<double> (System,density);
	
	//remove solvent amount
	System.setRemoveSolvent(removeSolventFraction);
	
	///This stores the new configuration
	std::string newName("");
	newName+=name;
	newName+="_wS";
	
	Script<double, Blob <double> > fileIO_new(newName.c_str(),std::ios::out,&System);
	fileIO_new.write();
	fileIO_new.close();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p, nParticles);
	newName+=".xyz";
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	return 0;
}

