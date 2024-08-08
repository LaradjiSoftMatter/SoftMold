//Adjusts bonding constants. Just uses Blob and script to read, modify, and write values.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <ctime>

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name bondLength\n";
		return 0;
	}
	
	char *name=argv[1];//Aaaaaaaargh, a V matey!
	std::stringstream cmdArg;
	cmdArg << argv[2];
	double bondLength;
	cmdArg >> bondLength;
	
	//the variables for the simulation
	Blob<double> System;
	
	//load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	//New values:
	for(int i=0;i<System.readNMolecules();i++)
	{
		molecule<double, fourVector<int> > *m=System.getMolecule();
		if(m[i].readType()==CHAIN)
		{
			if(m[i].readNBond()==1 && m[i].getBonds()[0].s[CHAINLENGTH]>4)
			{
				//interesting note: while the pointer is safe (you can't change it internally),
				//			the values the pointer points to are not.
				m[i].getConstants()[CHAINBOND+ABOND]=bondLength;
			}
		}
	}
	
	fileIO.open(name,std::ios::out);
	fileIO.write();
	fileIO.close();
		
	return 0;
}
