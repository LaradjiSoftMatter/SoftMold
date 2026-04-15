//Adjusts bonding constants. Just uses Blob and script to read, modify, and write values.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name\n";
		return 0;
	}
	
	char *name=argv[1];
	
	//the variables for the simulation
	Blob<double> System;
	
	//load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	molecule<double, fourVector<int> > *m=System.getMolecule();
	//New values:
	for(int i=0;i<System.readNMolecules();i++)
	{
		std::cout << "mIndex " << i << " is a ";
		switch(m[i].readType())
		{
			case CHAIN:
				std::cout << "CHAIN type  with " <<
				m[i].getBonds()[0].s[NCHAINS] << " number of chains, " <<
				m[i].getBonds()[0].s[CHAINLENGTH] << " number of monomers, " <<
				" start index=" << m[i].getBonds()[0].s[START] << 
				", bondlength=" << m[i].getConstants()[0] <<
				", k_bond=" << m[i].getConstants()[1] <<
				", -cos(theta_0)=" << m[i].getConstants()[2] <<
				", and k_bend=" << m[i].getConstants()[3] << std::endl;
				break;
			case BEND:
				std::cout << "BEND type  with " <<
				m[i].readNBond() << " bends" <<
				", -cos(theta_0)=" << m[i].getConstants()[0] <<
				", and k_bend=" << m[i].getConstants()[1] << std::endl;
				break;
			case BOND:
				std::cout << "BOND type  with " <<
				m[i].readNBond() << " bonds" <<
				", bondlength=" << m[i].getConstants()[0] <<
				", k_bond=" << m[i].getConstants()[1] << std::endl;
				break;
			case BEAD:
				std::cout << "BEAD type  with " <<
				m[i].readNBond() << " particles " << std::endl;
				break;
			case BALL:
				std::cout << "BALL type  with " <<
				m[i].readNBond() << " particles " << std::endl;
				break;
			case NANOCORE:
				std::cout << "NANOCORE type  with " <<
				m[i].readNBond() << " particles " << std::endl;
				break;
			case BOUNDARY:
				std::cout << "BOUNDARY type  with " <<
				m[i].readNBond() << " particles " << std::endl;
				break;
			default:
				std::cout << "m[i].readType()=" << m[i].readType() << std::endl;
				break;
		}
	}
	
	return 0;
}
