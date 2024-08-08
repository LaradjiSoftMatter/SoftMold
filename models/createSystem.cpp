//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#define NO_REQUIRED_COMMANDS

#include "../include/systemMD.h"

int main(int argc, char* argv[])
{
	if(argc!=6)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Usage: " << argv[0] << " name modelName pos.x pos.y pos.z\n";
		return 0;
	}
	
	
	//read in options
	char *name=argv[1];
	char *modelName=argv[2];
	threeVector<double> pos;
	std::stringstream cmdIn;
	cmdIn << argv[3] << ' ' << argv[4] << ' ' << argv[5] << '\n';
	cmdIn >> pos.x >> pos.y >> pos.z;
	
	//Configuration variables
	Blob<double> System;
	Blob<double> model;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	Script<double, Blob <double> > modelIO(modelName,std::ios::in,&model);
	modelIO.read();
	modelIO.close();
	
	//get the new number of particles
	int newNParticles=System.readNParticles()+model.readNParticles();
	System.allocParticle(newNParticles);
	
	//the offset for inserting particles, needed to remap molecules from model
	int newParticleOffset=System.readNParticles();
	
	//Locate center of mass, we are placing at pos
	position<double> comValue=com< position<double> >(model.getPositions(), model.readNParticles());
	
	//move it to position and insert it
	for(int i=0;i<model.readNParticles();i++)
	{
		position<double> newPos=model.getPositions()[i];
		newPos.x+=pos.x-comValue.x;
		newPos.y+=pos.y-comValue.y;
		newPos.z+=pos.z-comValue.z;
		
		System.addParticle(newPos,model.getVelocities()[i],model.getAccelerations()[i]);
	}
	
	//remap molecules
	for(int i=0;i<model.readNMolecules();i++)
	{
		//shorter names
		molecule< double, fourVector<int> > *m=model.getMolecule();
		int nBonded=m[i].readNBond();
		fourVector<int> *bond=m[i].getBonds();
		
		//our new molecule
		molecule<double,fourVector<int> > newMolecule;
		newMolecule=m[i];
		
		switch(m[i].readType())
		{
			case CHAIN:
				std::cout << "Copying molecule " << i << ", a chain type.\n";
				
				//remap all particles within molecule
				for(int k=0;k<nBonded;k++)
				{
					fourVector<int> newBond;
					newBond.s[START]=newParticleOffset+bond[k].s[START];
					newBond.s[NCHAINS]=bond[k].s[NCHAINS];
					newBond.s[LENGTH]=bond[k].s[LENGTH];
					
					//Actual remapped bond
					//This doesn't work
					newMolecule.setBond(k,newBond);
				}
				
				
				break;
			case BOND:
				std::cout << "Copying molecule " << i << ", a bond type.\n";
				
				//map all particles within molecule
				for(int k=0;k<nBonded;k++)
				{
					fourVector<int> newBond;
					newBond.s[0]=newParticleOffset+bond[k].s[0];
					newBond.s[1]=newParticleOffset+bond[k].s[1];
					
					//Actual remapped bond
					newMolecule.setBond(k,newBond);
				}
				break;
			case BEND:
				std::cout << "Copying molecule " << i << ", a bend type.\n";
				//map all particles within molecule
				for(int k=0;k<nBonded;k++)
				{
					fourVector<int> newBond;
					newBond.s[0]=newParticleOffset+bond[k].s[0];
					newBond.s[1]=newParticleOffset+bond[k].s[1];
					newBond.s[2]=newParticleOffset+bond[k].s[2];
					
					//Actual remapped bond
					newMolecule.setBond(k,newBond);
				}
				break;
			default:
				std::cout << "Molecule " << i << " isn't recognized.\n";
				break;
		}
		System.addMolecule(newMolecule);
	}
	
	
	//Save system
	fileIO.open(name,std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

