//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#define STAT_OUT

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name cutoff\n";
		std::cout << "Output is volume of each inner space consisting of only solvent.\n";
		std::cout << "Extra output is a volume.xyz file. The file shows the volume number\n";
		std::cout << "as a type of particle on a grid for the initial space.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	double cutoff=atof(argv[2]);
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	///Map and count the volumes
	int *excludeType=new int[System.readNTypes()];
	for(int i=1;i<System.readNTypes();i++)
		excludeType[i-1]=i;
	
	VolumeExtraction<double> volumize(System.getPositions(), System.readNParticles(),\
		System.readSize(), cutoff, excludeType, System.readNTypes()-1, System.readSeed());
	
	
	double start=omp_get_wtime();
	volumize.build();
	double end=omp_get_wtime();
	std::cout << "To run build: " << end-start << '\n';
	
	
	for(int i=0;i<1000;i++)
	{
		int inner=volumize.grabInner();
		volumize.moveToOuter(inner);
	}
	
	std::cout << "nExchanged: " << volumize.nExchanged() << '\n';
	
	volumize.build();
	std::cout << "nExchanged: " << volumize.nExchanged() << '\n';
	
	volumize.moveToOuter(volumize.grabInner());
	
	volumize.build();
	std::cout << "nExchanged: " << volumize.nExchanged() << '\n';
	
	
	/*
	std::fstream outFile;
	outFile.open("test.xyz",std::ios::out);
	
	position<double> *p=System.getPositions();
	outFile << System.readNParticles() << '\n' << "test\n";
	for(int i=0;i<System.readNParticles();i++)
		outFile << volumize.readVolumeIndex(i) << '\t' << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\n';
	
	outFile.close();
	*/
	volumize.exportMap("vMap.xyz",std::ios::out);
	
	delete excludeType;
	
	return 0;
}