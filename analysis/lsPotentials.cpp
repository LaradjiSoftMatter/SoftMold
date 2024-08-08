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
	if(argc<2 || argc>4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name cutoff rmin\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//std::string newName("frames_");
	//newName+=name;
	//newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	int nTypes=System.readNTypes();
	double *UC=System.getTwoBodyUconst();
	std::vector<std::vector<double>> Umax, Umin;
	//#define nTWOBODYUCONST 6
	for(int i=0;i<nTypes;i++)
	{
		std::vector<double> n,x;
		for(int j=0;j<nTypes;j++)
		{
			int index=(i+j*nTypes)*nTWOBODYUCONST;
			n.push_back(UC[index+2]);
			x.push_back(UC[index+1]*(UC[index]*UC[index])+UC[index+2]);
		}
		Umin.push_back(n);
		Umax.push_back(x);
	}
	
	std::cout << "Umin" << ' ';
	std::cout << nTypes << ' ' << nTypes << std::endl;
	for(auto &i:Umin)
	{
		for(auto &Um:i)
			std::cout << Um << '\t';
		std::cout << std::endl;
	}
	
	std::cout << "Umax" << ' ';
	std::cout << nTypes << ' ' << nTypes << std::endl;
	for(auto &i:Umax)
	{
		for(auto &Um:i)
			std::cout << Um << '\t';
		std::cout << std::endl;
	}
	
	//System.addTwoBodyUconst(RMIN);//C8
	//System.addTwoBodyUconst((Umax[g]-Umin[g])/(RMIN*RMIN));//C4
	//System.addTwoBodyUconst(Umin[g]);//C5,no index trick
	//System.addTwoBodyUconst(CUTOFF);//C7
	//System.addTwoBodyUconst((3.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C1
	//System.addTwoBodyUconst((2.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C0
	
	//xyzFormat<double> xyzFile(p,nParticles);
	//xyzFile.open(newName.c_str(), std::ios::in);
	
	//while(xyzFile.load())
	//{
	//	
	//}
	
	return 0;
}

