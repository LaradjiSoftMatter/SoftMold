//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

//For molecular dynamics forces and potentials
#include "include/MD.h"

//For the molecular dynamics variables
#include "include/system.h"

struct nonBonded {
	double Rmin,Umax,Umin,cutoff;
};

struct bond {
	double l,k;
};

struct bend {
	double cosTheta,k;
};

std::vector<std::vector<nonBonded> > twoBodyUConstMatrix(double *input, int nTypes)
{
	int nConst=6;
	std::vector<std::vector<nonBonded> Umatrix;
	for(int i=0;i<nTypes;i++)
	{
		std::vector<nonBonded> Urow;
		for(int j=0;j<nTypes;j++)
		{
			nonBonded U;
			int RminIndex=nConst*(i*nTypes+j);
			int UmaxIndex=RminIndex+1;
			int UminIndex=RminIndex+2;
			int cutoffIndex=RminIndex+3;
			U.Rmin=input[RminIndex];
			U.Umax=input[UmaxIndex];
			U.Umin=input[UminIndex];
			U.cutoff=input[cutoffIndex];
			U.Umax*=U.Rmin*U.Rmin;
			U.Umax+=U.Umin;
			Urow.push_back(U);
		}
		Umatrix.push_back(Urow);
	}
	return Umatrix;
}

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		std::cout << argv[0] << " name\n";
		return 0;
	}
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	int nTypes=System.readNTypes();
	double *Uconst=System.getTwoBodyUconst();
	std::vector<std::vector<nonBonded> > Umatrix=twoBodyUConstMatrix(Uconst,nTypes);
	std::cout << "Umax";
	for(int i=0;i<nTypes;i++)
		std::cout << "," << i;
	std::cout << std::endl;
	for(int i=0;i<nTypes;i++)
	{
		std::cout << i;
		for(int j=0;j<nTypes;j++)
			std::cout << "," << Umatrix[i][j].Umax;
		std::cout << std::endl;
	}
	
	
	int nMol=System.readNMolecules();
	molecule<double,fourVector<int> > *m;
	
	//close the file that it is reading from
	fileIO.close();
	
	
	return 0;
}
