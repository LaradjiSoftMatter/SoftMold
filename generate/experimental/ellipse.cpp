#include <iostream>
#include <cstdlib>
//Include files from the library:

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For certain lipid models (liposome, flat membrane, etc...)
#include "../models/lipidModels.h"
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

int main(int argc, char **argv)
{
	if(argc<5)
	{
		cout << "Usage: " << argv[0] << " name seed nLipids arialDensity" << endl;
		return 0;
	}

	char *name=argv[1];
	int nLipids=atoi(argv[3]);
	double arialDensity=atof(argv[4]);
	
	
	///the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(6);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(CUTOFF);
	
	System.setInitialTime(0);
	System.setFinalTime(50000);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(10);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
	///initialize constants (force and potential constants)
	//generate force constants
	double *Umin=new double[System.readNTypes()*System.readNTypes()];
	double *Umax=new double[System.readNTypes()*System.readNTypes()];
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			Umin[k]=0;
			Umax[k]=100;
		}
	}
	
	//umin and umax exceptions, tail types
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	//two body constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			int g=i+j*System.readNTypes();
			//U[r<=rmin]=((Umax-Umin)*(rmin-r)^2/rmin^2)+Umin
			//U[rmin<r<=rc]=(-2*Umin*(rc-r)^3/(rc-rmin)^3)+(3*Umin*(rc-r)^2/(rc-rm)^2)
			//F[r<=rmin]=(-2*(Umax-Umin)/rmin^2)*(rmin-r)/r
			//F[rmin<r<=rc]=(6*Umin/(rc-rmin)^3)*(rc-r)^2/r-(6*Umin/(rc-rmin)^2)*(rc-r)/r
			
			//constants[U[rmin<r<=rc]]: 0:-2*Umin/(rc-rmin)^3 1:3*Umin/(rc-rmin)^2
			//constants[F[rmin<r<=rc]]: 2:-6*Umin/(rc-rmin)^3  3:6*Umin/(rc-rmin)^2
			
			//constants[U[r<=rmin]]: 4:(Umax-Umin)/rmin^2  5:Umin
			//constants[F[r<=rmin]]: 6:2*(Umax-Umin)/rmin^2
			//constants[general]: 7:rc 8:rmin 9:rc^2 10:rmin^2
			
			//F constants, force constants
			System.addTwoBodyFconst(RMIN);//C8
			System.addTwoBodyFconst((2.0*(Umax[g]-Umin[g]))/(RMIN*RMIN));//C6
			System.addTwoBodyFconst(0);//part of index trick
			System.addTwoBodyFconst(CUTOFF);//C7
			System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C3
			System.addTwoBodyFconst((6.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C2
			
			//U constants, potential constants
			System.addTwoBodyUconst(RMIN);//C8
			System.addTwoBodyUconst((Umax[g]-Umin[g])/(RMIN*RMIN));//C4
			System.addTwoBodyUconst(Umin[g]);//C5,no index trick
			System.addTwoBodyUconst(CUTOFF);//C7
			System.addTwoBodyUconst((3.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)));//C1
			System.addTwoBodyUconst((2.0*Umin[g])/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN)));//C0
		}
	}
	
	//Make it big
	threeVector<double> size;
	size.x=400;
	size.y=400;
	size.z=400;
	System.setSize(size);
	
	//just place it far away
	threeVector<double> pos;
	pos.x=size.x/2.0;
	pos.y=size.y/2.0;
	pos.z=size.z/2.0;
	
	double lipidConstants[4];
	lipidConstants[0]=0.7;
	lipidConstants[1]=100;
	lipidConstants[2]=1.0;
	lipidConstants[3]=100;
	
	//Add liposome
	//double radius=liposome(System, nLipids, 3, pos, 0.7, arialDensity, lipidConstants, 4);
	for(int i=0;i<nLipids;i++)
	{
		
	}
	
	///Write configuration
	Script<double,Blob <double> > output(name,ios::out,&System);
	output.write();
	
	delete Umin,Umax;
	return 0;
}
