//Adjusts some constants. Just uses Blob and script to read, modify, and write values.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <ctime>

int main(int argc, char* argv[])
{
	if(argc!=6)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name type1 type2 Umax Umin\n";
		return 0;
	}
	
	char *name=argv[1];
	
	std::stringstream cmdStr;
	for(int i=2;i<argc;i++)
		cmdStr << argv[i] << ' ';
	int type1,type2;
	double Umax,Umin;
	cmdStr >> type1 >> type2 >> Umax >> Umin;
	
	//the variables for the simulation
	Blob<double> System;
	
	//load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	if(type1>=System.readNTypes() || type2>=System.readNTypes() || type1<0 || type2<0)
	{
		std::cout << "Error: Types don't exist, valid type range: 0 to " << System.readNTypes()-1 << '\n';
		return 0;
	}
	
	if(Umax<0)
		std::cout << "Warning: Not sure if Umax should be less than 0!\n";
	if(Umin>0)
		std::cout << "Warning: Not sure if Umin should be greater than 0!\n";
	
	//New values:
	int cindexF1=nTWOBODYFCONST*((type1*System.readNTypes())+type2);
	int cindexP1=nTWOBODYUCONST*((type1*System.readNTypes())+type2);
	int cindexF2=nTWOBODYFCONST*((type2*System.readNTypes())+type1);
	int cindexP2=nTWOBODYUCONST*((type2*System.readNTypes())+type1);
	
	//These are usually fixed
	double RMIN=1.0;
	double CUTOFF=2.0;
	
	System.getTwoBodyFconst()[cindexF1]=RMIN;
	System.getTwoBodyFconst()[cindexF1+1]=(2.0*(Umax-Umin))/(RMIN*RMIN);
	System.getTwoBodyFconst()[cindexF1+2]=0;
	System.getTwoBodyFconst()[cindexF1+3]=CUTOFF;
	System.getTwoBodyFconst()[cindexF1+4]=(6.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN));
	System.getTwoBodyFconst()[cindexF1+5]=(6.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN));
	
	System.getTwoBodyUconst()[cindexP1]=RMIN;
	System.getTwoBodyUconst()[cindexP1+1]=(Umax-Umin)/(RMIN*RMIN);
	System.getTwoBodyUconst()[cindexP1+2]=Umin;
	System.getTwoBodyUconst()[cindexP1+3]=CUTOFF;
	System.getTwoBodyUconst()[cindexP1+4]=(3.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN));
	System.getTwoBodyUconst()[cindexP1+5]=(2.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN));
	
	System.getTwoBodyFconst()[cindexF2]=RMIN;
	System.getTwoBodyFconst()[cindexF2+1]=(2.0*(Umax-Umin))/(RMIN*RMIN);
	System.getTwoBodyFconst()[cindexF2+2]=0;
	System.getTwoBodyFconst()[cindexF2+3]=CUTOFF;
	System.getTwoBodyFconst()[cindexF2+4]=(6.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN));
	System.getTwoBodyFconst()[cindexF2+5]=(6.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN));
	
	System.getTwoBodyUconst()[cindexP2]=RMIN;
	System.getTwoBodyUconst()[cindexP2+1]=(Umax-Umin)/(RMIN*RMIN);
	System.getTwoBodyUconst()[cindexP2+2]=Umin;
	System.getTwoBodyUconst()[cindexP2+3]=CUTOFF;
	System.getTwoBodyUconst()[cindexP2+4]=(3.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN));
	System.getTwoBodyUconst()[cindexP2+5]=(2.0*Umin)/((CUTOFF-RMIN)*(CUTOFF-RMIN)*(CUTOFF-RMIN));
	
	//How it was done originally (in setMDConstants.cpp):
	/*
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
	*/
	
	fileIO.open(name,std::ios::out);
	fileIO.write();
	fileIO.close();
		
	return 0;
}
