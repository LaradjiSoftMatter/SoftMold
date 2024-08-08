#include "include.h"

#define NANOPARTICLE 4

int main(int argc, char **argv)
{
	if(argc<14)
	{
		std::cout << "Usage: " << argv[0] << " name seed nLipids arealDensity akbend bkbend iRatio oRatio lipidLength initialTemp finalTemp tempStepInterval finalTime" << std::endl;
		return 0;
	}

	char *name=argv[1];
	int nLipids=atoi(argv[3]);
	double arealDensity=atof(argv[4]);
	double akbend=atof(argv[5]);
	double bkbend=atof(argv[6]);
	double iRatio=atof(argv[7]);
	double oRatio=atof(argv[8]);
	int lipidLength=atoi(argv[9]);

	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(5);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(atof(argv[13]));
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(atof(argv[10]));
	System.setFinalTemp(atof(argv[11]));
	System.setTempStepInterval(atof(argv[12]));
	
	//initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	
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
	Umin[TAILA+TAILA*System.readNTypes()]=-6;
	Umax[TAILA+TAILA*System.readNTypes()]=200;
	
	Umin[TAILB+TAILB*System.readNTypes()]=-6;
	Umax[TAILB+TAILB*System.readNTypes()]=200;
	
	Umin[TAILA+TAILB*System.readNTypes()]=-6;
	Umax[TAILA+TAILB*System.readNTypes()]=200;
	
	Umin[TAILB+TAILA*System.readNTypes()]=Umin[TAILA+TAILB*System.readNTypes()];
	Umax[TAILB+TAILA*System.readNTypes()]=Umax[TAILA+TAILB*System.readNTypes()];
	
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
	
	double radius=sqrt(((double)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=radius*4.0;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=System.readSize().x/2.0;
	pos.y=System.readSize().y/2.0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *aConstants=new double[4];
	double *bConstants=new double[4];
	aConstants[0]=bondLength;
	aConstants[1]=100;
	aConstants[2]=1;
	aConstants[3]=akbend;
	
	bConstants[0]=bondLength;
	bConstants[1]=100;
	bConstants[2]=1;
	bConstants[3]=bkbend;
	
	int *types=new int[4];
	types[0]=TAILA;
	types[1]=HEADA;
	types[2]=TAILB;
	types[3]=HEADB;
	
	//liposome with cytoskeleton
	radius=liposome<double>(System, nLipids, lipidLength, pos, 0.7, arealDensity, iRatio, oRatio, types, aConstants, bConstants, 4);
	
	//add solvent
	//solventFill<double>(System, 2.0, SOLVENT_FLAG);
	
	//Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	std::string newName("");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	delete aConstants,bConstants;
	delete types;
	return 0;
}


