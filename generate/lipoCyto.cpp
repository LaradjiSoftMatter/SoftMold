#include "include.h"

int main(int argc, char **argv)
{
	if(argc<10)
	{
		std::cout << "Usage: " << argv[0] << " name seed Umin_Anchor_Head fraction_solvent_removed nLipids arialDensity density nMonomers nTess" << std::endl;
		return 0;
	}

	char *name=argv[1];
	double removeSolventFraction=atof(argv[4]);
	int nLipids=atoi(argv[5]);
	double arialDensity=atof(argv[6]);
	double density=atof(argv[7]);
	int nMonomers=atoi(argv[8]);
	int nTess=atoi(argv[9]);
	
	
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
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(50000);
	System.setDeltaT(0.02);
	System.setStoreInterval(100);
	System.setMeasureInterval(10);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
	///initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	double *lbond=new double[System.readNTypes()*System.readNTypes()];
	double *kbond=new double[System.readNTypes()*System.readNTypes()*System.readNTypes()];
	double *abend=new double[System.readNTypes()*System.readNTypes()*System.readNTypes()];
	double *kbend=new double[System.readNTypes()*System.readNTypes()*System.readNTypes()];
	
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
			lbond[k]=0.7;
			kbond[k]=100;
		}
	}
	
	//umin and umax exceptions, tail types
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	//for polymer strand interactions
	Umin[MONOMER+MONOMER*System.readNTypes()]=0;
	Umax[MONOMER+MONOMER*System.readNTypes()]=100;
	
	Umin[HEAD+MONOMER*System.readNTypes()]=0;
	Umax[HEAD+MONOMER*System.readNTypes()]=100;
	
	Umin[MONOMER+HEAD*System.readNTypes()]=Umin[HEAD+MONOMER*System.readNTypes()];
	Umax[MONOMER+HEAD*System.readNTypes()]=Umax[HEAD+MONOMER*System.readNTypes()];

	Umin[TAIL+MONOMER*System.readNTypes()]=0;
	Umax[TAIL+MONOMER*System.readNTypes()]=100;

	Umin[MONOMER+TAIL*System.readNTypes()]=0;
	Umax[MONOMER+TAIL*System.readNTypes()]=100;
	
	Umin[HEAD+ANCHOR*System.readNTypes()]=atof(argv[3]);
	Umax[HEAD+ANCHOR*System.readNTypes()]=100;
	
	Umin[ANCHOR+HEAD*System.readNTypes()]=Umin[HEAD+ANCHOR*System.readNTypes()];
	Umax[ANCHOR+HEAD*System.readNTypes()]=Umax[HEAD+ANCHOR*System.readNTypes()];
	
//	Umin[HEAD+CYTO*System.readNTypes()]=-6;
//	Umax[HEAD+CYTO*System.readNTypes()]=200;
	
//	Umin[CYTO+HEAD*System.readNTypes()]=Umin[HEAD+CYTO*System.readNTypes()];
//	Umax[CYTO+HEAD*System.readNTypes()]=Umax[HEAD+CYTO*System.readNTypes()];
	
	//how to do it with one loop
	for(int i=0;i<System.readNTypes()*System.readNTypes()*System.readNTypes();i++)
	{
		kbend[i]=100;
		abend[i]=1;
	}
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//just place it far away
	threeVector<double> pos;
	pos.x=100;
	pos.y=100;
	pos.z=100;
	
	double lipidConstants[4];
	lipidConstants[0]=0.7;
	lipidConstants[1]=100;
	lipidConstants[2]=1.0;
	lipidConstants[3]=100;
	
	double cytoBondLength=1.0;
	double cytoConstants[4];
	cytoConstants[0]=cytoBondLength;
	cytoConstants[1]=100;
	cytoConstants[2]=1.0;
	cytoConstants[3]=100;
	
	//Add liposome
	double radius=liposome(System, nLipids, 3, pos, 0.7, arialDensity, lipidConstants, 4);
	
	//Add cytoskeleton
	sphericalCyto(System, nMonomers, pos, radius, nTess, cytoConstants, 4);
	
	///Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	delete lbond,kbond,kbend,abend;
	return 0;
}


