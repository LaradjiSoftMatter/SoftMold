#define SOLVENT_FLAG 0

#include "include.h"

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		//6
		std::cout << "Usage: " << argv[0] << " name seed solventDensity nLipids arealDensity\n";
		return 0;
	}

	char *name=argv[1];
	double solventDensity=atof(argv[3]);
	int nLipids=atoi(argv[4]);
	double arealDensity=atof(argv[5]);
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(4);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(10000);
	System.setDeltaT(0.02);
	System.setStoreInterval(1000);
	System.setMeasureInterval(100);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
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
	
	//tail types hydrophobicity
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	threeVector<double>size;
	size=0;//radius*4.0;
	size.z=20;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;//System.readSize().x/2.0;
	pos.y=0;//System.readSize().y/2.0;
	pos.z=10;//System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//bilayer with cytoskeleton, regular hexagonal cytoskeleton has an aspect ratio sqrt(3/2)=y/x
	bilayer<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants, 4, 1.0);
	
	//add solvent
	solventFill<double>(System, solventDensity, SOLVENT_FLAG);
	
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
	
	delete constants;
	return 0;
}


