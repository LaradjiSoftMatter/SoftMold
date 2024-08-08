#include "include.h"

int main(int argc, char **argv)
{
	if(argc<9)
	{
		//9, bilayer
		std::cout << "Usage: " << argv[0] << " name seed nLipids arealDensity UminAnchorHead nAnchorX nAnchorY nMonomers" << std::endl;
		return 0;
	}

	char *name=argv[1];
	int seed=atoi(argv[2]);
	int nLipids=atoi(argv[3]);
	double arealDensity=atof(argv[4]);
	double UminAnchorHead=atof(argv[5]);
	twoVector<int> nAnchors;
	nAnchors.x=atoi(argv[6]);
	nAnchors.y=atoi(argv[7]);
	int nMonomers=atoi(argv[8]);
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(6);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(1000);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
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
	
	Umin[ANCHOR+HEAD*System.readNTypes()]=UminAnchorHead;
	Umax[ANCHOR+HEAD*System.readNTypes()]=200;
	
	Umin[HEAD+ANCHOR*System.readNTypes()]=Umin[ANCHOR+HEAD*System.readNTypes()];
	Umax[HEAD+ANCHOR*System.readNTypes()]=Umax[ANCHOR+HEAD*System.readNTypes()];
	
	//Umin[TAIL2+TAIL2*System.readNTypes()]=-6;
	//Umax[TAIL2+TAIL2*System.readNTypes()]=200;
	
	//Umin[TAIL2+TAIL*System.readNTypes()]=-5.5;
	//Umax[TAIL2+TAIL*System.readNTypes()]=200;
	
	//Umin[TAIL+TAIL2*System.readNTypes()]=-5.5;
	//Umax[TAIL+TAIL2*System.readNTypes()]=200;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//double radius=sqrt(((double)nLipids/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;
	size.z=40;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;
	pos.y=0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//simple bilayer
	if(nAnchors.x>0 && nAnchors.y>0)
	{
		bilayer<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants, 4, (double(nAnchors.x)/double(nAnchors.y))/HEXAGONAL_ASPECT_RATIO);
		pos.z-=bondLength;
		flatHexagonalCyto<double>(System, nAnchors, nMonomers, pos, constants, 4);
	}
	else
	{
		bilayer<double>(System, nLipids, 3, pos, bondLength, arealDensity, constants, 4, 1.0);
	}

	std::cout << "Storing configuration...";
	
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
	
	std::cout << "Done.\nExiting...\n";
	
	delete constants;
	return 0;
}


