
#define SOLVENT_FLAG 0
#include "include.h"

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		std::cout << "Usage: " << argv[0] << " name seed solventDensity gamma\n";
		return 0;
	}

	char *name=argv[1];
	double solventDensity=atof(argv[3]);
	double gamma=atof(argv[4]);
	
	//the variables for the simulation
	Blob<double> System;
	
//	System.setGamma(gamma);
	System.setNTypes(2);
	System.setGammaType(0,gamma);
	System.setGammaType(1,1.0);
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
	System.setMeasureInterval(1);
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
	
	//Solvent Types don't interact
	Umin[SOLVENT_FLAG+SOLVENT_FLAG*System.readNTypes()]=0;
	Umax[SOLVENT_FLAG+SOLVENT_FLAG*System.readNTypes()]=0;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	threeVector<double>size;
	size=20;//radius*4.0;
	//size.z=20;
	System.setSize(size);
	
	//add solvent
	solventFill<double>(System, solventDensity, SOLVENT_FLAG);
	
	//change one particle to type=1
	//System.getPositions()[System.readNParticles()/2].type=SOLVENT_TRACER;
	
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
	
	return 0;
}
