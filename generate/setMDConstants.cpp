#include "include.h"

int main(int argc, char **argv)
{
	if(argc<7)
	{
		//7, bilayer, nanoparticle
		std::cout << "Usage: " << argv[0] << " name seed anchorHeadUmin nLipids lipidArealDensity nanoRadius nanoHeadUmin\n";
		return 0;
	}

	char *name=argv[1];
	double anchorHeadUmin=atof(argv[3]);
	int nLipids=atoi(argv[4]);
	double arealDensity=atof(argv[5]);
	double radius=atof(argv[6]);
	double UminNanoHead=atof(argv[7]);
	
	//the variables for the simulation
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
	System.setFinalTime(100);
	System.setDeltaT(0.01);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	//System.setSolventGamma(5.0);
	
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
	
	Umin[4+HEAD*System.readNTypes()]=UminNanoHead;
	Umax[4+HEAD*System.readNTypes()]=200;
	
	Umin[HEAD+4*System.readNTypes()]=Umin[4+HEAD*System.readNTypes()];
	Umax[HEAD+4*System.readNTypes()]=Umax[4+HEAD*System.readNTypes()];
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	threeVector<double>size;
	size=0;
	size.z=40;
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
	
	//Add a nanoparticle
	int nanoType=4;
	pos.x=System.readSize().x/2.0;
	pos.y=System.readSize().y/2.0;
	pos.z+=radius+1.0+3.0*2.0*0.7;
	
	position<double> *unitCell=new position<double>[4];
	unitCell[0].x=0;
	unitCell[0].y=0;
	unitCell[0].z=0;
	unitCell[0].type=nanoType;
	unitCell[1].x=0.5;
	unitCell[1].y=0.5;
	unitCell[1].z=0;
	unitCell[1].type=nanoType;
	unitCell[2].x=0.5;
	unitCell[2].y=0;
	unitCell[2].z=0.5;
	unitCell[2].type=nanoType;
	unitCell[3].x=0;
	unitCell[3].y=0.5;
	unitCell[3].z=0.5;
	unitCell[3].type=nanoType;
	
	position<double> *bonded=new position<double>[12];
	bonded[0].x=0.5;
	bonded[0].y=0.5;
	bonded[0].z=0;
	bonded[0].type=nanoType;
	bonded[1].x=0.5;
	bonded[1].y=0;
	bonded[1].z=0.5;
	bonded[1].type=nanoType;
	bonded[2].x=0;
	bonded[2].y=0.5;
	bonded[2].z=0.5;
	bonded[2].type=nanoType;
	bonded[3].x=-0.5;
	bonded[3].y=-0.5;
	bonded[3].z=0;
	bonded[3].type=nanoType;
	bonded[4].x=-0.5;
	bonded[4].y=0;
	bonded[4].z=-0.5;
	bonded[4].type=nanoType;
	bonded[5].x=0;
	bonded[5].y=-0.5;
	bonded[5].z=-0.5;
	bonded[5].type=nanoType;
	bonded[6].x=0.5;
	bonded[6].y=-0.5;
	bonded[6].z=0;
	bonded[6].type=nanoType;
	bonded[7].x=0.5;
	bonded[7].y=0;
	bonded[7].z=-0.5;
	bonded[7].type=nanoType;
	bonded[8].x=0;
	bonded[8].y=0.5;
	bonded[8].z=-0.5;
	bonded[8].type=nanoType;
	bonded[9].x=-0.5;
	bonded[9].y=0.5;
	bonded[9].z=0;
	bonded[9].type=nanoType;
	bonded[10].x=-0.5;
	bonded[10].y=0;
	bonded[10].z=0.5;
	bonded[10].type=nanoType;
	bonded[11].x=0;
	bonded[11].y=-0.5;
	bonded[11].z=0.5;
	bonded[11].type=nanoType;
	
	constants[0]=sqrt(2.0)/4.0;//bondlength=latticeLength*sqrt(2.0)/2.0 where latticeLength=0.5
	constants[1]=2800;//kbond
	
	nanoSphere<double>(System, pos, radius, 0.5*1.31345933134, unitCell, 4, bonded, 12, constants, 2);
	
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
	delete bonded,unitCell;
	return 0;
}


