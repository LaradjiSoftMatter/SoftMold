#include "include.h"

#define NANOPARTICLE 4

int main(int argc, char **argv)
{
	if(argc!=4)
	{
		std::cout << "Usage: " << argv[0] << " nLipids arial_density name" << std::endl;
		return 0;
	}

	char *name=argv[3];
	
	//load parameters
	int nLipids=atoi(argv[1]);
	double arialDensity=atof(argv[2]);
	
	///the variables for the simulation
	
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(5);
	System.setSeed(1234);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	threeVector<double> size;   //component sizes of system
	size.x=sqrt(nLipids/arialDensity);//pow((double)particles/density,1.0/3.0);
	size.y=size.x;
	size.z=100.0;//size.x;
	System.setSize(size);
	
	System.setInitialTime(0);
	System.setFinalTime(10);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	
	threeVector<int> latticeSize;
	threeVector<double> latticeLength;
	
	latticeSize.x=sqrt((double)nLipids/2.0)+1;//int(pow((double)particles,1.0/3.0))+1;
	latticeSize.y=latticeSize.x;
	latticeSize.z=latticeSize.y;
	
	latticeLength.x=size.x/latticeSize.x;
	latticeLength.y=size.y/latticeSize.y;
	latticeLength.z=size.z/latticeSize.z;
	
	//creating a bunch of molecules with specialized properties
	
	fourVector<int> bond;
	molecule<double,fourVector<int> > m;
	
	m.setType(CHAIN);
	//constants for bond
	m.addConstant(0.7);
	m.addConstant(100);
	m.addConstant(1.0);
	m.addConstant(100);
	//For a chain type bond
	bond.s[START]=0;
	bond.s[NCHAINS]=nLipids;
	bond.s[CHAINLENGTH]=3;
	m.addBond(bond);
	
	System.addMolecule(m);
	
	m.~molecule();
	
	MTRand *randNum=new MTRand(System.readSeed());
	
	//Velocities
	double Vrms=sqrt(3*System.readInitialTemp());
	
	///This is important to reduce excessive allocations
	System.allocParticle(System.readNParticles()+nLipids*3);
	
	//initialize lipid positions
	for(int i=0;i<latticeSize.x;i++)
	{
		for(int j=0;j<latticeSize.y;j++)
		{
			int lipidIndex=2*(j*latticeSize.x+i);
			double theta,phi;
			
			position<double> p;
			threeVector<double> v;
			threeVector<double> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			if(lipidIndex<nLipids)
			{
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=0.001+size.z/2.0-2.501;
				p.type=HEAD;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=1.001+size.z/2.0-2.501;//+latticeLength.z;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=2.001+size.z/2.0-2.501;//+latticeLength.z*2.0;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
			}
			if(lipidIndex+1<nLipids)
			{
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=3.001+size.z/2.0-2.501;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=4.001+size.z/2.0-2.501;//+latticeLength.z;
				p.type=TAIL;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
				
				//position
				p.x=i*latticeLength.x+0.001;
				p.y=j*latticeLength.y+0.001;
				p.z=5.001+size.z/2.0-2.501;//+latticeLength.z*2.0;
				p.type=HEAD;
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				//save it
				System.addParticle(p,v,a);
			}
		}
	}
	
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
	
	Umin[NANOPARTICLE+NANOPARTICLE*System.readNTypes()]=0;
	Umax[NANOPARTICLE+NANOPARTICLE*System.readNTypes()]=0;
	
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
	
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	
	output.write();
	
	delete randNum;
	delete lbond,kbond,kbend,abend;
	return 0;
}


