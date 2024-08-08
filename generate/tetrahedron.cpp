#include "include.h"

int main(int argc, char **argv)
{
	if(argc<3)
	{
		std::cout << "Usage: " << argv[0] << " name seed\n";
		return 0;
	}

	char *name=argv[1];
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(2);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(1000);
	System.setDeltaT(0.02);
	System.setStoreInterval(0.02);
	System.setMeasureInterval(0.02);
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
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	threeVector<double>size;
	size.x=50;
	size.y=50;
	size.z=50;
	System.setSize(size);
	
	//Adding a tetrahedron
	MTRand randNum(System.readSeed());
	double Vrms=sqrt(3.0*System.readInitialTemp());
	fourVector<int> bond;
	molecule<double,fourVector<int> > m;
	position<double> vertex;
	vertex.type=1;
	threeVector<double> v,a;
	a=0;
	//in plane
	vertex.x=25;
	vertex.y=25;
	vertex.z=25;
	double theta=M_PI*randNum.rand53();
	double phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	vertex.x=25;
	vertex.y=26;
	vertex.z=25;
	theta=M_PI*randNum.rand53();
	phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	vertex.x=26;
	vertex.y=25;
	vertex.z=25;
	theta=M_PI*randNum.rand53();
	phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	//out of plane
	vertex.x=25.5;
	vertex.y=25.5;
	vertex.z=25.5;
	theta=M_PI*randNum.rand53();
	phi=M_PI*2.0*randNum.rand53();
	v.x=Vrms*cos(phi)*sin(theta);
	v.y=Vrms*sin(phi)*sin(theta);
	v.z=Vrms*cos(theta);
	System.addParticle(vertex,v,a);
	
	m.addConstant(1);
	m.addConstant(100);
	m.setType(BOND);
	bond.s[0]=0;
	bond.s[1]=1;
	m.addBond(bond);
	bond.s[0]=0;
	bond.s[1]=2;
	m.addBond(bond);
	bond.s[0]=0;
	bond.s[1]=3;
	m.addBond(bond);
	bond.s[0]=1;
	bond.s[1]=2;
	m.addBond(bond);
	bond.s[0]=1;
	bond.s[1]=3;
	m.addBond(bond);
	bond.s[0]=2;
	bond.s[1]=3;
	m.addBond(bond);
	System.addMolecule(m);
	
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
