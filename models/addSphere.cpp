//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#include "../include/systemMD.h"
#include <ctime>
#define THRESHOLD 1.0

int main(int argc, char* argv[])
{
	if(argc!=7)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name nLipids posX posY posZ density\n";
		return 0;
	}
	
	threeVector<double> pos,norm;
	char *name=argv[1];
	int nLipids=atoi(argv[2]);
	pos.x=atoi(argv[3]);
	pos.y=atoi(argv[4]);
	pos.z=atoi(argv[5]);
	double density=atof(argv[6]);
	
	///Previous Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	std::string newName("");
	newName+=name;
	newName+="_wL";

	///This allocates new variables
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
	bond.s[LENGTH]=3;
	m.addBond(bond);
	System.addMolecule(m);
	
	m.~molecule();
	///This initializes the new variables
	//copypasta from an earlier version
	
	density=sqrt((double)nLipids/density);//sqrt(1/p)
	double radius=sqrt((density*density)/(4.0*M_PI));//density adjusted
	
	//adjust size of system if it is out of bounds
	threeVector<double> size=System.readSize();
	size.x=size.x<radius+pos.x+5.0?radius+pos.x+6.0:size.x;
	size.y=size.y<radius+pos.y+5.0?radius+pos.y+6.0:size.y;
	size.z=size.z<radius+pos.z+5.0?radius+pos.z+6.0:size.z;
	System.setSize(size);
	
	//with inner and outer surface area compensation
	int inner=((4.0*M_PI*radius*radius)/((4.0*M_PI*radius*radius)+(4.0*M_PI*(radius+5.0*1.0)*(radius+5.0*1.0))))*bond.s[NCHAINS];
	int outer=bond.s[NCHAINS]-inner;
	double s=3.6/sqrt(double(inner));
	double length=0;
	double dz=2.0/double(inner);
	double z=1.0-dz/2.0;
	
	std::cout << "Radius: " << radius << '\n';
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	double Vrms=sqrt(3*System.readInitialTemp());
	
	///This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+nLipids*3);
	
	for(int i=0;i<inner*bond.s[LENGTH];i+=bond.s[LENGTH])
	{
		double phi,theta;
		position <double> p;
		threeVector<double> v,a;
		//a is just a null vector
		a.x=0;
		a.y=0;
		a.z=0;
		
		double r=sqrt(1.0-z*z);
		//inner monomers to outer monomers
		p.x=pos.x+r*cos(length)*radius;
		p.y=pos.y+r*sin(length)*radius;
		p.z=pos.z+z*radius;
		p.type=HEAD;
		//velocity
		theta=M_PI*randNum->rand53();
		phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
		
		p.x=pos.x+r*cos(length)*(radius+1.0);
		p.y=pos.y+r*sin(length)*(radius+1.0);
		p.z=pos.z+z*(radius+1.0);
		p.type=TAIL;
		//velocity
		theta=M_PI*randNum->rand53();
		phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
		
		p.x=pos.x+r*cos(length)*(radius+1.0*2.0);
		p.y=pos.y+r*sin(length)*(radius+1.0*2.0);
		p.z=pos.z+z*(radius+1.0*2.0);
		p.type=TAIL;
		//velocity
		theta=M_PI*randNum->rand53();
		phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
		
		z=z-dz;
		length=length+s/r;
	}

	
	s=3.6/sqrt(double(outer));
	length=0;
	dz=2.0/double(outer);
	z=1.0-dz/2.0;
	
	for(int i=inner*bond.s[LENGTH];i<(outer+inner)*bond.s[LENGTH];i=i+bond.s[LENGTH])
	{
		double phi,theta;
		position <double> p;
		threeVector<double> v,a;
		//a is just a null vector
		a.x=0;
		a.y=0;
		a.z=0;
		
		double r=sqrt(1.0-z*z);
		//outer monomers to inner monomers
		p.x=pos.x+r*cos(length)*(radius+1.0*5.0);
		p.y=pos.y+r*sin(length)*(radius+1.0*5.0);
		p.z=pos.z+z*(radius+1.0*5.0);
		p.type=HEAD;
		//velocity
		theta=M_PI*randNum->rand53();
		phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
		
		p.x=pos.x+r*cos(length)*(radius+1.0*4.0);
		p.y=pos.y+r*sin(length)*(radius+1.0*4.0);
		p.z=pos.z+z*(radius+1.0*4.0);
		p.type=TAIL;
		//velocity
		theta=M_PI*randNum->rand53();
		phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
		
		p.x=pos.x+r*cos(length)*(radius+1.0*3.0);
		p.y=pos.y+r*sin(length)*(radius+1.0*3.0);
		p.z=pos.z+z*(radius+1.0*3.0);
		p.type=TAIL;
		//velocity
		theta=M_PI*randNum->rand53();
		phi=M_PI*2*randNum->rand53();
		v.x=Vrms*cos(phi)*sin(theta);
		v.y=Vrms*sin(phi)*sin(theta);
		v.z=Vrms*cos(theta);
		System.addParticle(p,v,a);
		
		z=z-dz;
		length=length+s/r;
	}
	
	
	///This stores the new configuration
	Script<double, Blob <double> > fileIO_new(newName.c_str(),std::ios::out,&System);
	fileIO_new.write();
	fileIO_new.close();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p, nParticles);
	newName+=".xyz";
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	return 0;
}

