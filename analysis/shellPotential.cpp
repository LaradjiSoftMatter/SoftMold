/**
 * \brief Creates a bilayer membrane with nanoparticles.
 * Nanoparticle chain placement.
 */

#include <iostream>
#include <cstdlib>
#include <map>
#include <list>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include "../models/lipidModels.h"
#include "../models/nanoModels.h"

#define HEAD 2
#define TAIL 3

int main(int argc, char **argv)
{
	if(argc!=10)
	{
		std::cerr << "Usage: " << argv[0] << " Umax Umin rc rmin nanoRadius nTesselations dr trials seed\n";
		return 0;
	}
	std::stringstream cmdArg;
	double UmaxNano, UminNano, rc, rmin, nanoRadius, density, dr;
	int trials,seed, nTesselations;
	
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> UmaxNano >> UminNano >> rc >> rmin >> nanoRadius >> nTesselations >> dr >> trials >> seed;
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(6);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(100);
	System.setDeltaT(0.01);
	System.setStoreInterval(100);
	System.setMeasureInterval(10);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	//System.setSolventGamma(5.0);
	threeVector<double>size;
	size=100;
	System.setSize(size);
	
	//initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	
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
	
	int nanoType=4;
	
	Umin[nanoType+HEAD*System.readNTypes()]=UminNano;
	Umax[nanoType+HEAD*System.readNTypes()]=UmaxNano;
	
	Umin[HEAD+nanoType*System.readNTypes()]=Umin[nanoType+HEAD*System.readNTypes()];
	Umax[HEAD+nanoType*System.readNTypes()]=Umax[nanoType+HEAD*System.readNTypes()];
	
	//two body constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			int g=i+j*System.readNTypes();
			
			//F constants, force constants
			System.addTwoBodyFconst(rmin);//C8
			System.addTwoBodyFconst((2.0*(Umax[g]-Umin[g]))/(rmin*rmin));//C6
			System.addTwoBodyFconst(0);//part of index trick
			System.addTwoBodyFconst(rc);//C7
			System.addTwoBodyFconst((6.0*Umin[g])/((rc-rmin)*(rc-rmin)));//C3
			System.addTwoBodyFconst((6.0*Umin[g])/((rc-rmin)*(rc-rmin)*(rc-rmin)));//C2
			
			//U constants, potential constants
			System.addTwoBodyUconst(rmin);//C8
			System.addTwoBodyUconst((Umax[g]-Umin[g])/(rmin*rmin));//C4
			System.addTwoBodyUconst(Umin[g]);//C5,no index trick
			System.addTwoBodyUconst(rc);//C7
			System.addTwoBodyUconst((3.0*Umin[g])/((rc-rmin)*(rc-rmin)));//C1
			System.addTwoBodyUconst((2.0*Umin[g])/((rc-rmin)*(rc-rmin)*(rc-rmin)));//C0
		}
	}
	
	//std::vector<double> constants(2);
	//constants[0]=sqrt(2.0)/4.0;//bondlength=latticeLength*sqrt(2.0)/2.0 where latticeLength=0.5
	//constants[1]=2800;//kbond
	
	threeVector<double> nanoPos=System.readSize();
	nanoPos.x/=2.0;
	nanoPos.y/=2.0;
	nanoPos.z/=2.0;
	MTRand randNum(System.readSeed());
	
	int start=System.readNParticles();
	
	//hollow np
	std::cerr << "Generating a hollow nanoparticle!" << std::endl;
	shellNano(System, nanoPos, nanoRadius, nTesselations, nanoType);
	density=(System.readNParticles()-start)/(4.0*M_PI*nanoRadius*nanoRadius);
	std::cerr << "density = " << density << std::endl;
	int end=System.readNParticles();
	
	position<double> *p=System.getPositions();
	
	threeVector<double> nanoCom=0;
	for(int i=start;i<end;i++)
	{
		nanoCom.x+=p[i].x;
		nanoCom.y+=p[i].y;
		nanoCom.z+=p[i].z;
	}
	if(end-start!=0)
	{
		nanoCom.x/=(end-start);
		nanoCom.y/=(end-start);
		nanoCom.z/=(end-start);
	}
	
	double cutoffSquared=rc*rc;
	std::map<int,std::list<double>> potentials;
	
	double *pConst=System.getTwoBodyUconst();
	int cindex=nTWOBODYUCONST*((nanoType*System.readNTypes())+HEAD);
	
	for(int i=0;i<trials;i++)
	{
		if(i%(trials/10)==0)
			std::cerr << i << " trials" << std::endl;
		double theta=M_PI*randNum.rand53();
		double phi=M_PI*2.0*randNum.rand53();
		double radius=nanoRadius+randNum.rand53()*rc;
		
		position<double> testP;
		testP.x=radius*cos(phi)*sin(theta)+nanoCom.x;
		testP.y=radius*sin(phi)*sin(theta)+nanoCom.y;
		testP.z=radius*cos(theta)+nanoCom.z;
		testP.type=HEAD;
		
		double totalPotential=0;
		#pragma omp parallel for reduction(+:totalPotential)
		for(int l=start;l<end;l++)
		{
			threeVector<double> d;
			d.x=p[l].x-testP.x;
			d.y=p[l].y-testP.y;
			d.z=p[l].z-testP.z;
			totalPotential+=nonBondedP<double>(d, &pConst[cindex], cutoffSquared);
		}
		
		int binPos=floor(radius/dr);
		potentials[binPos].push_back(totalPotential);
	}
	for(auto potentialBin:potentials)
	{
		double radius=(dr*potentialBin.first)+dr/2.0;
		double avgPotential=0;
		for(auto potential:potentialBin.second)
			avgPotential+=potential;
		if(potentialBin.second.size()>0)
			avgPotential/=potentialBin.second.size();
		std::cout << radius << '\t' << avgPotential << std::endl;
	}
	
	//std::cerr << "Storing configuration...";
	
	//Write configuration
	//Script<double,Blob <double> > output(name,ios::out,&System);
	//output.write();
	
	//This stores the xyz file to check new configuration
	int nParticles=System.readNParticles();
	
	std::string newName("");
	newName+="out";
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	std::cerr << "Done.\nExiting...\n";
	
	return 0;
}


