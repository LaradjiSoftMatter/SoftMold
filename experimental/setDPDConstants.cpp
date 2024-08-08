/**
 * \brief Simple dissipative particle dynamic system. Just a bilayer in various geometries.
 * The steps to make a system are pretty simple, especially given the lipidModels.h header.
 * Just put together the parameters for initial conditions, generate the geometry, and
 * then run the system. That's it! If you want more, less, or something different, just
 * modify this program. There is a section that takes the command line arguments (just
 * after the main function definition) and checks that there are enough. After that,
 * atoi() and atof() (see some cstdlib documentation) convert the text to values. Those
 * are basically your initial conditions. The Blob class template is also useful for this.
 * Just take a parameter you want to modify and either use setParameter, addParameter, or
 * delParameter members in Blob to modify it. I've used the variable definition
 * 'Blob\<double\> System' for convienience. Functions in lipidModels.h accept System as
 * a reference, and can modify parameters as one would in main. For convienience, this
 * program outputs a name.xyz file to quickly view the initial geometry; this is useful
 * for debugging. Note that this is very similar to the setMDConstants.cpp program, but
 * it uses systemDPD.h instead. Also, solvent is now a requirement. It is recomended that
 * one generate the geometry first, then use the solventFill() function from lipidModels.h.
 */

#include <iostream>
#include <cstdlib>
#include "../include/systemDPD.h"
#include "../models/lipidModels.h"
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=7)
	{
		cout << "Usage: " << argv[0] << " name seed nLipids arealDensity density deltaT" << endl;
		return 0;
	}

	char *name=argv[1];
	int nLipids=atoi(argv[3]);
	double arealDensity=atof(argv[4]);
	double density=atof(argv[5]);
	
	///the variables for the simulation
	
	Blob<double> System;
	
	System.setNTypes(7);
	System.setSeed(atoi(argv[2]));
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(1.0);
	
	System.setInitialTime(0);
	System.setFinalTime(10);
	System.setDeltaT(atof(argv[6]));
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	System.setInitialTemp(1);
	System.setFinalTemp(1);
	
	//initialize constants (force and potential constants)
	//generate force constants
	double *chi=new double[System.readNTypes()*System.readNTypes()];//conservative
	double *gamma=new double[System.readNTypes()*System.readNTypes()];//dissipative
	double *sigma=new double[System.readNTypes()*System.readNTypes()];//random
	
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
			chi[k]=0;//25.0;
			sigma[k]=3.0;
			gamma[k]=(sigma[k]*sigma[k])/(2.0*System.readInitialTemp());
			lbond[k]=0.45;
			kbond[k]=100;
		}
	}
	
	//chi exceptions:
	chi[HEAD+TAIL*System.readNTypes()]=200;
	chi[TAIL+HEAD*System.readNTypes()]=chi[HEAD+TAIL*System.readNTypes()];
	
	chi[SOLVENT+TAIL*System.readNTypes()]=200;
	chi[TAIL+SOLVENT*System.readNTypes()]=chi[SOLVENT+TAIL*System.readNTypes()];
	
	//two body constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			int g=i+j*System.readNTypes();
			
			//fC[cindex]==cutoff
			//fC[cindex+1]==chi
			//fC[cindex+2]==gamma
			//fC[cindex+3]==sigma*deltaT^0.5
			System.addTwoBodyFconst(System.readCutoff());
			System.addTwoBodyFconst(chi[g]);
			System.addTwoBodyFconst(gamma[g]);
			System.addTwoBodyFconst(sigma[g]*sqrt(3.0/System.readDeltaT()));
			//System.addTwoBodyFconst(sqrt(3)*sigma[g]*sqrt(System.readDeltaT()));
		}
	}
	
	threeVector<double> size;
	size.x=sqrt((double)nLipids/arealDensity);
	size.y=sqrt((double)nLipids/arealDensity);
	size.z=10;
	System.setSize(size);
	threeVector<double> pos=size;
	pos.z/=2.0;

	double *constants=new double[2];
	constants[0]=0.45;
	constants[1]=100;
	
	bilayer<double>(System, nLipids, 4, pos, 0.45, arealDensity,constants,2,1.0);
	solventFill<double>(System, density, SOLVENT);
	
	///Write configuration
	Script<double,Blob <double> > output(name,ios::out,&System);
	output.write();
	output.close();
	
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
	
	delete chi,gamma,sigma,lbond,kbond,kbend,abend;
	return 0;
}


