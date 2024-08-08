//This simply runs a system through the timesteps. It can use MPI.

#include "include/system.h"
#include <ctime>
#include <mpi.h>

int main(int argc, char* argv[])
{
	int myrank;
	MPI_Init(&argc, &argv);
	
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name\n";
		return 0;
	}
	
	char *name=argv[1];
	
	
	///the variables for the simulation
	//physical part
	position<double> *p=NULL;
	threeVector<double> *v=NULL;
	threeVector<double> *a=NULL;
	molecule<double> *m=NULL;
	
	//constants loaded at run time
	double gamma;
	double tempI;
	double tempF;
	double *twoBodyFconst=NULL;
	double *twoBodyUconst=NULL;
	double *twoBodyMFconst=NULL;
	double *twoBodyMUconst=NULL;
	double *threeBodyMFconst=NULL;
	double *threeBodyMUconst=NULL;
	
	//basic properties of the particle system
	int seed;      //random number seed
	int types;     //number of particle types
	int molecules; //number of molecules
	int particles; //number of particles
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	double cutoff;              //cutoff for algorithms
	threeVector<double> size;   //component sizes of system
	double ti;              //initial time
	double end;             //end of run
	double deltat;          //time step size, depends on algorithm in process
	double storeInterval;   //period for storing data for continuation, float?
	double measureInterval; //measure and gather data at time steps, float?
	
	
	///load variables, then initialize them
	
	Blob<double> System(p, v, a, m, gamma, tempI, tempF, twoBodyFconst, twoBodyUconst, twoBodyMFconst,
		twoBodyMUconst, threeBodyMFconst, threeBodyMUconst, seed, types, molecules, particles, 
		wrap, cutoff, size, ti, end, deltat, storeInterval, measureInterval);
	
	Script<double,Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	//initialize the variables, objects, algorithms, data collection, etc...
	double tempStep=(storeInterval*(tempF-tempI))/(end-ti);
	
	Verlet<double> integrate(p, a, v, particles, size, deltat, wrap);
	Langevin<double> thermostat(a, v, particles, gamma, deltat, seed);
	Cell<double, Potential<double>, Force <double> > pairInteractions(p, a, twoBodyFconst,
		twoBodyUconst, particles, types, size, wrap, cutoff);
	Molecules<double, HarmonicPotential<double>, HarmonicForce<double>, BendPotential<double>,
			  BendForce<double> > moleculeInteractions(p, a,
			  twoBodyMFconst, twoBodyMUconst, threeBodyMFconst, threeBodyMUconst, NULL,
			  NULL, NULL, NULL, m, particles, types, molecules, size);
	
	//initialize data collection objects
	Kinetic<double> kinetic(v, particles);
	
	std::string framesFilename("frames_");
	framesFilename+=argv[1];
	framesFilename+=".xyz";
	xyzFormat<double> xyzFile(p, particles);
	
	///this is how you call the pair interactions, build then compute...
	pairInteractions.build();
	
	pairInteractions.computeForce();
	moleculeInteractions.computeForce();
	
	//using integer indexing, the 0.0000001 fixes an accuracy issue with the gnu c++ compiler.
	//Don't believe it affects anything else...
	int endInt=int(end/deltat+0.0000001);//end
	int startInt=int(ti/deltat+0.0000001);//start
	
	int storeint=int(storeInterval/deltat+0.0000001);//when to store
	int measureint=int(measureInterval/deltat+0.0000001);//when to measure
	
	//this corrects an issue where an extra data point is added when the system is restarted
	if(ti==0)
	{
		//I'm too lazy to put a proper potential calculation here, it's only one step out of millions..., but it's close
		pairInteractions.outputPotential(ti,argv[1]);
		kinetic.output(ti,argv[1]);
		xyzFile.open(framesFilename,std::ios::out | std::ios::app);
		xyzFile.store();
		xyzFile.close();
	}
	
	
	
	std::cout << "starting main loop: \n";
	time_t current=time(NULL);
	
	///the molecular dynamics loop, the "running" of the system
	for(int i=startInt;i<=endInt;i++)
	{
		ti=i*deltat;
		integrate.first();
		
		for(int k=0;k<particles;k++)
		{
			a[k].x=0;
			a[k].y=0;
			a[k].z=0;
		}
		thermostat.compute(tempI);
		pairInteractions.build();
		pairInteractions.computeForce();
		moleculeInteractions.computeForce();
		integrate.second();
		
		if(i%storeint==0 && i!=startInt)
		{
			//double last=current;
			time_t last=current;
			//current=omp_get_wtime();
			current=time(NULL);
			//time since last storage step, good for benchmarking
			std::cout << ti << '\t' << current-last << std::endl;
			tempI+=tempStep;
			fileIO.open(name,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
		}

		if(i%measureint==0 && i!=startInt)
		{
			//this is here because one needs both potentials to gauge the system accurately
			//this really isn't needed because the second part of the integration doesn't modify position
			//pairInteractions.build();
			double potential=pairInteractions.computePotential();
			potential+=moleculeInteractions.computePotential();
			std::fstream dataFile;
			std::string buf("potential_");
			buf+=argv[1];
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			dataFile << ti << '\t' << potential << std::endl;
			dataFile.close();
			kinetic.output(ti,argv[1]);
			//if you want to add a new measure when you do a measure step, add it here
		}
		//someMeasureInterval would be an integer, like every 10 steps rather than 10 tau
		//if(i%someMeasureInterval && i!=0)
		//{
		//	//you would put the measure here
		//}
		//you could potentially add more measures here
	}
	
	MPI_Finalize();
	return 0;
}

