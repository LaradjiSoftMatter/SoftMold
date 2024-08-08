/**
 * \brief This simply runs a system through the timesteps. It can utilize OpenMP.
 * It is currently setup to use the implicit solvent molecular dynamic model used by Laradji's Computational Soft Matter lab. 
 * The particular system being studied is a phospholipid system. 
 * Various additions, such as cytoskeletons and nanoparticles, have been studied as well.
 * Multicomponent lipid systems have been studied.
 */

//Comment these flags out if you don't want them:

//Enable error readout and halting
//#define ERRORS_ENABLED

//Enable warning readout and halting
//#define WARNINGS_ENABLED

//For anchor data
//#define ANCHOR_DATA

//Flag to recenter the mass in the system
//#define RECENTER_MASS

//only use the filled boxes, I'm going to change this to only use particle to thread decomposition later
//#define LOW_DENSITY

//If the local density is high, the threads conflict, this fixes that by skipping some cells when low density is active
//#define LOW_DENSITY_SKIP 7

//When to start computing diffusion, in tau
//#define DIFFUSION_START 100.0

//For molecular dynamics forces and potentials
#include "include/MD.h"

//For the molecular dynamics variables
#include "include/system.h"

//For quick single processor MD routines
#include "include/MDroutines.h"

//an object that handles memory accross MPI threads
#include "include/MPImem.h"

//extractData must be located after system
//It is a header dependent on the system in the include directory, but is not part of the program
#include "dataExtraction.h"
#include <ctime>


//frequency that the system is resized, 1=every time step, 2=every other time step, 3=every 3rd time step, etc...
#define resizeRate 1

int main(int argc, char* argv[])
{
	MPI_Init (&argc, &argv);
	int rank, numProc;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &numProc);
	
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name\n";
		MPI_Finalize();
		return 0;
	}
	
	char *name=argv[1];
	
	//the variables for the simulation, remember that for some reason constructor isn't explicit when nothing is in it
	Blob<double> System;
	position <double> *pos=NULL;
	threeVector<double> *acc=NULL;
	threeVector<double> *vel=NULL;
	threeVector<double> size;
	double dT;
	double gamma;
	int nParticles;
	double temperature;
	double cutoff;
	double *tBFC;
	int nTypes;
	int seed;
	double time;
	
	
	//for unrolling bonding and bending lists
	std::vector< threeVector<int> > bendList;//index is the key retrieved from bendConstIndex
	std::vector< twoVector<int> > bondList;//index is the key retrieved from bondConstIndex
	std::vector< twoVector<double> > bendConst;//same index as bendConstIndex
	std::vector< twoVector<double> > bondConst;//same index as bondConstIndex
	std::vector<int> bendConstIndex;//ordered index key pairs
	std::vector<int> bondConstIndex;//ordered index key pairs
	
	//these never change during simulation
	int endInt, startInt, storeInt, measureInt, tempStepInt;
	double tempStep;
	
	if(rank==0)
	{
		//load variables, then initialize them, Script requires some functions from Blob
		Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
		
		fileIO.read();
		
		fileIO.close();
		
		//using integer indexing, the 0.0000001 fixes an accuracy issue with the gnu c++ compiler.
		//Don't believe it affects anything else...
		endInt=int(System.readFinalTime()/System.readDeltaT()+0.0000001);//end
		startInt=int(System.readInitialTime()/System.readDeltaT()+0.0000001);//start
		
		storeInt=int(System.readStoreInterval()/System.readDeltaT()+0.0000001);//when to store
		measureInt=int(System.readMeasureInterval()/System.readDeltaT()+0.0000001);//when to measure
		
		tempStepInt=0;
		tempStep=0;
		
		if(System.readTempStepInterval()>0)
		{
			tempStep=System.readTempStepInterval()*(System.readFinalTemp()-System.readInitialTemp())/(System.readFinalTime()-System.readInitialTime());
			tempStepInt=(endInt-startInt)/int((System.readFinalTime()-System.readInitialTime())/(System.readTempStepInterval()));
		}
		
		//extract the bonds and bends
		molecule<double, fourVector<int> > *m=System.getMolecule();
		for(int i=0;i<System.readNMolecules();i++)
		{
			switch(m[i].readType())
			{
				case BOND:
				{
					twoVector<double> constants;
					constants.s[ABOND]=m[i].getConstants()[ABOND];
					constants.s[KBOND]=m[i].getConstants()[KBOND];
					bondConst.push_back(constants);
					
					bondConstIndex.push_back(bondList.size());
					
					for(int j=0;j<m[i].readNBond();j++)
					{
						twoVector<int> bond;
						bond.s[0]=m[i].getBonds()[j].s[0];
						bond.s[1]=m[i].getBonds()[j].s[1];
						bondList.push_back(bond);
					}
				}
				break;
				
				case BEND:
				{
					twoVector<double> constants;
					constants.s[ABEND]=m[i].getConstants()[ABEND];
					constants.s[KBEND]=m[i].getConstants()[KBEND];
					bendConst.push_back(constants);
					
					bendConstIndex.push_back(bendList.size());
					
					for(int j=0;j<m[i].readNBond();j++)
					{
						threeVector<int> bend;
						bend.s[0]=m[i].getBonds()[j].s[0];
						bend.s[1]=m[i].getBonds()[j].s[1];
						bend.s[2]=m[i].getBonds()[j].s[2];
						bendList.push_back(bend);
					}
				}
				break;
				
				case CHAIN:
				{
					twoVector<double> bondConstants, bendConstants;
					bondConstants.s[ABOND]=m[i].getConstants()[CHAINBOND+ABOND];
					bondConstants.s[KBOND]=m[i].getConstants()[CHAINBOND+KBOND];
					
					bendConstants.s[ABEND]=m[i].getConstants()[CHAINBEND+ABEND];
					bendConstants.s[KBEND]=m[i].getConstants()[CHAINBEND+KBEND];
					
					bondConst.push_back(bondConstants);
					bendConst.push_back(bendConstants);
					
					bondConstIndex.push_back(bondList.size());
					bendConstIndex.push_back(bendList.size());
					
					for(int j=0;j<m[i].readNBond();j++)
					{
						int start=m[i].getBonds()[j].s[START];
						int nChains=m[i].getBonds()[j].s[NCHAINS];
						int length=m[i].getBonds()[j].s[LENGTH];
						for(int k=start;k<start+nChains*length;k+=length)
						{
							for(int l=k;l<k+length-1;l++)
							{
								twoVector<int> bond;
								bond.s[0]=l;
								bond.s[1]=l+1;
								bondList.push_back(bond);
							}
							for(int l=k;l<k+length-2;l++)
							{
								threeVector<int> bend;
								bend.s[0]=l;
								bend.s[1]=l+1;
								bend.s[2]=l+2;
								bendList.push_back(bend);
							}
						}
					}
				}
				break;
				
				default:
					break;
			}
		}
		
		//reorder bonds and bends to make searching easier
		for(int i=0;i<bondConstIndex.size();i++)
		{
			int start=bondConstIndex[i];
			int end=bondList.size();
			if(i<bondConstIndex.size()-1)
				end=bondConstIndex[i+1];
			std::sort(bondList.begin()+start, bondList.begin()+end, sortByDim<twoVector<int>, 0>);
		}
		for(int i=0;i<bendConstIndex.size();i++)
		{
			int start=bendConstIndex[i];
			int end=bendList.size();
			if(i<bendConstIndex.size()-1)
				end=bendConstIndex[i+1];
			std::sort(bendList.begin()+start, bendList.begin()+end, sortByDim<threeVector<int>, 0>);
		}
	}
	mpiMemory<double> memoryHandler;
	
	memoryHandler.initComm(System.getPositions(), System.getVelocities(), System.getAccelerations(),
			      System.readNParticles(), &bondList[0], bondList.size(), &bendList[0], bendList.size(),
			      &bondConst[0], &bondConstIndex[0], bondConst.size(),
			      &bendConst[0], &bendConstIndex[0], bendConst.size(),
			      System.readSize());
	
	//For testing
	MPI_Finalize();
	return 0;
	
	
	//memory handler
	//mpiMemory memoryHandler(System.getPositions(), System.getVelocities(), System.getAccelerations(),
	//		System.readNParticles(), System.getMolecule(), System.readNMolecules(),
	//		System.getTwoBodyFconst(), System.getTwoBodyUconst(), System.readNTypes(), 
	//		System.readSize(), System.readCutoff(), nThreads)
	if(rank==0 && false)//this will change when all threads are active
	{
		dataExtraction<double, Blob <double> > dataCollection(&System,name);
		
		pos=System.getPositions();
		vel=System.getVelocities();
		acc=System.getAccelerations();
		nParticles=System.readNParticles();
		dT=System.readDeltaT();
		gamma=System.readGamma();
		temperature=System.readInitialTemp();
		cutoff=System.readCutoff();
		tBFC=System.getTwoBodyFconst();
		nTypes=System.readNTypes();
		size=System.readSize();
		seed=System.readSeed()+rank*139;
		time=System.readInitialTime();
		
		//for random size fluctuations
		MTRand randNum(seed);
		
		//Other initializations
		std::string framesFilename("frames_");
		framesFilename+=name;
		framesFilename+=".xyz";
		
		for(int k=0;k<nParticles;k++)
			acc[k]=0;
		
		langevinThermostat(vel, acc, nParticles, temperature, dT, gamma, &randNum);
		
		//Force is computed here regardless of previous state.
		computeForce<double,Force<double> > (pos, vel, acc, nParticles, cutoff, tBFC, nTypes
						     #ifdef SIZE_BOUND
						     ,size
						     #endif
		);
		
		//Cuda should work for this
		for(int k=0;k<bondConstIndex.size();k++)
		{
			twoVector<double> bC=bondConst[k];
			int end=(k==bondConstIndex.size()-1)?bondList.size():bondConstIndex[k+1];
			for(int j=bondConstIndex[k];j<end;j++)
			{
				twoVector<int> bonded=bondList[j];
				threeVector<double> d;
				d.x=pos[bonded.s[0]].x-pos[bonded.s[1]].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=pos[bonded.s[0]].y-pos[bonded.s[1]].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=pos[bonded.s[0]].z-pos[bonded.s[1]].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				harmonicF(d, acc[bonded.s[0]], acc[bonded.s[1]], bC.s);
			}
		}
		
		//Cuda should work for this
		for(int k=0;k<bendConstIndex.size();k++)
		{
			twoVector<double> bC=bendConst[k];
			int end=(k==bendConstIndex.size()-1)?bendList.size():bendConstIndex[k+1];
			for(int j=bendConstIndex[k];j<end;j++)
			{
				threeVector<int> bended=bendList[j];
				threeVector<double> da,db;
				da.x=pos[bended.s[0]].x-pos[bended.s[1]].x;
				if(da.x>size.x/2.0) da.x-=size.x;
				if(da.x<-size.x/2.0) da.x+=size.x;
				da.y=pos[bended.s[0]].y-pos[bended.s[1]].y;
				if(da.y>size.y/2.0) da.y-=size.y;
				if(da.y<-size.y/2.0) da.y+=size.y;
				da.z=pos[bended.s[0]].z-pos[bended.s[1]].z;
				if(da.z>size.z/2.0) da.z-=size.z;
				if(da.z<-size.z/2.0) da.z+=size.z;
				
				db.x=pos[bended.s[1]].x-pos[bended.s[2]].x;
				if(db.x>size.x/2.0) db.x-=size.x;
				if(db.x<-size.x/2.0) db.x+=size.x;
				db.y=pos[bended.s[1]].y-pos[bended.s[2]].y;
				if(db.y>size.y/2.0) db.y-=size.y;
				if(db.y<-size.y/2.0) db.y+=size.y;
				db.z=pos[bended.s[1]].z-pos[bended.s[2]].z;
				if(db.z>size.z/2.0) db.z-=size.z;
				if(db.z<-size.z/2.0) db.z+=size.z;
				
				bendF(da, db, acc[bended.s[0]], acc[bended.s[1]], acc[bended.s[2]],bC.s);
			}
		}
		
		dataCollection.initialize();
		
		//this corrects an issue where an extra data point is added when the system is restarted
		if(time==0)
		{
			dataCollection.compute();
			xyzFormat<double> xyzFile(pos, nParticles);
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
		}
		else
		{
			//Surprise! This is done because a previously run configuration doesn't do this upon exit
			secondIntegration(pos, vel, acc, nParticles, dT);
		}
		
		std::cout << "starting main loop: \n";
		
		time_t current=std::time(NULL);
		
		bool exitFlag=false;//for premature exits!
		
		//the molecular dynamics loop, the "running" of the system
		for(int i=startInt;i<=endInt && !exitFlag;i++)
		{
			pos=System.getPositions();
			vel=System.getVelocities();
			acc=System.getAccelerations();
			nParticles=System.readNParticles();
			dT=System.readDeltaT();
			gamma=System.readGamma();
			temperature=System.readInitialTemp();
			cutoff=System.readCutoff();
			tBFC=System.getTwoBodyFconst();
			nTypes=System.readNTypes();
			size=System.readSize();
			
			time=time+static_cast<double>(i)*dT;
			firstIntegration(pos, vel, acc, nParticles, dT);
			
			for(int k=0;k<nParticles;k++)
			{
				acc[k].x=0;
				acc[k].y=0;
				acc[k].z=0;
			}
			
			//The system is stored here because force is updated here, but velocities are updated next.
			//It causes a problem when it reenters the loop from a previously run configuration.
			if(tempStepInt>0)
				if(i%tempStepInt==0 && tempStepInt!=0 && i<endInt)
					temperature+=tempStep;
			
			if(i%storeInt==0 && i!=startInt)
			{
				//fileIO.open(name,std::ios::out);
				//fileIO.write();
				//fileIO.close();
				
				xyzFormat<double> xyzFile(pos, nParticles);
				xyzFile.open(framesFilename,std::ios::out | std::ios::app);
				xyzFile.store();
				xyzFile.close();
			}
			
			//Cuda works for these
			langevinThermostat(vel, acc, nParticles, temperature, dT, gamma, &randNum);
			
			//Calculate Force
			computeForce<double,Force<double> > (pos, vel, acc, nParticles, cutoff, tBFC, nTypes
							     #ifdef SIZE_BOUND
							     ,size
							     #endif
			);
			
			//Cuda should work for this
			for(int k=0;k<bondConstIndex.size();k++)
			{
				twoVector<double> bC=bondConst[k];
				int end=(k==bondConstIndex.size()-1)?bondList.size():bondConstIndex[k+1];
				for(int j=bondConstIndex[k];j<end;j++)
				{
					twoVector<int> bonded=bondList[j];
					threeVector<double> d;
					d.x=pos[bonded.s[0]].x-pos[bonded.s[1]].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=pos[bonded.s[0]].y-pos[bonded.s[1]].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=pos[bonded.s[0]].z-pos[bonded.s[1]].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					harmonicF(d, acc[bonded.s[0]], acc[bonded.s[1]], bC.s);
				}
			}
			
			//Cuda should work for this
			for(int k=0;k<bendConstIndex.size();k++)
			{
				twoVector<double> bC=bendConst[k];
				int end=(k==bendConstIndex.size()-1)?bendList.size():bendConstIndex[k+1];
				for(int j=bendConstIndex[k];j<end;j++)
				{
					threeVector<int> bended=bendList[j];
					threeVector<double> da,db;
					da.x=pos[bended.s[0]].x-pos[bended.s[1]].x;
					if(da.x>size.x/2.0) da.x-=size.x;
					if(da.x<-size.x/2.0) da.x+=size.x;
					da.y=pos[bended.s[0]].y-pos[bended.s[1]].y;
					if(da.y>size.y/2.0) da.y-=size.y;
					if(da.y<-size.y/2.0) da.y+=size.y;
					da.z=pos[bended.s[0]].z-pos[bended.s[1]].z;
					if(da.z>size.z/2.0) da.z-=size.z;
					if(da.z<-size.z/2.0) da.z+=size.z;
					
					db.x=pos[bended.s[1]].x-pos[bended.s[2]].x;
					if(db.x>size.x/2.0) db.x-=size.x;
					if(db.x<-size.x/2.0) db.x+=size.x;
					db.y=pos[bended.s[1]].y-pos[bended.s[2]].y;
					if(db.y>size.y/2.0) db.y-=size.y;
					if(db.y<-size.y/2.0) db.y+=size.y;
					db.z=pos[bended.s[1]].z-pos[bended.s[2]].z;
					if(db.z>size.z/2.0) db.z-=size.z;
					if(db.z<-size.z/2.0) db.z+=size.z;
					
					bendF(da, db, acc[bended.s[0]], acc[bended.s[1]], acc[bended.s[2]],bC.s);
				}
			}
			
			secondIntegration(pos, vel, acc, nParticles, dT);
			
			//Measurements are output here
			if(i%measureInt==0 && i!=startInt)
			{
				//double last=current;
				time_t last=current;
				//current=omp_get_wtime();
				current=std::time(NULL);
				//time since last storage step, good for benchmarking
				std::cout << time << '\t' << current-last << std::endl;
				
				//distributedMemory.gather();
				
				//Data calculations
				dataCollection.compute();
				
				//Change execution conditions due to calculations
				#ifdef ANCHOR_DATA
				if(dataCollection.readNAnchors()>0)
					if(dataCollection.readNBrokenAnchors()/dataCollection.readNAnchors()>0.5)
						exitFlag=true;
				#endif
				current=std::time(NULL);
			}
			//end=omp_get_wtime();
			//std::cout << "measure: " << end-start << '\n';
			//start=omp_get_wtime();		
			
			//someMeasureInterval would be an integer, like every 10 steps rather than 10 tau
			//if(i%someMeasureInterval && i!=0)
			//{
			//	//you would put the measure here
			//}
			//you could potentially add more measures here
			//end=omp_get_wtime();
			//std::cout << "Resize: " << end-start << '\n';
			//std::cin.get();
		}
	}
	MPI_Finalize();
	return 0;
}
