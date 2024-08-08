/**
 * \brief This simply runs a system through the timesteps. It can utilize OpenMP.
 * It is currently setup to use the implicit solvent molecular dynamic model used by Laradji's Computational Soft Matter lab. 
 * The particular system being studied is a phospholipid system. 
 * Various additions, such as cytoskeletons and nanoparticles, have been studied as well.
 * Multicomponent lipid systems have been studied.
 */

//no solvent within frames
#define noSolventFrames

//For molecular dynamics forces and potentials
#include "include/MD.h"

//For the molecular dynamics variables
#include "include/system.h"

//extractData must be located after system
//It is a header dependent on the system in the include directory, but is not part of the program
#include "dataExtraction.h"
#include <ctime>

#define resizeRate 1

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name\n";
		return 0;
	}
	
	char *name=argv[1];
	
	//the variables for the simulation, remember that for some reason constructor isn't explicit when nothing is in it
	Blob<double> System;
	
	//load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	//initialize the variables, objects, algorithms, data collection, etc...
	Verlet<double> integrate(System.getPositions(), System.getAccelerations(), System.getVelocities(), System.readNParticles(),
				 System.readSize(), System.readDeltaT(), System.readPeriodic());
	Langevin<double> thermostat(System.getAccelerations(), System.getVelocities(), System.readNParticles(), System.readGamma(),
				    System.readDeltaT(), System.readSeed());
	
	//for data
	dataExtraction<double, Blob <double> > dataCollection(&System,name);
	
	//for xyz frames
	position <double> *p=System.getPositions();
	int xyzNParticles=System.readNParticles();
	xyzFormat<double> xyzFile(p, xyzNParticles);
	
	//for random size fluctuations
	MTRand randNum(System.readSeed());
	
	//these two do pretty much the same thing, but CellOpt is much much much faster
	CellOpt<double, Potential<double>, Force <double> > pairInteractions(System.getPositions(), System.getAccelerations(), 
		System.getTwoBodyFconst(), System.getTwoBodyUconst(), System.readNParticles(), System.readNTypes(), System.readSize(),
		System.readPeriodic(), System.readCutoff());
	//Cell<double> neighbors(System.p, System.nParticles, System.cutoff, System.size);
	
	//Map and count the volumes
	int *excludeType=new int[System.readNTypes()];
	for(int i=1;i<System.readNTypes();i++)
		excludeType[i-1]=i;
	
	VolumeExtraction<double> volumize(System.getPositions(), System.readNParticles(),\
		System.readSize(), System.readCutoff(), excludeType, System.readNTypes()-1, System.readSeed());
	delete excludeType;
	
	//Other initializations
	std::string framesFilename("frames_");
	framesFilename+=name;
	framesFilename+=".xyz";
	
	int nSolvent=0;
	if(System.readRemoveSolvent()>0)
		for(int i=0;i<System.readNParticles();i++)
			#ifndef SOLVENT_FLAG
				if(p[i].type==SOLVENT)
					nSolvent++;
			#else
				if(p[i].type==SOLVENT_FLAG)
					nSolvent++;
			#endif
	
	//This whole section builds counter interactions for BOND structures
	/*
	int nBondMol=0;
	for(int i=0;i<System.readNMolecules();i++)
		if(System.getMolecule()[i].readType()==BOND)
			nBondMol++;
	
	CellOpt<double, CounterPotential<double>, CounterForce<double> > *bondCounterInteractions=NULL;
	int **bondList=NULL;
	if(nBondMol!=0)
	{
		bondList=new int*[nBondMol];//each list
		bondCounterInteractions= new CellOpt<double, CounterPotential<double>, CounterForce<double> > [nBondMol];
	}
	
	for(int i=0,listIndex=0;i<System.readNMolecules();i++)
	{
		if(System.getMolecule()[i].readType()==BOND)
		{
			int nBondList=(System.getMolecule()[i].readNBond())*2;
			bondList[listIndex]=new int[nBondList];
			for(int j=0,bondIndex=0;j<System.getMolecule()[i].readNBond();j++)
			{
				bondList[listIndex][bondIndex++]=System.getMolecule()[i].getBonds()[j].s[0];
				bondList[listIndex][bondIndex++]=System.getMolecule()[i].getBonds()[j].s[1];
			}
			std::sort(&bondList[listIndex][0],&bondList[listIndex][nBondList]);
			int *end=std::unique(&bondList[listIndex][0],&bondList[listIndex][nBondList]);
			//Supposedly this works with pointer arithmetic?
			nBondList=(end-bondList[listI2.3607759e+08ndex]);
			bondCounterInteractions[listIndex].initialize(System.getPositions(), System.getAccelerations(), 
				System.getTwoBodyFconst(), System.getTwoBodyUconst(), System.readNParticles(), 
				System.readNTypes(), System.readSize(), System.readPeriodic(), System.readCutoff(),
				nBondList,bondList[listIndex]);
			listIndex++;
		}
	}
	*/
	//End of building bond counter interactions, 
	// the actual execution is below, but commented out
	
	for(int k=0;k<System.readNParticles();k++)
		System.getAccelerations()[k]=0;
	
	thermostat.compute(System.readInitialTemp());

	//this is how you call the pair interactions, build then compute...
	//Cuda works for these
	//Force is computed here regardless of previous state.
	pairInteractions.build();
	pairInteractions.computeForce();
	
	//Cuda doesn't work for this loop yet
	for(int k=0;k<System.readNMolecules();k++)
	{
		switch(System.getMolecule()[k].readType())
		{
			case BOND:
			{
				System.doBondForce(k);
				break;
			}
			case BEND:
			{
				System.doBendForce(k);
				break;
			}
			case CHAIN:
			{
				System.doChainForce(k);
				break;
			}
			default:
			{
				//does nothing
				break;
			}
		}
	}
	
	//counter-force
	//for(int k=0;k<nBondMol;k++)
	//{
	//	bondCounterInteractions[k].build();
	//	bondCounterInteractions[k].computeForce();
	//}
	
	//this corrects an issue where an extra data point is added when the system is restarted
	if(System.readInitialTime()==0)
	{
		dataCollection.initialize();
		
		xyzFile.open(framesFilename,std::ios::out | std::ios::app);
		xyzFile.store();
		xyzFile.close();
	}
	else
	{
		//Surprise! This is done because a previously run configuration doesn't do this upon exit
		integrate.second();
	}
	
	//using integer indexing, the 0.0000001 fixes an accuracy issue with the gnu c++ compiler.
	//Don't believe it affects anything else...
	int endInt=int(System.readFinalTime()/System.readDeltaT()+0.0000001);//end
	int startInt=int(System.readInitialTime()/System.readDeltaT()+0.0000001);//start
	
	int storeint=int(System.readStoreInterval()/System.readDeltaT()+0.0000001);//when to store
	int measureint=int(System.readMeasureInterval()/System.readDeltaT()+0.0000001);//when to measure
	
	double tempStep=(System.readStoreInterval()*(System.readInitialTemp()-System.readFinalTemp()))/(System.readFinalTime()-System.readInitialTime());
	
	
	std::cout << "starting main loop: \n";
	
	time_t current=time(NULL);
	
	//the molecular dynamics loop, the "running" of the system
	for(int i=startInt;i<=endInt;i++)
	{
		threeVector<double> *acc=System.getAccelerations();
		System.setInitialTime((double)i*System.readDeltaT());
		
		integrate.first();
		
		for(int k=0;k<System.readNParticles();k++)
		{
			//This version is really slow!
			//System.getAccelerations()[k]=0;
			
			//Direct access is much faster
			acc[k].x=0;
			acc[k].y=0;
			acc[k].z=0;
		}
		
		//The system is stored here because force is updated here, but velocities are updated next.
		//It causes a problem when it reenters the loop from a previously run configuration.
		if(i%storeint==0 && i!=startInt)
		{
			System.setInitialTemp(System.readInitialTemp()+tempStep);
			fileIO.open(name,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
		}
		
		//Cuda works for these
		thermostat.compute(System.readInitialTemp());
		
		//Build linked lists
		pairInteractions.build();
		
		pairInteractions.computeForce();
		
		//Cuda doesn't work for this loop yet
		for(int k=0;k<System.readNMolecules();k++)
		{
			switch(System.getMolecule()[k].readType())
			{
				case BOND:
				{
					System.doBondForce(k);
					break;
				}
				case BEND:
				{
					System.doBendForce(k);
					break;
				}
				case CHAIN:
				{
					System.doChainForce(k);
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
		}
		
		//counter-force
		//for(int k=0;k<nBondMol;k++)
		//{
		//	bondCounterInteractions[k].build();
		//	bondCounterInteractions[k].computeForce();
		//}
		
		integrate.second();
		
		//this just needs to be placed before measure, it just moves inner solvent particles to outer
		if(i<System.readRemoveSolvent()*nSolvent)
		{
			//These 2 lines are very slow, there might be a better way
			volumize.build();
			volumize.moveToOuter(volumize.grabInner());
		}
		
		//Measurements are done here
		if(i%measureint==0 && i!=startInt)
		{
			//double last=current;
			time_t last=current;
			//current=omp_get_wtime();
			current=time(NULL);
			//time since last storage step, good for benchmarking
			std::cout << System.readInitialTime() << '\t' << current-last << std::endl;
			
			dataCollection.compute();
			
			current=time(NULL);
		}
		
		//short section to resize system, note that it only works when deltaLXY is something other than 0, it flags execution.
		//This needs to be put in it's own object.
		if(i%resizeRate==0 && i!=0 && System.readDeltaLXY()!=0)
		{
			threeVector<double> size=System.readSize();
			threeVector<double> oldSize=System.readSize();
			
			threeVector<double> fluctuation;
			fluctuation.x=System.readDeltaLXY()*(2.0*randNum.rand53()-1.0);
			fluctuation.y=System.readDeltaLXY()*(2.0*randNum.rand53()-1.0);
			//constrain z to maintain volume, but vary projected area
			//Note to anyone who sees this later, this method is very anistropic.
			// If you have any issue, try lowering deltaLXY or exclude any new types
			// from the interaction.
			fluctuation.z=(size.x*size.y)/((size.x+fluctuation.x)*(size.y+fluctuation.y));
			//System.readDeltaLXY()*(2.0*randNum.rand53()-1.0);
			size.x+=fluctuation.x;
			size.y+=fluctuation.y;
			size.z*=fluctuation.z;
			
			threeVector<double> aSize=size;
			aSize.x/=oldSize.x;
			aSize.y/=oldSize.y;
			aSize.z/=oldSize.z;
			double oldVolume=oldSize.x*oldSize.y*oldSize.z;
			double newVolume=size.x*size.y*size.z;
			
			double dPotential=pairInteractions.computeDPotential(aSize);
			for(int k=0;k<System.readNMolecules();k++)
			{
				//pick a structure by type
				switch(System.getMolecule()[k].readType())
				{
					case BOND:
					{
						dPotential+=System.doBondDPotential(k,aSize);
						break;
					}
					case BEND:
					{
						dPotential+=System.doBendDPotential(k,aSize);
						break;
					}
					case CHAIN:
					{
						dPotential+=System.doChainDPotential(k,aSize);
						break;
					}
					default:
					{
						//does nothing
						break;
					}
				}
			}
			
			//The next line was the old way of doing this, the one below it is the new way.
			//Combined with the changes in cellOpt and systemMD, this is about 15% faster than the old one.
			//double D=(oldPotential-newPotential)/System.readInitialTemp();//this is a really small number sometimes
			double D=(-dPotential)/System.readInitialTemp();//this is a really small number sometimes
			
			//If volume isn't changed, then this doesn't contribute.
			//D-=(double(System.readNParticles())*log(oldVolume/newVolume));
			D=exp(-D);
			double randNumber=randNum.rand53();
			
			//accept change
			if(D<=randNum.rand53())
			{
				for(int k=0;k<System.readNParticles();k++)
				{
					//This is kind of a template for excluding a type if needed
					//Don't forget to do the same thing to the above section as well
					//if(Sys.getP()[k].type!=excludedType)
					//{
						p[k].x*=(size.x/oldSize.x);
						p[k].y*=(size.y/oldSize.y);
						p[k].z*=(size.z/oldSize.z);
						//}
			}
			
			pairInteractions.resize(size);
			integrate.resize(size);
			System.setSize(size);
		}
		}
		
		//someMeasureInterval would be an integer, like every 10 steps rather than 10 tau
		//if(i%someMeasureInterval && i!=0)
		//{
		//	//you would put the measure here
		//}
		//you could potentially add more measures here
	}
	
	//Clear any memory allocated in main:
	/*
	if(nBondMol!=0)
	{
		delete[] bondCounterInteractions;
		for(int i=0;i<nBondMol;i++)
		{
			delete bondList[i];
		}
		delete bondList;
	}
	*/
	return 0;
}

