/**
 * \brief This simply runs a system through the timesteps. It can utilize OpenMP.
 * It is currently setup to use the implicit solvent molecular dynamic model used by Laradji's Computational Soft Matter lab. 
 * The particular system being studied is a phospholipid system. 
 * Various additions, such as cytoskeletons and nanoparticles, have been studied as well.
 * Multicomponent lipid systems have been studied.
 */


//Macros that can be overloaded in the library (comment them out if you don't want them):

//turning this flag off saves about 5% of runtime, but could also result in failure
//#define CELL_SIZE_FAILURE

//related to the above, if you turn on CELL_SIZE_FAILURE and receive a failure, make this bigger
//#define MAX_CELL_SIZE 8192

//sort optimization, testing, potentially better
//#define SORT_OPTIMIZATION

//omp molecules, testing, bullshit so far
//#define OMP_MOLECULES

//Enable error readout and halting
//#define ERRORS_ENABLED

//Enable warning readout and halting
//#define WARNINGS_ENABLED

//For anchor data
//#define ANCHOR_DATA

//only use the filled boxes
#define LOW_DENSITY

//If the local density is high, the threads conflict, this fixes that by skipping some cells when low density is active
#define LOW_DENSITY_SKIP 7

//a debugging flag
//#define FORCE_COUNT


//Include files from the library:

//For molecular dynamics forces and potentials
#include "include/MD.h"

//For the molecular dynamics variables
#include "include/system.h"

//For data extraction that is already enabled
#include "dataExtraction.h"
#include <ctime>




//Macros that are used in this file:

//Flag to recenter the mass in the system
//#define RECENTER_MASS

//When to start computing diffusion, in tau
#define DIFFUSION_START 100.0

//frequency that the system is resized, 1=every time step, 2=every other time step, 3=every 3rd time step, etc...
#define resizeRate 4




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
	
	threeVector<double> *acc=System.getAccelerations();
	
	//for diffusion
	position <double> *aP=new position<double>[System.readNParticles()];
	
	//for xyz frames
	position <double> *p=System.getPositions();
	int xyzNParticles=System.readNParticles();
	xyzFormat<double> xyzFile(p, xyzNParticles);
	
	//for diffusion
	for(int i=0;i<xyzNParticles;i++)
		aP[i]=p[i];
	
	double resizeHistInterval=0.00001;
	std::vector<double> resizeHist;
	if(System.readDeltaLXY()!=0)
	{
		//plus 1 for the end inclusion
		//resizeHist.alloc(static_cast<int>(System.readDeltaLXY()/resizeHistInterval)+1);
		int nIntervals=static_cast<int>(2.0*System.readDeltaLXY()/resizeHistInterval)+1;
		for(int i=0;i<nIntervals;i++)
			resizeHist.push_back(0);
	}


	//initialize the variables, objects, algorithms, data collection, etc...
	Verlet<double> integrate(System.getPositions(), System.getAccelerations(), System.getVelocities(), System.readNParticles(),
				 System.readSize(), System.readDeltaT(), System.readPeriodic(), aP);
	
	Langevin<double> thermostat;
	
	if(System.readGamma()>0)
	{
		thermostat.initialize(System.getAccelerations(), System.getVelocities(), System.readNParticles(), System.readGamma(),
				      System.readDeltaT(), System.readSeed());
	}
	else if(System.getGammaType()!=NULL)
	{
		thermostat.initialize(System.getAccelerations(), System.getVelocities(), System.getPositions(), System.readNParticles(),
				      System.getGammaType(), System.readDeltaT(), System.readSeed(), System.readNTypes());
	}
	else
	{
		std::cout << "Error(main): No gamma available!\n";
		return 0;
	}
	
	//for data
	dataExtraction<double, Blob <double> > dataCollection(&System,name,aP);
	
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
	
	//VolumeExtraction<double> volumize(System.getPositions(), System.readNParticles(),
	//	System.readSize(), System.readCutoff(), excludeType, System.readNTypes()-1, System.readSeed());
	delete excludeType;
	
	//Other initializations
	std::string framesFilename("frames_");
	framesFilename+=name;
	framesFilename+=".xyz";
	
	int nSolvent=0;
	if(System.readRemoveSolvent()>0)
	{
		#ifdef SOLVENT_FLAG
		for(int i=0;i<System.readNParticles();i++)
			if(p[i].type==SOLVENT_FLAG)
				nSolvent++;
		#endif
		#ifdef SOLVENT
		for(int i=0;i<System.readNParticles();i++)
			if(p[i].type==SOLVENT)
				nSolvent++;
		#endif
	}
	
	for(int k=0;k<System.readNParticles();k++)
		System.getAccelerations()[k]=0;
	
	//this is how you call the pair interactions, build then compute...
	//Cuda works for these
	//Force is computed here regardless of previous state.
	
	pairInteractions.build();
	pairInteractions.computeForce();
	
	//#ifdef FORCE_COUNT
	/*
	for(int i=0;i<System.readNParticles();i++)
	{
		threeVector<double> *a=System.getAccelerations();
		std::cout << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\t' << a[i].x << '\t' << a[i].y << '\t' << a[i].z << std::endl;
	}
	throw 0;
	*/
	//#endif
	
	thermostat.compute(System.readInitialTemp());
	
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
			/*
			case BEAD:
			{
				for(int j=0;j<System.getMolecule()[k].readNBonds();j++)
				{
					int beadIndex=System.getMolecule()[k].getBonds()[j].s[0];
					int beadType=System.getMolecule()[k].getBonds()[j].s[1];
					double r0=System.getMolecule()[k].getConstants()[j];
					int nTypes=System.readNTypes();
					double *constants=System.getConstants();
					//constants[3] is usually the cutoff of the system
					double cutoffSqr=constants[3]*constants[3];//System.readCutoff()*System.readCutoff();
					
					for(int l=0;i<System.readNParticles();l++)
					{
						if(l!=beadIndex)
						{
							threeVector<double> d;
							d.x=p[l].x-p[beadIndex].x-r0;
							d.y=p[l].y-p[beadIndex].y-r0;
							d.z=p[l].z-p[beadIndex].z-r0;
							
							int cIndex=p[l].type*nTypes+beadType;
							
							nonBondedF(d, acc[l], acc[beadIndex] , &constants[cindex] , cutoffSqr);
						}
					}
				}
			}
			*/
			default:
			{
				//does nothing
				break;
			}
		}
	}
	
	dataCollection.initialize();
	
	//this corrects an issue where an extra data point is added when the system is restarted
	if(System.readInitialTime()==0)
	{
		dataCollection.compute();
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
	
	
	int tempStepInt=0;
	double tempStep=0;
	if(System.readTempStepInterval()>0)
	{
		tempStep=System.readTempStepInterval()*(System.readFinalTemp()-System.readInitialTemp())/(System.readFinalTime()-System.readInitialTime());
		tempStepInt=(endInt-startInt)/int((System.readFinalTime()-System.readInitialTime())/(System.readTempStepInterval()));
	}
	
	double NNUmin=-0.5;
	
	for(int k=4;k<System.readNTypes();k++)
	{
		for(int l=k+1;l<System.readNTypes();l++)
		{
			int offsetIndex=nTWOBODYUCONST*(k+l*System.readNTypes());
			
			//F constants, force constants
			//System.getTwoBodyFconst()[0+offsetIndex]=1.0;//C8
			System.getTwoBodyFconst()[1+offsetIndex]=2.0*(200.0+NNUmin);//C6
			//System.getTwoBodyFconst()[2+offsetIndex]=0;//part of index trick
			//System.getTwoBodyFconst()[3+offsetIndex]=2.0;//C7
			System.getTwoBodyFconst()[4+offsetIndex]=6.0*NNUmin;//C3
			System.getTwoBodyFconst()[5+offsetIndex]=6.0*NNUmin;//C2
			
			//U constants, potential constants
			//System.getTwoBodyUconst()[0+offsetIndex]=1.0;//C8
			System.getTwoBodyUconst()[1+offsetIndex]=200-NNUmin;//C4
			System.getTwoBodyUconst()[2+offsetIndex]=NNUmin;//C5,no index trick
			//System.getTwoBodyUconst()[3+offsetIndex]=2.0;//C7
			System.getTwoBodyUconst()[4+offsetIndex]=3.0*NNUmin;//C1
			System.getTwoBodyUconst()[5+offsetIndex]=2.0*NNUmin;//C0
			
			offsetIndex=nTWOBODYUCONST*(l+k*System.readNTypes());
			
			//F constants, force constants
			//System.getTwoBodyFconst()[0+offsetIndex]=1.0;//C8
			System.getTwoBodyFconst()[1+offsetIndex]=2.0*(200.0+NNUmin);//C6
			//System.getTwoBodyFconst()[2+offsetIndex]=0;//part of index trick
			//System.getTwoBodyFconst()[3+offsetIndex]=2.0;//C7
			System.getTwoBodyFconst()[4+offsetIndex]=6.0*NNUmin;//C3
			System.getTwoBodyFconst()[5+offsetIndex]=6.0*NNUmin;//C2
			
			//U constants, potential constants
			//System.getTwoBodyUconst()[0+offsetIndex]=1.0;//C8
			System.getTwoBodyUconst()[1+offsetIndex]=200-NNUmin;//C4
			System.getTwoBodyUconst()[2+offsetIndex]=NNUmin;//C5,no index trick
			//System.getTwoBodyUconst()[3+offsetIndex]=2.0;//C7
			System.getTwoBodyUconst()[4+offsetIndex]=3.0*NNUmin;//C1
			System.getTwoBodyUconst()[5+offsetIndex]=2.0*NNUmin;//C0
		}
	}
	
	//std::cout << tempStep << '\n' << tempStepInt << '\n';
	//std::cout << System.readTempStepInterval() << '\n'; 
	std::cout << "starting main loop: \n";
	
	time_t current=time(NULL);
	
	bool exitFlag=false;//for premature exits!
	
	//the molecular dynamics loop, the "running" of the system
	for(int i=startInt;i<=endInt && !exitFlag;i++)
	{
//double start=omp_get_wtime();
		System.setInitialTime((double)i*System.readDeltaT());
		
//double end=omp_get_wtime();
//std::cout << "Set: " << end-start << '\n';
//start=omp_get_wtime();
		integrate.first();
//end=omp_get_wtime();
//std::cout << "Integrate: " << end-start << '\n';
//start=omp_get_wtime();		
		for(int k=0;k<System.readNParticles();k++)
		{
			//This version is really slow!
			//System.getAccelerations()[k]=0;
			
			//Direct access is much faster
			acc[k].x=0;
			acc[k].y=0;
			acc[k].z=0;
		}
//end=omp_get_wtime();
//std::cout << "reset A: " << end-start << '\n';
//start=omp_get_wtime();		
		//The system is stored here because force is updated here, but velocities are updated next.
		//It causes a problem when it reenters the loop from a previously run configuration.
		if(System.readTempStepInterval()>0)
			if(i%tempStepInt==0 && tempStepInt!=0 && i<endInt)
				System.setInitialTemp(System.readInitialTemp()+tempStep);
		
		if(i%storeint==0 && i!=startInt)
		{
			fileIO.open(name,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
			
			//Modify the nano-nano potential
			
			std::string buf("NNUmin_");
			buf+=name;
			buf+=".dat";
			
			std::fstream NNUminFile;
			NNUminFile.open(buf.c_str(),std::ios::out | std::ios::app);
			NNUminFile << System.readInitialTime() << '\t' << NNUmin << '\n';
			NNUminFile.close();
			
			NNUmin+=0.001;
			
			for(int k=4;k<System.readNTypes();k++)
			{
				for(int l=k+1;l<System.readNTypes();l++)
				{
					int offsetIndex=nTWOBODYUCONST*(k+l*System.readNTypes());
					
					//F constants, force constants
					//System.getTwoBodyFconst()[0+offsetIndex]=1.0;//C8
					System.getTwoBodyFconst()[1+offsetIndex]=2.0*(200.0+NNUmin);//C6
					//System.getTwoBodyFconst()[2+offsetIndex]=0;//part of index trick
					//System.getTwoBodyFconst()[3+offsetIndex]=2.0;//C7
					System.getTwoBodyFconst()[4+offsetIndex]=6.0*NNUmin;//C3
					System.getTwoBodyFconst()[5+offsetIndex]=6.0*NNUmin;//C2
					
					//U constants, potential constants
					//System.getTwoBodyUconst()[0+offsetIndex]=1.0;//C8
					System.getTwoBodyUconst()[1+offsetIndex]=200-NNUmin;//C4
					System.getTwoBodyUconst()[2+offsetIndex]=NNUmin;//C5,no index trick
					//System.getTwoBodyUconst()[3+offsetIndex]=2.0;//C7
					System.getTwoBodyUconst()[4+offsetIndex]=3.0*NNUmin;//C1
					System.getTwoBodyUconst()[5+offsetIndex]=2.0*NNUmin;//C0
					
					offsetIndex=nTWOBODYUCONST*(l+k*System.readNTypes());
					
					//F constants, force constants
					//System.getTwoBodyFconst()[0+offsetIndex]=1.0;//C8
					System.getTwoBodyFconst()[1+offsetIndex]=2.0*(200.0+NNUmin);//C6
					//System.getTwoBodyFconst()[2+offsetIndex]=0;//part of index trick
					//System.getTwoBodyFconst()[3+offsetIndex]=2.0;//C7
					System.getTwoBodyFconst()[4+offsetIndex]=6.0*NNUmin;//C3
					System.getTwoBodyFconst()[5+offsetIndex]=6.0*NNUmin;//C2
					
					//U constants, potential constants
					//System.getTwoBodyUconst()[0+offsetIndex]=1.0;//C8
					System.getTwoBodyUconst()[1+offsetIndex]=200-NNUmin;//C4
					System.getTwoBodyUconst()[2+offsetIndex]=NNUmin;//C5,no index trick
					//System.getTwoBodyUconst()[3+offsetIndex]=2.0;//C7
					System.getTwoBodyUconst()[4+offsetIndex]=3.0*NNUmin;//C1
					System.getTwoBodyUconst()[5+offsetIndex]=2.0*NNUmin;//C0
				}
			}
			
		}
//end=omp_get_wtime();
//std::cout << "Stupid Thing: " << end-start << '\n';
//start=omp_get_wtime();		
		//Cuda works for these
		thermostat.compute(System.readInitialTemp());
//end=omp_get_wtime();
//std::cout << "ThermoStat: " << end-start << '\n';
//start=omp_get_wtime();		
		//Build linked lists
		pairInteractions.build();
//end=omp_get_wtime();
//std::cout << "Build: " << end-start << '\n';
//start=omp_get_wtime();		
		pairInteractions.computeForce();
//end=omp_get_wtime();
//std::cout << "Compute: " << end-start << '\n';
//start=omp_get_wtime();		
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
//end=omp_get_wtime();
//std::cout << "Mol: " << end-start << '\n';
//start=omp_get_wtime();		
		integrate.second();
//end=omp_get_wtime();
//std::cout << "integrate 2: " << end-start << '\n';
//start=omp_get_wtime();		
		//this just needs to be placed before measure, it just moves inner solvent particles to outer
		if(i<System.readRemoveSolvent()*nSolvent)
		{
			//These 2 lines are very slow, there might be a better way
			//volumize.build();
			//volumize.moveToOuter(volumize.grabInner());
		}
//end=omp_get_wtime();
//std::cout << "if solvent: " << end-start << '\n';
//start=omp_get_wtime();
		
		//fast measurements are done in the main loop
		//if(i%20==0)
			//dataCollection.computeFast();
		
		//Measurements are output here
		if(i%measureint==0 && i!=startInt)
		{
			//double last=current;
			time_t last=current;
			//current=omp_get_wtime();
			current=time(NULL);
			//time since last storage step, good for benchmarking
			std::cout << System.readInitialTime() << '\t' << current-last << std::endl;
			
			
			//dataCollection.computeFast();
			
			//Data calculations that we are interested in starting
				if(System.readInitialTime()>DIFFUSION_START && System.readDeltaLXY()==0)
					dataCollection.startDiffusion();
			
			//Data calculations
				dataCollection.compute();
			
			//Change execution conditions due to calculations
			#ifdef ANCHOR_DATA
				if(dataCollection.readNAnchors()>0)
					if(dataCollection.readNBrokenAnchors()/dataCollection.readNAnchors()>0.5)
						exitFlag=true;
			#endif
			
				//Reason to remove net displacement:
				//	For systems where the scale changes, the membrane tends to wander upward.
				//	While there is no net momentum (net(velocity)~0), the act of resizing moves the bilayer.
				//	This recompensates for that by keeping the center of mass in the center of the system.
				//This doesn't remove net momentum! Don't forget the diffusion parameter aP doesn't rescale.
				
			#ifdef RECENTER_MASS
				if(System.readDeltaLXY()!=0)
				{
					//removing any net motion
					position<double> displaceSystem=com< position<double> >(System.getPositions(),System.readNParticles());
					displaceSystem.x-=System.readSize().x/2.0;
					displaceSystem.y-=System.readSize().y/2.0;
					displaceSystem.z-=System.readSize().z/2.0;
					
					for(int j=0;j<System.readNParticles();j++)
					{
						//move a particle by the displacement
						p[j].x-=displaceSystem.x;
						p[j].y-=displaceSystem.y;
						p[j].z-=displaceSystem.z;
						
						//check the boundaries
						p[j].x-=(p[j].x>=System.readSize().x)?System.readSize().x:0;
						p[j].y-=(p[j].y>=System.readSize().y)?System.readSize().y:0;
						p[j].z-=(p[j].z>=System.readSize().z)?System.readSize().z:0;
						
						p[j].x+=(p[j].x<0)?System.readSize().x:0;
						p[j].y+=(p[j].y<0)?System.readSize().y:0;
						p[j].z+=(p[j].z<0)?System.readSize().z:0;
					}
				}
			#endif
				
			current=time(NULL);
		}
//end=omp_get_wtime();
//std::cout << "measure: " << end-start << '\n';
//start=omp_get_wtime();		
		//short section to resize system, note that it only works when deltaLXY is something other than 0, it flags execution.
		//This needs to be put in it's own object.
		if(i%resizeRate==0 && i!=0 && System.readDeltaLXY()!=0)
		{
			threeVector<double> size=System.readSize();
			threeVector<double> oldSize=System.readSize();
			
			threeVector<double> fluctuation;
			fluctuation.x=System.readDeltaLXY()*(2.0*randNum.rand53()-1.0);
			fluctuation.y=System.readDeltaLXY()*(2.0*randNum.rand53()-1.0);
			//fluctuation.y=fluctuation.x;
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
			double D=(dPotential)/System.readInitialTemp();//this is a really small number sometimes
			
			//If volume isn't changed, then this doesn't contribute.
			//D-=(double(System.readNParticles())*log(oldVolume/newVolume));
			D=exp(-D);
			double randNumber=randNum.rand53();

//std::cout << dPotential << '\n';			
//std::cout << D << '\n' << fluctuation.x << '\t' << fluctuation.y << '\n';
//std::cin.get();

			//accept change
			if(D<=randNumber)
			{
				for(int k=0;k<System.readNParticles();k++)
				{
					//This is kind of a template for excluding a type if needed
					//Don't forget to do the same thing to the above section as well
					//if(Sys.getP()[k].type!=excludedType)
					//{
						p[k].x/=aSize.x;
						p[k].y/=aSize.y;
						p[k].z/=aSize.z;
					//}
				}
				resizeHist[(fluctuation.x+System.readDeltaLXY())/resizeHistInterval]+=0.5;
				resizeHist[(fluctuation.y+System.readDeltaLXY())/resizeHistInterval]+=0.5;
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
//end=omp_get_wtime();
//std::cout << "Resize: " << end-start << '\n';
//std::cin.get();
	}
	
	if(System.readDeltaLXY()!=0)
	{	
		std::fstream resizeHistFile;
		resizeHistFile.open("resizeHist.dat",std::ios::out);
		if(resizeHistFile.is_open())
		{
			for(int i=0;i<resizeHist.size();i++)
				resizeHistFile << (static_cast<double>(i)*resizeHistInterval)-System.readDeltaLXY()\
					 << '\t' << resizeHist[i] << '\n';
			resizeHistFile.close();
		}
	}
	
	if(aP!=NULL)
		delete aP;
	
	return 0;
}

