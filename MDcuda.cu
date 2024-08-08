/**
 * \brief This simply runs a system through the timesteps. It can utilize OpenMP.
 * It is currently setup to use the implicit solvent molecular dynamic model used by Laradji's Computational Soft Matter lab. 
 * The particular system being studied is a phospholipid system. 
 * Various additions, such as cytoskeletons and nanoparticles, have been studied as well.
 * Multicomponent lipid systems have been studied.
 * This is the GPU (CUDA) version of the MD.cpp program.
**/

//These macros won't work on the GPU directly (as of Jun 14, 2022)
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
//#define LOW_DENSITY

//If the local density is high, the threads conflict, this fixes that by skipping some cells when low density is active
//#define LOW_DENSITY_SKIP 7

//a debugging flag
//#define FORCE_COUNT


//Include files from the library:

//For molecular dynamics forces and potentials
#include "include/MD.h"

//For the molecular dynamics variables
#include "include/system.h"

//For data extraction that is already enabled
//#include "dataExtraction.h"

//For GPU
#include "include/cuda/mpdCuda.h"

#include <ctime>

//#define VERLET_HALF_CARRY


//Macros that are used in this file:

//Flag to recenter the mass in the system
//#define RECENTER_MASS

//When to start computing diffusion, in tau
#define DIFFUSION_START 100.0

//frequency that the system is resized, 1=every time step, 2=every other time step, 3=every 3rd time step, etc...
#define resizeRate 8




int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "usage: command name\n";
		return 0;
	}
	
	char *name=argv[1];
	std::string nameA(name);
	nameA+=".tmp";
	
	//the variables for the simulation, remember that for some reason constructor isn't explicit when nothing is in it
	Blob<float> System;
	
	//load variables, then initialize them, Script requires some functions from Blob
	Script<float, Blob <float> > fileIO(name,std::ios::in,&System);
	try {
		fileIO.read();
		fileIO.close();
	} catch(int e)
	{
		if(e==0)
		{
			fileIO.close();
			fileIO.open(nameA.c_str(),std::ios::in);
			fileIO.read();
			fileIO.close();
		}
		else
		{
			std::cerr << "Cannot open temp file! System load failure!" << std::endl;
			return -1;
		}
	}
	
	threeVector<float> *acc=System.getAccelerations();
	
	//for diffusion
	position <float> *aP=new position<float>[System.readNParticles()];
	
	//for xyz frames
	position <float> *p=System.getPositions();
	int nParticles=System.readNParticles();
	xyzFormat<float> xyzFile(p, nParticles);
	
	float trial=0, accepted=0;
	
	float resizeHistInterval=0.00001, maxInterval=System.readDeltaLXY();
	std::vector<float> resizeHist,rejectHist,dPotentialVal;
	if(System.readDeltaLXY()!=0)
	{
		//plus 1 for the end inclusion
		//resizeHist.alloc(static_cast<int>(System.readDeltaLXY()/resizeHistInterval)+1);
		int nIntervals=static_cast<int>(2.0*maxInterval/resizeHistInterval)+1;
		for(int i=0;i<nIntervals;i++)
			resizeHist.push_back(0);
		for(int i=0;i<nIntervals;i++)
			rejectHist.push_back(0);
	}
	
	
	//initialize the variables, objects, algorithms, data collection, etc...
	//for random size fluctuations
	std::mt19937 randGen(System.readSeed());
	threeVector<float> dL(System.readDeltaLXY(),System.readDeltaLXY(),0);
	mpd::barostat<float> bStat(dL, nParticles);
	mpd::randomAdaptor random(System.readSeed());
	mpd::interaction<float,3,2> bendList(nParticles);
	mpd::interaction<float,2,2> bondList(nParticles);
	//mpd::interaction<float,2,2> ballList(nParticles);
	std::vector<mpd::interaction<float,2,2>> ballLists;
	mpd::cell<float> cData(nParticles,System.readCutoff(),System.readSize(),System.readDeltaLXY());
	std::vector<float> mass(nParticles,1.0);
	
	//for data and resizing
	mpd::dataCollection<float> dataCollection(nParticles);
	
	//Add to the bond and bend lists
	for(int k=0;k<System.readNMolecules();k++)
	{
		auto mol=System.getMolecule()[k];
		switch(mol.readType())
		{
			case BOND:
			{
				auto b=mol.getBonds();
				auto c=mol.getConstants();
				for(int j=0;j<mol.readNBond();j++)
				{
					bondList.addInteraction(b[j].s[0],0,b[j].s,c);
					bondList.addInteraction(b[j].s[1],1,b[j].s,c);
				}
				break;
			}
			case BEND:
			{
				auto b=mol.getBonds();
				auto c=mol.getConstants();
				for(int j=0;j<mol.readNBond();j++)
				{
					bendList.addInteraction(b[j].s[0],0,b[j].s,c);
					bendList.addInteraction(b[j].s[1],1,b[j].s,c);
					bendList.addInteraction(b[j].s[2],2,b[j].s,c);
				}
				break;
			}
			case CHAIN:
			{
				auto b=mol.getBonds();
				auto c=mol.getConstants();
				for(int j=0;j<mol.readNBond();j++)
				{
					int start=b[j].s[START];
					int nChains=b[j].s[NCHAINS];
					int length=b[j].s[CHAINLENGTH];
					for(int i=start;i<start+length*nChains;i+=length)
					{
						std::vector<int> cha;
						for(int l=i;l<i+length;l++)
							cha.push_back(l);
						for(int l=0;l<length-2;l++)
						{
						bendList.addInteraction(cha[l],0,&(cha[l]),c+2);
						bendList.addInteraction(cha[l+1],1,&(cha[l]),c+2);
						bendList.addInteraction(cha[l+2],2,&(cha[l]),c+2);
						}
						for(int l=0;l<length-1;l++)
						{
						bondList.addInteraction(cha[l],0,&(cha[l]),c);
						bondList.addInteraction(cha[l+1],1,&(cha[l]),c);
						}
					}
				}
				
				break;
			}
			case BALL:
			{
				mpd::interaction<float,2,2> ballList(nParticles);
				auto b=mol.getBonds();
				auto c=mol.getConstants();
				mass[b[0].s[0]]=(4.0)*M_PI*c[0]*c[0];
				for(int j=0;j<mol.readNBond();j++)
				{
					//ballList.addInteraction(b[j].s[0],0,b[j].s,c);
					ballList.addInteraction(b[j].s[1],1,b[j].s,c);
				}
				ballLists.emplace_back(ballList);
				break;
			}
			default:
			{
				std::cerr << "Molecule type " << mol.readType() 
					<< " not yet implemented!" << std::endl;
				break;
			}
		}
	}
	
	mpd::state<float> state(System.getPositions(),System.getVelocities(),System.getAccelerations(),
				System.getTwoBodyFconst(),System.getTwoBodyUconst(),
				nParticles,System.readNTypes(),System.readDeltaT(),System.readGamma(),
				System.readInitialTemp(),System.readSize(),mass.data());
	
	//Other initializations
	std::string framesFilename("frames_");
	framesFilename+=name;
	framesFilename+=".xyz";
	
	std::string potentialFileName("potential_");
	potentialFileName+=name;
	potentialFileName+=".dat";
	
	std::string kineticFileName("kinetic_");
	kineticFileName+=name;
	kineticFileName+=".dat";
	
	std::string sizeFileName("size_");
	sizeFileName+=name;
	sizeFileName+=".dat";
	
	//Send to device
	bondList.toDevice();
	for(auto &ballList:ballLists)
		ballList.toDevice();
	bendList.toDevice();
	state.toDevice();
	
	//molecular dynamics forces
	mpd::zeroRange_device(state.deviceState().a,nParticles);
	mpd::bondForces_device(bondList.deviceInteraction(),state.deviceState());
	for(auto &ballList:ballLists)
		mpd::ballForces_device(ballList.deviceInteraction(),state.deviceState());
	mpd::bendForces_device(bendList.deviceInteraction(),state.deviceState());
	mpd::cellComputeForce_device(cData.deviceCell(), state.deviceState());
	mpd::applyMass_device(state.deviceState());
	
	
	//this corrects an issue where an extra data point is added when the system is restarted
	if(System.readInitialTime()==0)
	{
		mpd::zeroRange_device(dataCollection.deviceState().kineticEnergy,nParticles);
		mpd::kinetic_device(state.deviceState(),dataCollection.deviceState());
		float kinetic=
			mpd::reduceRange_device(dataCollection.deviceState().kineticEnergy,nParticles);
		
		mpd::zeroRange_device(dataCollection.deviceState().potentialEnergy,nParticles);
		mpd::bondPotential_device(bondList.deviceInteraction(),state.deviceState(),
					  dataCollection.deviceState());
		for(auto &ballList:ballLists)
			mpd::ballPotential_device(ballList.deviceInteraction(),state.deviceState(),
					  dataCollection.deviceState());
		mpd::bendPotential_device(bendList.deviceInteraction(),state.deviceState(),
					  dataCollection.deviceState());
		mpd::cellComputePotential_device(cData.deviceCell(),state.deviceState(),
					  dataCollection.deviceState());
		float potential=
			mpd::reduceRange_device(dataCollection.deviceState().potentialEnergy,nParticles);
		
		std::fstream potentialFile(potentialFileName, std::ios::out | std::ios::app);
		potentialFile << "0\t" << potential << std::endl;
		
		std::fstream kineticFile(kineticFileName, std::ios::out | std::ios::app);
		kineticFile << "0\t" << kinetic << std::endl;
		
		xyzFile.open(framesFilename,std::ios::out | std::ios::app);
		xyzFile.store();
		xyzFile.close();
	}
	else
	{
		//Surprise! This is done because a previously run configuration don't do this upon exit
		mpd::verletSecond_device(state.deviceState());
	}
	//using integer indexing, the 0.0000001 fixes an accuracy issue with the gnu c++ compiler.
	//Don't believe it affects anything else...
	int endInt=int(System.readFinalTime()/System.readDeltaT()+0.0000001);//end
	int startInt=int(System.readInitialTime()/System.readDeltaT()+0.0000001);//start
	
	int storeint=int(System.readStoreInterval()/System.readDeltaT()+0.0000001);//when to store
	int measureint=int(System.readMeasureInterval()/System.readDeltaT()+0.0000001);//when to measure
	
	
	int tempStepInt=0;
	float tempStep=0;
	if(System.readTempStepInterval()>0)
	{
		tempStep=System.readTempStepInterval()*(System.readFinalTemp()-System.readInitialTemp())/
			 (System.readFinalTime()-System.readInitialTime());
		tempStepInt=(endInt-startInt)/int((System.readFinalTime()-System.readInitialTime())/
			    (System.readTempStepInterval()));
	}
	
	std::cerr << "starting main loop: \n";
	
	time_t current=time(NULL);
	
	bool exitFlag=false;//for premature exits!
	
	//the molecular dynamics loop, the "running" of the system
	for(int i=startInt;i<=endInt && !exitFlag;i++)
	{
		System.setInitialTime((float)i*System.readDeltaT());
		
		mpd::verletFirst_device(state.deviceState());
		mpd::zeroRange_device(state.deviceState().a,nParticles);
		//The system is stored here because force is updated here, but velocities are updated next.
		//It causes a problem when it reenters the loop from a previously run configuration.
		if(System.readTempStepInterval()>0)
			if(i%tempStepInt==0 && tempStepInt!=0 && i<endInt)
				System.setInitialTemp(System.readInitialTemp()+tempStep);
		
		if(i%storeint==0 && i!=startInt)
		{
			state.toHost();
			fileIO.open(nameA.c_str(),std::ios::out);
			fileIO.write();
			fileIO.close();
			
			fileIO.open(name,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
		}
		state.temperature=System.readInitialTemp();
		mpd::langevin_device(state.deviceState(),random.deviceState());
		mpd::cellComputeForce_device(cData.deviceCell(), state.deviceState());
		mpd::bondForces_device(bondList.deviceInteraction(),state.deviceState());
		for(auto &ballList:ballLists)
			mpd::ballForces_device(ballList.deviceInteraction(),state.deviceState());
		mpd::bendForces_device(bendList.deviceInteraction(),state.deviceState());
		mpd::applyMass_device(state.deviceState());
		mpd::verletSecond_device(state.deviceState());
		
		//Measurements are output here
		if(i%measureint==0 && i!=startInt)
		{
			//float last=current;
			time_t last=current;
			current=time(NULL);
			//time since last storage step, good for benchmarking
			std::cerr << System.readInitialTime() << '\t' << current-last << std::endl;
			
			//Data calculations that we are interested in starting
			
			//Data calculations
			mpd::zeroRange_device(dataCollection.deviceState().kineticEnergy,nParticles);
			mpd::kinetic_device(state.deviceState(),dataCollection.deviceState());
			float kinetic=
				mpd::reduceRange_device(dataCollection.deviceState().kineticEnergy,nParticles);
			
			mpd::zeroRange_device(dataCollection.deviceState().potentialEnergy,nParticles);
			mpd::bondPotential_device(bondList.deviceInteraction(),state.deviceState(),
						  dataCollection.deviceState());
			for(auto &ballList:ballLists)
				mpd::ballPotential_device(ballList.deviceInteraction(),state.deviceState(),
						  dataCollection.deviceState());
			mpd::bendPotential_device(bendList.deviceInteraction(),state.deviceState(),
						  dataCollection.deviceState());
			mpd::cellComputePotential_device(cData.deviceCell(),state.deviceState(),
						  dataCollection.deviceState());
			float potential=
				mpd::reduceRange_device(dataCollection.deviceState().potentialEnergy,nParticles);
			
			//dataCollection.toHost();
			//state.toHost();
			
			std::fstream potentialFile(potentialFileName, std::ios::out | std::ios::app);
			potentialFile << System.readInitialTime() << "\t" << potential << std::endl;
			
			std::fstream kineticFile(kineticFileName, std::ios::out | std::ios::app);
			kineticFile << System.readInitialTime() << "\t" << kinetic << std::endl;
			
			threeVector<float> size=System.readSize();
			std::fstream sizeFile(sizeFileName, std::ios::out | std::ios::app);
			sizeFile << System.readInitialTime() << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			current=time(NULL);
		}
		//short section to resize system, note that it only works when deltaLXY 
		// is something other than 0, it flags execution.
		if(i%resizeRate==0 && i!=0 && System.readDeltaLXY()!=0)
		{
			threeVector<float> oldSize=System.readSize();
			//threeVector<float> fluctuation=bStat.dependentFluctuationConstantVolume(randGen, oldSize);
			threeVector<float> fluctuation=bStat.independentFluctuationConstantVolume(randGen, oldSize);
			
			threeVector<float> newSize=oldSize;
			newSize.x+=fluctuation.x;
			newSize.y+=fluctuation.y;
			newSize.z+=fluctuation.z;
			threeVector<float> scale=newSize;
			scale.x/=oldSize.x;
			scale.y/=oldSize.y;
			scale.z/=oldSize.z;
			
			mpd::zeroRange_device(bStat.deviceState().dPotential,nParticles);
			mpd::bondDPotential_device(bondList.deviceInteraction(),state.deviceState(),
						  bStat.deviceState(),scale);
			for(auto &ballList:ballLists)
				mpd::ballDPotential_device(ballList.deviceInteraction(),state.deviceState(),
						  bStat.deviceState(),scale);
			mpd::bendDPotential_device(bendList.deviceInteraction(),state.deviceState(),
						  bStat.deviceState(),scale);
			mpd::cellComputeDPotential_device(cData.deviceCell(),state.deviceState(),
						  bStat.deviceState(),scale);
			float dPotential=
				mpd::reduceRange_device(bStat.deviceState().dPotential,nParticles);
			if(bStat.MCtest(dPotential, System.readInitialTemp(), 0.0, 
				System.readTension(),oldSize,newSize, randGen))
			{
				//this can't perform the state resize as of yet, but can rescale the 
				//device positions
				mpd::rescale_device(state.deviceState(),scale);
				//this should be the only resize for everything
				state.resize(newSize);
				//this could be deleted if the above is universal
				cData.resize(newSize,System.readDeltaLXY());
				System.setSize(newSize);
				resizeHist[(fluctuation.x+maxInterval)/resizeHistInterval]+=0.5;
				resizeHist[(fluctuation.y+maxInterval)/resizeHistInterval]+=0.5;
				//int rIndex=(fluctuation.x+System.readDeltaLXY())/resizeHistInterval;
				
				accepted++;
			}
			else
			{
				rejectHist[(fluctuation.x+maxInterval)/resizeHistInterval]+=0.5;
				rejectHist[(fluctuation.y+maxInterval)/resizeHistInterval]+=0.5;
			}
			trial++;
		}
	}
	
	if(System.readDeltaLXY()!=0)
	{	
		std::fstream resizeHistFile;
		std::string resizeHistFileName="resizeHist_";
		resizeHistFileName+=argv[1];
		resizeHistFileName+=".dat";
		
		resizeHistFile.open(resizeHistFileName.c_str(),std::ios::out);
		if(resizeHistFile.is_open())
		{
			for(int i=0;i<resizeHist.size();i++)
				resizeHistFile << (static_cast<float>(i)*resizeHistInterval)-maxInterval\
					 << '\t' << resizeHist[i] << '\t' << rejectHist[i] << std::endl;
			resizeHistFile.close();
		}
		std::cerr << "Resize acceptance ratio: " << accepted/trial << std::endl;
	}
	
	return 0;
}

