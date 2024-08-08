//This simply runs a system through the timesteps. It can use OpenMP.
//There was an issue with setting a pointer to the Blob type and allocation, so
// Blob type contains it's own pointers now. Everything is being rewritten for this.

//For molecular dynamics forces and potentials
#include "include/DPD.h"

//For the molecular dynamics variables
#include "include/system.h"
#include <ctime>

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: command name\n";
		return 0;
	}
	
	char *name=argv[1];
	
	///the variables for the simulation, remember that for some reason constructor isn't explicit when nothing is in it
	Blob<double> System;
	
	///load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	///initialize the variables, objects, algorithms, data collection, etc...
	Verlet<double> integrate(System.getPositions(), System.getAccelerations(), System.getVelocities(), System.readNParticles(),
				 System.readSize(), System.readDeltaT(), System.readPeriodic());
	Kinetic<double> kinetic(System.getVelocities(), System.readNParticles());
	position <double> *xyzPositions=System.getPositions();
	int xyzNParticles=System.readNParticles();
	xyzFormat<double> xyzFile(xyzPositions, xyzNParticles);
	MTRand randNum(System.readSeed());
	
	CellOpt<double, ForceDPD <double> > pairInteractions(System.getPositions(), System.getVelocities(),
		System.getAccelerations(), System.getTwoBodyFconst(), System.readNParticles(),
		System.readNTypes(), System.readSize(), System.readPeriodic(), System.readCutoff(), System.readSeed());
	
	///Other initializations
	std::string framesFilename("frames_");
	framesFilename+=name;
	framesFilename+=".xyz";
	
	for(int i=0;i<System.readNParticles();i++)
	{
		System.getAccelerations()[i].x=0;
		System.getAccelerations()[i].y=0;
		System.getAccelerations()[i].z=0;
	}
	
	///this is how you call the pair interactions, build then compute...
	///Cuda works for these
	pairInteractions.build();
	pairInteractions.computeForce();
	
	
	///Cuda doesn't work for this loop yet
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
	
	//this corrects an issue where an extra data point is added when the system is restarted
	if(System.readInitialTime()==0)
	{
		///Cuda works for this
		double potential=0;
		
		///Cuda does not work on this loop yet
		//Better version of molecule interactions
		double lBond=0;
		int nBond=0;
		//Go through all molecule structures
		for(int k=0;k<System.readNMolecules();k++)
		{
			//pick a structure by type
			switch(System.getMolecule()[k].readType())
			{
				case BOND:
				{
					potential+=System.doBondPotential(k);
					//Go through all in bond list
					double lBondPrivate=0;
					//#pragma omp parallel for reduction(+:lBondPrivate)
					for(int l=0;l<System.getMolecule()[k].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=System.getMolecule()[k].getBonds()[l].s[0];
						int secondParticle=System.getMolecule()[k].getBonds()[l].s[1];
						
						//calculate bond length using minimum image
						double dx=System.getPositions()[firstParticle].x-System.getPositions()[secondParticle].x;
						double dy=System.getPositions()[firstParticle].y-System.getPositions()[secondParticle].y;
						double dz=System.getPositions()[firstParticle].z-System.getPositions()[secondParticle].z;
						
						if(dx>=System.readSize().x/2.0) dx-=System.readSize().x;
						if(dx<=-System.readSize().x/2.0) dx+=System.readSize().x;
						if(dy>=System.readSize().y/2.0) dy-=System.readSize().y;
						if(dy<=-System.readSize().y/2.0) dy+=System.readSize().y;
						if(dz>=System.readSize().z/2.0) dz-=System.readSize().z;
						if(dz<=-System.readSize().z/2.0) dz+=System.readSize().z;
						
						lBondPrivate+=sqrt(dx*dx+dy*dy+dz*dz);
					}
					lBond+=lBondPrivate;
					nBond+=System.getMolecule()[k].readNBond();
					break;
				}
				case BEND:
				{
					potential+=System.doBendPotential(k);
					break;
				}
				case CHAIN:
				{
					potential+=System.doChainPotential(k);
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
		}
		
		std::fstream dataFile;
		std::string buf("potential_");
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << System.readInitialTime() << '\t' << potential << std::endl;
		dataFile.close();
		
		//Just copy this one and use it as a template
		buf.clear();
		buf="lBond_";
		buf+=name;
		buf+=".dat";
		dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
		dataFile << System.readInitialTime() << '\t' << (lBond/(double)nBond) << std::endl;
		dataFile.close();
		
		kinetic.output(System.readInitialTime(),name);
		
		xyzFile.open(framesFilename,std::ios::out | std::ios::app);
		xyzFile.store();
		xyzFile.close();
	}
	
	//using integer indexing, the 0.0000001 fixes an accuracy issue with the gnu c++ compiler.
	//Don't believe it affects anything else...
	int endInt=int(System.readFinalTime()/System.readDeltaT()+0.0000001);//end
	int startInt=int(System.readInitialTime()/System.readDeltaT()+0.0000001);//start
	
	int storeint=int(System.readStoreInterval()/System.readDeltaT()+0.0000001);//when to store
	int measureint=int(System.readMeasureInterval()/System.readDeltaT()+0.0000001);//when to measure	
	
	std::cout << "starting main loop: \n";
	time_t current=time(NULL);
	
	///the dissipative particle dynamics loop, the "running" of the system
	for(int i=startInt;i<=endInt;i++)
	{
		System.setInitialTime((double)i*System.readDeltaT());
		
		integrate.first();
		
		for(int k=0;k<System.readNParticles();k++)
		{
			System.getAccelerations()[k].x=0;
			System.getAccelerations()[k].y=0;
			System.getAccelerations()[k].z=0;
		}
		
		///Cuda works for these
		
		//Build linked lists
		pairInteractions.build();
		
		pairInteractions.computeForce();
		
		///Cuda doesn't work for this loop yet
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
		
		integrate.second();
		
		//threeVector<double> momentum=pairwiseSum< threeVector<double> >(System.getVelocities(), 0, System.readNParticles());
		//for(int k=0;k<System.readNParticles();k++)
		//	System.getVelocities()[i]-=momentum;
		
		threeVector<double> momentum;
		momentum.x=0;
		momentum.y=0;
		momentum.z=0;
		for(int k=0;k<System.readNParticles();k++)
		{
			momentum.x+=System.getVelocities()[k].x;
			momentum.y+=System.getVelocities()[k].y;
			momentum.z+=System.getVelocities()[k].z;
		}
		
		momentum.x/=(double)System.readNParticles();
		momentum.y/=(double)System.readNParticles();
		momentum.z/=(double)System.readNParticles();
		
		//std::cout << momentum.x << ' ' << momentum.y << ' ' << momentum.z << '\n';
		//std::cin.get();
		
		for(int k=0;k<System.readNParticles();k++)
		{
			System.getVelocities()[k].x-=momentum.x;
			System.getVelocities()[k].y-=momentum.y;
			System.getVelocities()[k].z-=momentum.z;
		}
		
		if(i%storeint==0 && i!=startInt)
		{
			fileIO.open(name,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
		}
		
		if(i%measureint==0 && i!=startInt)
		{
			//double last=current;
			time_t last=current;
			//current=omp_get_wtime();
			current=time(NULL);
			//time since last storage step, good for benchmarking
			std::cout << System.readInitialTime() << '\t' << current-last << std::endl;
			
			double potential=0;
			
			///Cuda does not work on this loop yet
			//Better version of molecule interactions
			double lBond=0;
			int nBond=0;
			//Go through all molecule structures
			for(int k=0;k<System.readNMolecules();k++)
			{
				//pick a structure by type
				switch(System.getMolecule()[k].readType())
				{
					case BOND:
					{
						potential+=System.doBondPotential(k);
						//Go through all in bond list in parallel using openmp
						double lBondPrivate=0;
						//#pragma omp parallel for reduction(+:lBondPrivate)
						for(int l=0;l<System.getMolecule()[k].readNBond();l++)
						{
							//These are the first and second particles of the bond
							int firstParticle=System.getMolecule()[k].getBonds()[l].s[0];
							int secondParticle=System.getMolecule()[k].getBonds()[l].s[1];
							
							//calculate bond length using minimum image
							double dx=System.getPositions()[firstParticle].x-System.getPositions()[secondParticle].x;
							double dy=System.getPositions()[firstParticle].y-System.getPositions()[secondParticle].y;
							double dz=System.getPositions()[firstParticle].z-System.getPositions()[secondParticle].z;
							
							if(dx>=System.readSize().x/2.0) dx-=System.readSize().x;
							if(dx<=-System.readSize().x/2.0) dx+=System.readSize().x;
							if(dy>=System.readSize().y/2.0) dy-=System.readSize().y;
							if(dy<=-System.readSize().y/2.0) dy+=System.readSize().y;
							if(dz>=System.readSize().z/2.0) dz-=System.readSize().z;
							if(dz<=-System.readSize().z/2.0) dz+=System.readSize().z;
							
							lBondPrivate+=sqrt(dx*dx+dy*dy+dz*dz);
						}
						lBond+=lBondPrivate;
						nBond+=System.getMolecule()[k].readNBond();
						break;
					}
					case BEND:
					{
						potential+=System.doBendPotential(k);
						break;
					}
					case CHAIN:
					{
						potential+=System.doChainPotential(k);
						break;
					}
					default:
						//does nothing
						break;
				}
			}
			
			std::fstream dataFile;
			std::string buf("potential_");
			buf+=name;
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			dataFile << System.readInitialTime() << '\t' << potential << std::endl;
			dataFile.close();
			
			//Just copy this one and use it as a template
			buf.clear();
			buf="lBond_";
			buf+=name;
			buf+=".dat";
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			dataFile << System.readInitialTime() << '\t' << (lBond/(double)nBond) << std::endl;
			dataFile.close();
			
			kinetic.output(System.readInitialTime(),name);
			
			//This really isn't accurate without this here
			current=time(NULL);
			//if you want to add a new measure when you do a measure step, add it here
		}
	}
	
	return 0;
}
