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

//Include files from the library:

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <ctime>




//frequency that the system is resized, 1=every time step, 2=every other time step, 3=every 3rd time step, etc...
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
	
	//for easy (and quick) access to certain variables
	position <double> *p=System.getPositions();
	threeVector<double> *a=System.getAccelerations();
	threeVector<double> *v=System.getVelocities();
	int nParticles=System.readNParticles();
	
	double *cF=System.getTwoBodyFconst();//constants for pairwise force interactions
	
	int nTypes=System.readNTypes();
	
	double cutoff=System.readCutoff();
	double cutoffSquared=cutoff*cutoff;
	threeVector<double> size=System.readSize();
	
	molecule<double, fourVector<int> > *m=System.getMolecule();
	int nMolecules=System.readNMolecules();
	
	double temperature=System.readInitialTemp();
	double deltaT=System.readDeltaT();
	double gamma=System.readGamma();
	double sigma=sqrt((6.0*temperature*gamma)/deltaT);
	double halfDt=deltaT*0.5;
	
	
	
	
	//for xyz frames
	xyzFormat<double> xyzFile(p, nParticles);
	
	
	//for random numbers
	MTRand randNum(System.readSeed());
	
	//Other initializations
	std::string framesFilename("frames_");
	framesFilename+=name;
	framesFilename+=".xyz";
	
	
	
	
	
	
	
	
	
	
	//Begin pre-computations to obtain the initial accelerations
	
	//Initialize all accelerations to 0
	for(int k=0;k<nParticles;k++)
		a[k]=0;
	
	//Compute the forces due to the implicit solvent bath.
	//Simplified version of Langevin::compute()
	for(int i=0;i<nParticles;i++)
	{
		double psi;
		psi=2.0*randNum.rand53()-1.0;
		a[i].x+=(-gamma*v[i].x+sigma*psi);
		psi=2.0*randNum.rand53()-1.0;
		a[i].y+=(-gamma*v[i].y+sigma*psi);
		psi=2.0*randNum.rand53()-1.0;
		a[i].z+=(-gamma*v[i].z+sigma*psi);
	}
	
	//Compute the pairwise forces for all pairs (very, very slow).
	//Simplified version of cellOpt::computeForce()
	for(int i=0;i<nParticles-1;i++)
	{
		//this is for summing the constituent accelerations of a particle
		threeVector<double> ai=0;
		
		for(int j=i+1;j<nParticles;j++)
		{
			//temporary position containers
			position<double> pi=p[i];
			position<double> pj=p[j];
			
			//check for the minimum image for each coordinate in each direction
			if(pi.x-pj.x>size.x/2.0)
				pj.x+=size.x;
			if(pi.x-pj.x<-size.x/2.0)
				pj.x-=size.x;
			
			if(pi.y-pj.y>size.y/2.0)
				pj.y+=size.y;
			if(pi.y-pj.y<-size.y/2.0)
				pj.y-=size.y;
			
			if(pi.z-pj.z>size.z/2.0)
				pj.z+=size.z;
			if(pi.z-pj.z<-size.z/2.0)
				pj.z-=size.z;
			
			//Computes the forces, note that the acceleration of the i-th particle (ai)
			// is a variable, while the acceleration of the j-th particle (a[j]) is an 
			// index of the pointer a. This allows newton's third law to be computed
			// explicitly, rather than implicitly via an additional aj variable.
			//The definition (near the end of systemMD) is:
			//	template <typename T>
			//	inline void Force(
			//		T cutoffSquared, 
			//		int nT, 
			//		position<T> &p1, //note the reference
			//		position<T> &p2, //note the reference
			//		threeVector<T> &a1, //note the reference
			//		threeVector<T> &a2, //note the reference
			//		T *fC) //note the pointer
			//Note that the force definition doesn't include the minimum image, 
			// so we will implement it here. Also, note that references are modifiable.
			Force(cutoffSquared, nTypes, pi, pj, ai, a[j], cF);
		}
		
		//add to the acceleration of the i-th particle
		a[i]+=ai;
	}
	
	//Compute the molecular interaction forces.
	for(int k=0;k<nMolecules;k++)
	{
		switch(m[k].readType())
		{
			case BOND:
			{
				//This is present in the System variable
				//System.doBondForce(k);
				
				// But it explicitly does this:
				fourVector<int> *bond=m[k].getBonds();
				
				for(int j=0;j<m[k].readNBond();j++)
				{
					//The indices for the first and second particles.
					int first=bond[j].s[0];
					int second=bond[j].s[1];
					
					//Distance between particles
					threeVector<double> d;
					d.x=p[first].x-p[second].x;
					d.y=p[first].y-p[second].y;
					d.z=p[first].z-p[second].z;
					
					//Minimum image distances. Note they are compacted slightly.
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					
					//Computes the harmonic force, note that this one uses
					// a slightly different variation of force. Instead of allowing 
					// p[first] and p[second] to vary, we only calculate the absolute distance
					// between the two. This is quick and easy for a "System" function,
					// but if we wanted a more complex interaction between "first" and 
					// "second", we would need to change the definition. The added flexibility
					// of having p[i] and p[j] for the force variant of pairwise 
					// interactions is required if you wanted to add or remove dimensions or
					// degrees of freedom.
					//The definition is:
					//	template <typename T>
					//	void Blob::harmonicF(threeVector<T> d, 
					//		threeVector<T> &a1, //note the reference
					//		threeVector<T> &a2, //note the reference
					//		T *constants) //note the pointer
					System.harmonicF(d, a[first], a[second], m[k].getConstants());
				}
				
				break;
			}
			case BEND:
			{
				//This is present in the System variable
				//System.doBendForce(k);
				
				// But it explicitly does this:
				fourVector<int> *bond=m[k].getBonds();
				
				for(int j=0;j<m[k].readNBond();j++)
				{
					//The indices for the first, second, and third particles.
					int first=bond[j].s[0];
					int second=bond[j].s[1];
					int third=bond[j].s[2];
					
					//Distances between particles. Notice the excessive repitition.
					threeVector<double> da,db;
					
					da.x=p[first].x-p[second].x;
					da.y=p[first].y-p[second].y;
					da.z=p[first].z-p[second].z;
					
					db.x=p[second].x-p[third].x;
					db.y=p[second].y-p[third].y;
					db.z=p[second].z-p[third].z;
					
					//Minimum image distances.
					if(da.x>size.x/2.0) da.x-=size.x;
					if(da.x<-size.x/2.0) da.x+=size.x;
					if(da.y>size.y/2.0) da.y-=size.y;
					if(da.y<-size.y/2.0) da.y+=size.y;
					if(da.z>size.z/2.0) da.z-=size.z;
					if(da.z<-size.z/2.0) da.z+=size.z;
					
					if(db.x>size.x/2.0) db.x-=size.x;
					if(db.x<-size.x/2.0) db.x+=size.x;
					if(db.y>size.y/2.0) db.y-=size.y;
					if(db.y<-size.y/2.0) db.y+=size.y;
					if(db.z>size.z/2.0) db.z-=size.z;
					if(db.z<-size.z/2.0) db.z+=size.z;
					
					//Computes the bending force. Note the added complexity 
					// due to tha additional particles.
					//The definition is:
					//	template <typename T>
					//	void Blob<T>::bendF(threeVector<T> da, 
					//		threeVector<T> db, 
					//		threeVector<T> &a1, //note the reference
					//		threeVector<T> &a2, //note the reference
					//		threeVector<T> &a3, //note the reference
					//		T *constants) //note the pointer
					System.bendF(da, db, a[first], a[second], a[third], m[k].getConstants());
				}
				
				break;
			}
			case CHAIN:
			{
				//This is present in the System variable
				//System.doChainForce(k);
				
				// But it explicitly does this:
				fourVector<int> *bond=m[k].getBonds();
				double *cBond=m[k].getConstants()+CHAINBOND;
				double *cBend=m[k].getConstants()+CHAINBEND;
				
				for(int j=0;j<m[k].readNBond();j++)
				{
					//bonding information
					int start=bond[j].s[START];
					int nChains=bond[j].s[NCHAINS];
					int length=bond[j].s[LENGTH];
					
					//go through all chain lengths
					for(int k=start; k<start+length*nChains; k+=length)
					{
						threeVector<double> da,db;
						//go through all particles in chain
						for(int l=k;l<k+length-3;l++)
						{
							//Getting the indices via simple additions.
							int first=l;
							int second=l+1;
							int third=l+2;
							
							//Distances between particles
							da.x=p[first].x-p[second].x;
							da.y=p[first].y-p[second].y;
							da.z=p[first].z-p[second].z;
							
							db.x=p[second].x-p[third].x;
							db.y=p[second].y-p[third].y;
							db.z=p[second].z-p[third].z;
							
							//Minimum image distances.
							if(da.x>size.x/2.0) da.x-=size.x;
							if(da.x<-size.x/2.0) da.x+=size.x;
							if(da.y>size.y/2.0) da.y-=size.y;
							if(da.y<-size.y/2.0) da.y+=size.y;
							if(da.z>size.z/2.0) da.z-=size.z;
							if(da.z<-size.z/2.0) da.z+=size.z;
							
							if(db.x>size.x/2.0) db.x-=size.x;
							if(db.x<-size.x/2.0) db.x+=size.x;
							if(db.y>size.y/2.0) db.y-=size.y;
							if(db.y<-size.y/2.0) db.y+=size.y;
							if(db.z>size.z/2.0) db.z-=size.z;
							if(db.z<-size.z/2.0) db.z+=size.z;
							
							//These are exactly the same as the above definitions of bonded forces.
							System.harmonicF(da, a[first], a[second], cBond);
							System.bendF(da, db, a[first], a[second], a[third], cBend);
						}
						//last three are a special case, requires bond between last two particles in chain
						int first=k+length-3;
						int second=k+length-2;
						int third=k+length-1;
						
						//Distances between particles
						da.x=p[first].x-p[second].x;
						da.y=p[first].y-p[second].y;
						da.z=p[first].z-p[second].z;
						
						db.x=p[second].x-p[third].x;
						db.y=p[second].y-p[third].y;
						db.z=p[second].z-p[third].z;
						
						//Minimum image distances.
						if(da.x>size.x/2.0) da.x-=size.x;
						if(da.x<-size.x/2.0) da.x+=size.x;
						if(da.y>size.y/2.0) da.y-=size.y;
						if(da.y<-size.y/2.0) da.y+=size.y;
						if(da.z>size.z/2.0) da.z-=size.z;
						if(da.z<-size.z/2.0) da.z+=size.z;
						
						if(db.x>size.x/2.0) db.x-=size.x;
						if(db.x<-size.x/2.0) db.x+=size.x;
						if(db.y>size.y/2.0) db.y-=size.y;
						if(db.y<-size.y/2.0) db.y+=size.y;
						if(db.z>size.z/2.0) db.z-=size.z;
						if(db.z<-size.z/2.0) db.z+=size.z;
						
						//More of the same definitions of the bonded forces.
						System.harmonicF(da, a[first], a[second], cBond);
						System.harmonicF(db, a[second], a[third], cBond);
						System.bendF(da, db, a[first], a[second], a[third], cBend);
					}
				}
				
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
		xyzFile.open(framesFilename,std::ios::out | std::ios::app);
		xyzFile.store();
		xyzFile.close();
	}
	else
	{
		//Surprise! This is done because a previously run configuration doesn't do this upon exit
		//Simplified version of Verlet::second()
		for(int i=0;i<nParticles;i++)
		{
			v[i].x+=(a[i].x*halfDt);
			v[i].y+=(a[i].y*halfDt);
			v[i].z+=(a[i].z*halfDt);
		}
	}
	
	//using integer indexing, the 0.0000001 fixes an accuracy issue with the gnu c++ compiler.
	//Don't believe it affects anything else...
	int endInt=int(System.readFinalTime()/System.readDeltaT()+0.0000001);//end
	int startInt=int(System.readInitialTime()/System.readDeltaT()+0.0000001);//start
	
	int storeint=int(System.readStoreInterval()/System.readDeltaT()+0.0000001);//when to store
	int measureint=int(System.readMeasureInterval()/System.readDeltaT()+0.0000001);//when to measure
	
	
	
	std::cout << "starting main loop: \n";
	
	time_t current=time(NULL);

	//The molecular dynamics loop, the "running" of the system. Notice the shear number of
	// repeated lines! If you see the unsimplified versions of all of these you will notice that
	// they would make this much, much longer.
	for(int i=startInt;i<=endInt;i++)
	{
		//Increment the simulation time by a time step.
		System.setInitialTime(static_cast<double>(i)*System.readDeltaT());
		
		//Compute the first half of the integration
		//Simplified version of Verlet::first()
		for(int k=0;k<nParticles;k++)
		{
			//increment velocity
			v[k].x+=(a[k].x*halfDt);
			v[k].y+=(a[k].y*halfDt);
			v[k].z+=(a[k].z*halfDt);
			
			//increment position
			p[k].x+=v[k].x*deltaT;
			p[k].y+=v[k].y*deltaT;
			p[k].z+=v[k].z*deltaT;
			
			//If it crosses the boundaries, wrap it to the other side.
			// "Periodic boundaries"
			if(p[k].x>size.x) p[k].x-=size.x;
			if(p[k].x<0) p[k].x+=size.x;
			
			if(p[k].y>size.y) p[k].y-=size.y;
			if(p[k].y<0) p[k].y+=size.y;
			
			if(p[k].z>size.z) p[k].z-=size.z;
			if(p[k].z<0) p[k].z+=size.z;
		}
		
		//Zero out all the accelerations.
		for(int k=0;k<nParticles;k++)
			a[k]=0;
		
		//Check to see if we need to store the current configuration.
		if(i%storeint==0 && i!=startInt)
		{
			fileIO.open(name,std::ios::out);
			fileIO.write();
			fileIO.close();
			
			xyzFile.open(framesFilename,std::ios::out | std::ios::app);
			xyzFile.store();
			xyzFile.close();
		}
		
		//Compute the forces due to the implicit solvent bath.
		//Simplified version of Langevin::compute()
		for(int k=0;k<nParticles;k++)
		{
			double psi;
			psi=2.0*randNum.rand53()-1.0;
			a[k].x+=(-gamma*v[k].x+sigma*psi);
			psi=2.0*randNum.rand53()-1.0;
			a[k].y+=(-gamma*v[k].y+sigma*psi);
			psi=2.0*randNum.rand53()-1.0;
			a[k].z+=(-gamma*v[k].z+sigma*psi);
		}
		
		//Compute the pairwise forces for all pairs (very, very slow).
		//Simplified version of cellOpt::computeForce()
		for(int k=0;k<nParticles-1;k++)
		{
			//this is for summing the constituent accelerations of a particle
			threeVector<double> ak=0;
			
			for(int j=k+1;j<nParticles;j++)
			{
				//temporary position containers
				position<double> pk=p[k];
				position<double> pj=p[j];
				
				//check for the minimum image for each coordinate in each direction
				if(pk.x-pj.x>size.x/2.0) pj.x+=size.x;
				if(pk.x-pj.x<-size.x/2.0) pj.x-=size.x;
				
				if(pk.y-pj.y>size.y/2.0) pj.y+=size.y;
				if(pk.y-pj.y<-size.y/2.0) pj.y-=size.y;
				
				if(pk.z-pj.z>size.z/2.0) pj.z+=size.z;
				if(pk.z-pj.z<-size.z/2.0) pj.z-=size.z;
				
				//Computes the force.
				Force(cutoffSquared, nTypes, pk, pj, ak, a[j], cF);
			}
			
			//add to the acceleration of the i-th particle
			a[k]+=ak;
		}
		
		//Compute the molecular interaction forces.
		for(int k=0;k<System.readNMolecules();k++)
		{
			switch(System.getMolecule()[k].readType())
			{
				case BOND:
				{
					//This is present in the System variable
					//System.doBondForce(k);
					
					// But it explicitly does this:
					fourVector<int> *bond=m[k].getBonds();
					
					for(int j=0;j<m[k].readNBond();j++)
					{
						//The indices for the first and second particles.
						int first=bond[j].s[0];
						int second=bond[j].s[1];
						
						//Distance between particles
						threeVector<double> d;
						d.x=p[first].x-p[second].x;
						d.y=p[first].y-p[second].y;
						d.z=p[first].z-p[second].z;
						
						//Minimum image distances. Note they are compacted slightly.
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						
						//Computes the harmonic force.
						System.harmonicF(d, a[first], a[second], m[k].getConstants());
					}
					
					break;
				}
				case BEND:
				{
					//This is present in the System variable
					//System.doBendForce(k);
					
					// But it explicitly does this:
					fourVector<int> *bond=m[k].getBonds();
					
					for(int j=0;j<m[k].readNBond();j++)
					{
						//The indices for the first, second, and third particles.
						int first=bond[j].s[0];
						int second=bond[j].s[1];
						int third=bond[j].s[2];
						
						//Distances between particles. Notice the excessive repitition.
						threeVector<double> da,db;
						
						da.x=p[first].x-p[second].x;
						da.y=p[first].y-p[second].y;
						da.z=p[first].z-p[second].z;
						
						db.x=p[second].x-p[third].x;
						db.y=p[second].y-p[third].y;
						db.z=p[second].z-p[third].z;
						
						//Minimum image distances.
						if(da.x>size.x/2.0) da.x-=size.x;
						if(da.x<-size.x/2.0) da.x+=size.x;
						if(da.y>size.y/2.0) da.y-=size.y;
						if(da.y<-size.y/2.0) da.y+=size.y;
						if(da.z>size.z/2.0) da.z-=size.z;
						if(da.z<-size.z/2.0) da.z+=size.z;
						
						if(db.x>size.x/2.0) db.x-=size.x;
						if(db.x<-size.x/2.0) db.x+=size.x;
						if(db.y>size.y/2.0) db.y-=size.y;
						if(db.y<-size.y/2.0) db.y+=size.y;
						if(db.z>size.z/2.0) db.z-=size.z;
						if(db.z<-size.z/2.0) db.z+=size.z;
						
						//Computes the bending force.
						System.bendF(da, db, a[first], a[second], a[third], m[k].getConstants());
					}
					
					break;
				}
				case CHAIN:
				{
					//This is present in the System variable
					//System.doChainForce(k);
					
					// But it explicitly does this:
					fourVector<int> *bond=m[k].getBonds();
					double *cBond=m[k].getConstants()+CHAINBOND;
					double *cBend=m[k].getConstants()+CHAINBEND;
					
					for(int j=0;j<m[k].readNBond();j++)
					{
						//bonding information
						int start=bond[j].s[START];
						int nChains=bond[j].s[NCHAINS];
						int length=bond[j].s[LENGTH];
						
						
						//go through all chain lengths
						for(int k=start; k<start+length*nChains; k+=length)
						{
							threeVector<double> da,db;
							//go through all particles in chain
							for(int l=k;l<k+length-3;l++)
							{
								//Getting the indices via simple additions.
								int first=l;
								int second=l+1;
								int third=l+2;
								
								//Distances between particles
								da.x=p[first].x-p[second].x;
								da.y=p[first].y-p[second].y;
								da.z=p[first].z-p[second].z;
								
								db.x=p[second].x-p[third].x;
								db.y=p[second].y-p[third].y;
								db.z=p[second].z-p[third].z;
								
								//Minimum image distances.
								if(da.x>size.x/2.0) da.x-=size.x;
								if(da.x<-size.x/2.0) da.x+=size.x;
								if(da.y>size.y/2.0) da.y-=size.y;
								if(da.y<-size.y/2.0) da.y+=size.y;
								if(da.z>size.z/2.0) da.z-=size.z;
								if(da.z<-size.z/2.0) da.z+=size.z;
								
								if(db.x>size.x/2.0) db.x-=size.x;
								if(db.x<-size.x/2.0) db.x+=size.x;
								if(db.y>size.y/2.0) db.y-=size.y;
								if(db.y<-size.y/2.0) db.y+=size.y;
								if(db.z>size.z/2.0) db.z-=size.z;
								if(db.z<-size.z/2.0) db.z+=size.z;
								
								//These are exactly the same as the above definitions of bonded forces.
								System.harmonicF(da, a[first], a[second], cBond);
								System.bendF(da, db, a[first], a[second], a[third], cBend);
							}
							//last three are a special case, requires bond between last two particles in chain
							int first=k+length-3;
							int second=k+length-2;
							int third=k+length-1;
							
							//Distances between particles
							da.x=p[first].x-p[second].x;
							da.y=p[first].y-p[second].y;
							da.z=p[first].z-p[second].z;
							
							db.x=p[second].x-p[third].x;
							db.y=p[second].y-p[third].y;
							db.z=p[second].z-p[third].z;
							
							//Minimum image distances.
							if(da.x>size.x/2.0) da.x-=size.x;
							if(da.x<-size.x/2.0) da.x+=size.x;
							if(da.y>size.y/2.0) da.y-=size.y;
							if(da.y<-size.y/2.0) da.y+=size.y;
							if(da.z>size.z/2.0) da.z-=size.z;
							if(da.z<-size.z/2.0) da.z+=size.z;
							
							if(db.x>size.x/2.0) db.x-=size.x;
							if(db.x<-size.x/2.0) db.x+=size.x;
							if(db.y>size.y/2.0) db.y-=size.y;
							if(db.y<-size.y/2.0) db.y+=size.y;
							if(db.z>size.z/2.0) db.z-=size.z;
							if(db.z<-size.z/2.0) db.z+=size.z;
							
							//More of the same definitions of the bonded forces.
							System.harmonicF(da, a[first], a[second], cBond);
							System.harmonicF(db, a[second], a[third], cBond);
							System.bendF(da, db, a[first], a[second], a[third], cBend);
						}
					}
					
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
		}
		
		//Simplified version of Verlet::second()
		for(int k=0;k<nParticles;k++)
		{
			v[k].x+=(a[k].x*halfDt);
			v[k].y+=(a[k].y*halfDt);
			v[k].z+=(a[k].z*halfDt);
		}
		
		//Measurements are output here
		if(i%measureint==0 && i!=startInt)
		{
			time_t last=current;
			current=time(NULL);
			
			//time since last storage step, good for benchmarking
			std::cout << System.readInitialTime() << '\t' << current-last << std::endl;
			
			current=time(NULL);
		}
		
		//someMeasureInterval would be an integer, like every 10 steps rather than 10 tau
		//if(i%someMeasureInterval && i!=0)
		//{
		//	//you would put the measure here
		//}
		//you could potentially add more measures here
	}
	
	return 0;
}

