//This simply runs a system through the timesteps. It can use OpenMP.
//There was an issue with setting a pointer to the Blob type and allocation, so
// Blob type contains it's own pointers now. Everything is being rewritten for this.

#include "../include/systemMD.h"
#include <ctime>

int main(int argc, char* argv[])
{
	if(argc!=6)
	{
		std::cout << "Usage: " << argv[0] << " N X Y Z CUTOFF\n";
		std::cout << "Benchmark program to show the usefulness of this library.\n";
		std::cout << "Look at benchmark.cpp to see some common optimizations.\n";
		std::cout << "Currently uses OpenMP.\n";
		return 0;
	}

	//Our System object, see include/systemMD.h for more information.
	Blob<double> System;

	//preallocates the number of particles we need
	int nParticles=atoi(argv[1]);
	System.allocParticle(nParticles);
	
	//set the size of our system
	threeVector<double> size;
	size.x=atof(argv[2]);
	size.y=atof(argv[3]);
	size.z=atof(argv[4]);
	System.setSize(size);

	//set cutoff
	double cutoff=atof(argv[5]);
	System.setCutoff(cutoff);

	//Isn't this a nice way to add particles!
	for(int k=0;k<nParticles;k++)
	{
		//normally you would set each of these to some value
		position <double> aPosition;
		threeVector <double> aVelocity;
		threeVector <double> anAcceleration;

		//Now we add the particle information!
		System.addParticle(aPosition,aVelocity,anAcceleration);
	}

	//Random number generator, mersenne twister, just using current time for seed
	MTRand randNum(time(NULL));



	///Setting some values
	//Terrible method, using the object members to individually access data.
	double start=omp_get_wtime();
	for(int k=0;k<System.readNParticles();k++)
	{
		//Each time we access this member (getAccelerations), there is overhead!
		System.getAccelerations()[k].x=0;
		System.getAccelerations()[k].y=0;
		System.getAccelerations()[k].z=0;
		//Same for positions, we are using readsize as well.
		System.getPositions()[k].x=randNum.rand53()*System.readSize().x;
		System.getPositions()[k].y=randNum.rand53()*System.readSize().y;
		System.getPositions()[k].z=randNum.rand53()*System.readSize().z;
	}
	double end=omp_get_wtime();	
	std::cout << "Time to set accelerations using terrible call: " << end-start << std::endl;
	


	//This uses a pointer to the what getAccelerations returns.
	position<double> *p=System.getPositions();
	threeVector<double> *a=System.getAccelerations();
	start=omp_get_wtime();
	for(int k=0;k<nParticles;k++)
	{
		//No overhead, direct access to data!
		a[k].x=0;
		a[k].y=0;
		a[k].z=0;
		//Same for positions, we are using readsize as well.
		p[k].x=randNum.rand53()*size.x;
		p[k].y=randNum.rand53()*size.y;
		p[k].z=randNum.rand53()*size.z;
	}
	end=omp_get_wtime();
	std::cout << "Time to set accelerations using simple call: " << end-start << std::endl;



	///Local pair access
	//Typical method, "all pairs"
	
	//Eliminating the need to multiply cutoff:
	double cutoffSqr=cutoff*cutoff;

	start=omp_get_wtime();
	#pragma omp parallel for 
	for(int k=0;k<nParticles;k++)
	{
		for(int j=0;j<nParticles;j++)
		{
			if(k!=j)
			{
				threeVector<double> d;
				d.x=p[k].x-p[j].x;
				d.y=p[k].y-p[j].y;
				d.z=p[k].z-p[j].z;
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
				{
					//we are just counting neighbors
					a[k].x+=1.0;
					a[k].y+=1.0;
					a[k].z+=1.0;
				}
			}
		}
	}
	end=omp_get_wtime();
	std::cout << "Time to count neighbors using \"all pairs\" method: " << end-start << std::endl;

	//resetting our counts
	for(int i=0;i<nParticles;i++)
	{
		a[i].x=0;
		a[i].y=0;
		a[i].z=0;
	}

	//optimized "cell list" method
	



	return 0;
}
