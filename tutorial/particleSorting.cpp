#include "../include/algorithms/cell.h"//for Cell, position, threeVector
#include <time.h>//for time
#include <stdlib.h>//for rand, RAND_MAX
#include <iostream>//for cout

int main()
{
	//initialize the important stuff
	int nParticles=100000;
	position<double> *p=new position<double>[nParticles];
	threeVector<double> size;
	size.x=50;
	size.y=50;
	size.z=50;
	double cutoff=2.0;
	
	//position particles randomly
	for(int i=0;i<nParticles;i++)
	{
		p[i].x=((double)rand()/(double)RAND_MAX)*size.x;
		p[i].y=((double)rand()/(double)RAND_MAX)*size.y;
		p[i].z=((double)rand()/(double)RAND_MAX)*size.z;
	}
	
	
	///This is the slow (all pairs) method:
	
	//set nNeighbors to 0
	int nNeighbors=0;
	
	//start the clock for all pairs method
	time_t start=time(NULL);
	
	//locate neighbors for particle i
	for(int i=0;i<nParticles;i++)
	{
		//check against all particles
		for(int j=0;j<nParticles;j++)
		{
			//Make sure you don't include yourself
			if(i!=j)
			{
				//compute distance between particle i and j
				double dx=p[i].x-p[j].x;
				double dy=p[i].y-p[j].y;
				double dz=p[i].z-p[j].z;
				double dr=dx*dx+dy*dy+dz*dz;
				
				//increment nNeighbors if within cutoff^2
				if(dr<cutoff*cutoff)
					nNeighbors++;
			}
		}
	}
	//clock stops when all neighbors are found
	time_t end=time(NULL);
	
	//output results using all pair method
	std::cout << "Time taken to locate neighbors using all pair method: " << end-start << " seconds.\n";
	std::cout << "Found " << ((double)nNeighbors/(double)nParticles) << " neighbors per particle.\n";
	
	
	
	
	
	
	
	
	
	
	
	
	
	///This is the fast method:
	
	
	//create class that contains the particle lists, this is needed
	Cell<double> neighbors(p,nParticles,cutoff,size);
	
	//reset nNeighbors to 0 to start with neighbor method
	nNeighbors=0;
	
	//start the clock for neighbor method
	start=time(NULL);
	
	//Build linked lists, this must be done before query, this is needed
	neighbors.build();
	
	//locate neighbors for particle i
	for(int i=0;i<nParticles;i++)
	{
		//query j from neighbor lists based on proximity to particle i, this is needed
		for(int j=neighbors.query(i);j!=-1;j=neighbors.query(i))
		{
			//make sure you don't include yourself
			if(i!=j)
			{
				//compute distance between particle i and j
				double dx=p[i].x-p[j].x;
				double dy=p[i].y-p[j].y;
				double dz=p[i].z-p[j].z;
				double dr=dx*dx+dy*dy+dz*dz;
				
				//increment nNeighbors if it is within cutoff^2
				if(dr<cutoff*cutoff)
					nNeighbors++;
			}
		}
	}
	
	//clock stops when all neighbors are found
	end=time(NULL);
	
	//output results
	std::cout << "Time taken to locate neighbors using Cell method: " << end-start << " seconds.\n";
	std::cout << "Found " << ((double)nNeighbors/(double)nParticles) << " neighbors per particle.\n";
	
	delete p;
	return 0;
}
