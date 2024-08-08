//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#define STAT_OUT

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <vector>

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name cutoff\n";
		std::cout << "Output is volume of each inner space consisting of only solvent.\n";
		std::cout << "Extra output is a volume.xyz file. The file shows the volume number\n";
		std::cout << "as a type of particle on a grid for the initial space.\n";
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	double cutoff=atof(argv[2]);
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//some shortcuts
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	threeVector<double> s=System.readSize();
	
	//faster neighbor searching, I will have to double check that cutoff
	Cell<double> neighbors(p, nParticles, cutoff+1.0, s);
	
	///Map and count the clusters
	//create a stack
	int *stack=new int[nParticles];
	//for flagging already used particles
	bool *flag=new bool[nParticles];
	
	//set all flags to false
	for(int i=0;i<nParticles;i++)
		flag[i]=false;
	
	//keep track of clusters
	std::vector<std::vector<int> > cluster;
	
	//build lists
	neighbors.build();
	
	//go through all particles
	for(int i=0;i<nParticles;i++)
	{
		//if current particle isn't flagged
		if(!flag[i] && p[i].type==CYTO)
		{
			std::vector<int> currentCluster;
			currentCluster.push_back(i);
			//recursively add particles to stack
			stack[0]=i;
			flag[i]=true;
			for(int j=0;j>-1;j--)
			{
				threeVector<double> minImg;
				
				//check all neighbors
				///The commented out loop here is really slow
				//go through all particles after current one
				for(int k=neighbors.query(stack[j], minImg);k!=-1;k=neighbors.query(stack[j],minImg))
				{
					//if type is nanoparticles and it hasn't been flagged
					if(p[k].type==CYTO && !flag[k] && k>i)
					{
						//determine distance from k particle
						//minimum image is taken care of with minImg
						threeVector<double> d;
						d.x=p[k].x-p[stack[j]].x+minImg.x;
						d.y=p[k].y-p[stack[j]].y+minImg.y;
						d.z=p[k].z-p[stack[j]].z+minImg.z;
						
						double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
						if(dr<cutoff)
						{
							flag[k]=true;
							stack[j++]=k;
							//add particle to cluster
							currentCluster.push_back(k);
						}
					}
				}
			}
			//pushes back current cluster
			cluster.push_back(currentCluster);
		}
	}
	
	//number of total clusters
	std::cout << "NClusters: " << cluster.size() << '\n';
	//size of cluster "a", you can loop through integer "a" "nClusters" times
	for(int i=0;i<cluster.size();i++)
		std::cout << "Size of cluster " << i << ": " << cluster[i].size() << '\n';
	
	delete stack,flag;
	
	return 0;
}