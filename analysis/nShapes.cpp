//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include <ctime>
#include <vector>


class dist {
	public:
		dist(position<double> *pos, threeVector<double> size)
		{
			this->p=pos;
			this->s=size;
		};
		
		void setCenter(threeVector<double> i)
		{
			center=i;
		};
		
		void setShell(fourVector<double> i, double thickness)
		{
			shell=i;
			inner=shell.t*shell.t;
			outer=shell.t+thickness;
			outer*=outer;
		};
		
		//binary operator
		bool operator () (int i, int j)
		{
			threeVector<double> di,dj;
			di.x=p[i].x-center.x;
			di.y=p[i].y-center.y;
			di.z=p[i].z-center.z;
			
			//Minimum image for i, this makes things a little wierd...
			//things get further away until s/2.0, then they get closer...
			di.x+=(di.x>s.x/2.0)?-s.x:0;
			di.y+=(di.y>s.y/2.0)?-s.y:0;
			di.z+=(di.z>s.z/2.0)?-s.z:0;
			
			di.x+=(di.x<-s.x/2.0)?s.x:0;
			di.y+=(di.y<-s.y/2.0)?s.y:0;
			di.z+=(di.z<-s.z/2.0)?s.z:0;
			
			dj.x=p[j].x-center.x;
			dj.y=p[j].y-center.y;
			dj.z=p[j].z-center.z;
			
			//Minimum image for j
			dj.x+=(dj.x>s.x/2.0)?-s.x:0;
			dj.y+=(dj.y>s.y/2.0)?-s.y:0;
			dj.z+=(dj.z>s.z/2.0)?-s.z:0;
			
			dj.x+=(dj.x<-s.x/2.0)?s.x:0;
			dj.y+=(dj.y<-s.y/2.0)?s.y:0;
			dj.z+=(dj.z<-s.z/2.0)?s.z:0;
			
			//if the distance's are equal, then it will be at the end of the list (same particle)
			if((di.x*di.x+di.y*di.y+di.z*di.z)<(dj.x*dj.x+dj.y*dj.y+dj.z*dj.z))
				return true;
			///I'm not sure which is more efficient... Compiler probably knows!
			//else
				return false;
		};
		
		//unary operator
		bool operator () (int i)
		{
			threeVector<double> d;
			d.x=p[i].x-shell.x;
			d.y=p[i].y-shell.y;
			d.z=p[i].z-shell.z;
			
			//Minimum image for i, this makes things a little wierd...
			//things get further away until s/2.0, then they get closer...
			d.x+=(d.x>s.x/2.0)?-s.x:0;
			d.y+=(d.y>s.y/2.0)?-s.y:0;
			d.z+=(d.z>s.z/2.0)?-s.z:0;
			
			d.x+=(d.x<-s.x/2.0)?s.x:0;
			d.y+=(d.y<-s.y/2.0)?s.y:0;
			d.z+=(d.z<-s.z/2.0)?s.z:0;
			
			double dr=d.x*d.x+d.y*d.y+d.z*d.z;
			if(dr<outer)
				return false;
			//else
				return true;
		}
		
	private:
		position<double> *p;
		threeVector<double> s;
		threeVector<double> center;
		fourVector<double> shell;
		double inner,outer;
};

int main(int argc, char* argv[])
{
#ifdef MD_SYSTEM
	if(argc!=4)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name cutoff cytoType\n";
		return 0;
	}
#else
	if(argc!=8)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << " name cutoff cytoType anchorDistance sizeX sizeY sizeZ\n";
		return 0;
	}
#endif

	///read in options
	char *name=argv[1];
	double cutoff=atof(argv[2]);
	int cytoType=atoi(argv[3]);
	double anchorDistance=atof(argv[4]);
	
	//registration thresholds, each one defines a new behavior for that threshold
	//int nThreshold=argc-5;
	//double *regThreshold=new double[nThreshold];
	//for(int i=0;i<nThreshold;i++)
	//	regThreshold[i]=atof(argv[i+5]);
#ifdef MD_SYSTEM
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//general initializations
	MTRand randNum(System.readSeed());
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	threeVector<double> s=System.readSize();
#else
	MTRand randNum(time(NULL));
	position<double> *p=NULL;
	int nParticles=0;
	threeVector<double> s;
	s.x=atof(argv[5]);
	s.y=atof(argv[6]);
	s.z=atof(argv[7]);
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	xyzFile.load();
#endif
	
	int *stack=new int[nParticles];
	int *index=new int[nParticles];//Sortable indices
	bool *exclude=new bool[nParticles];
	
	int nIndex=nParticles;
	//std::cout << nIndex << '\n';
	Cell<double> neighbors(p, nParticles, cutoff, s, index, nIndex);
	Cell<double> largeSweep(p, nParticles, anchorDistance, s);
	
	
	if(xyzFile.loadLast())//load one frame
	//while(xyzFile.load())//load all frames
	{
		//exclude membrane near cytoskeleton
		for(int i=0;i<nParticles;i++)
		{
			exclude[i]=false;
			///if index is used in largeSweep, this isn't needed!
			index[i]=i;
		}
		
		largeSweep.build();
		
		double anchorDistanceSqr=anchorDistance*anchorDistance;
		for(int i=0;i<nParticles;i++)
		{
			if(p[i].type==cytoType)
			{
				largeSweep.resetState();
				exclude[i]=true;
				threeVector<double> minImg(0);
				for(int j=largeSweep.query(i,minImg);j!=-1;j=largeSweep.query(i,minImg))
				{
					if(i!=j)
					{
						threeVector<double> d;
						d.x=p[i].x-p[j].x+minImg.x;
						d.y=p[i].y-p[j].y+minImg.y;
						d.z=p[i].z-p[j].z+minImg.z;
						
						double dr=d.x*d.x+d.y*d.y+d.z*d.z;
						
						//is it within the shell?
						if(dr<anchorDistanceSqr)
							exclude[j]=true;
					}
				}
			}
		}
		
		//swap exclusions to the end
		int pivot=nParticles-1;
		for(int i=0;i<pivot;i++)
		{
			if(exclude[i])
			{
				for(;exclude[pivot];pivot--);
				if(pivot>=0)
				{
					bool bBuf;
					int iBuf;
					iBuf=index[i];
					index[i]=index[pivot];
					index[pivot]=iBuf;
					
					bBuf=exclude[i];
					exclude[i]=exclude[pivot];
					exclude[pivot]=bBuf;
				}
			}
		}
		
		if(pivot>0)
		{
			neighbors.changeN(pivot,index);
			
			neighbors.build();
			
			for(int i=0;i<pivot;i++)
				exclude[i]=false;
			
			std::vector< std::vector <int> > blebs;
			
			double cutoffSqr=cutoff*cutoff;
			///This might be the best stack implementation I've made. All others need to use this.
			for(int i=0;i<pivot;i++)
			{
				//if it hasn't been excluded
				if(!exclude[i])
				{
					std::vector<int> bleb;
					//push it on stack and begin neighbor search
					stack[0]=i;
					exclude[i]=true;
					for(int j=0;j>-1;)
					{
						//pop it off stack
						int current=stack[j--];
						
						bleb.push_back(current);
						
						neighbors.resetState();
						//std::cout << current << '\t' << pivot << '\n';
						
						//find nearest neighbors
						threeVector<double> minImg(0);
						for(int k=neighbors.queryIndex(current,minImg);k!=-1;k=neighbors.queryIndex(current,minImg))
						{
							//only check range if it isn't excluded
							if(!exclude[k])
							{
								int iIndex=index[current];
								int kIndex=index[k];
								threeVector<double> d;
								d.x=p[iIndex].x-p[kIndex].x+minImg.x;
								d.y=p[iIndex].y-p[kIndex].y+minImg.y;
								d.z=p[iIndex].z-p[kIndex].z+minImg.z;
								
								double dr=d.x*d.x+d.y*d.y+d.z*d.z;
								
								//is it within the range?
								if(dr<cutoffSqr)
								{
									//exclude it and push it on the stack
									exclude[k]=true;
									stack[++j]=k;
								}
							}
						}
					}
					blebs.push_back(bleb);
				}
			}
			
			std::cout << name << '\t' << blebs.size() << '\n';
		}
		else
		{
			//no blebs found
			std::cout << name << '\t' << 0 << '\n';
		}
		//std::cout << "Number of Blebs: " << blebs.size() << std::endl;
		//for(int i=0;i<blebs.size();i++)
		//	std::cout << "Size of bleb " << i << ": " << blebs[i].size() << std::endl;
		
		//std::cout << pivot << "\nList\n";
		
		//for(int i=0;i<blebs.size();i++)
		//{
		//	for(int j=0;j<blebs[i].size();j++)
		//	{
		//		int k=index[blebs[i][j]];
		//		std::cout << i << '\t' << p[k].x << '\t' << p[k].y << '\t' << p[k].z << '\n';
		//	}
		//}
		
	}
	
	
	/*
	///Determine shape filling
	//while(xyzFile.load())
	{
		std::vector< fourVector<double> > shell;
		std::vector<bool> cytoContent;
		bool cytoInside=false;
		
		//center of mass object
		Com<double> com(p, nParticles, s);
		
		//object to compare distances
		dist comCompare(p,System.readSize());
				
		//set center at center of mass
		threeVector<double> centerOfMass=com.compute();
		comCompare.setCenter(centerOfMass);
		
		//initialize indices
		//for(int i=0;i<nParticles;i++)
		//	index[i]=i;
		
		//sort by distance from center of mass
		std::sort(&index[0],&index[nParticles],comCompare);
		
		int nBlebs=0;
		time_t begin_time=time(NULL);
		//endPart is the last particle in the list, go until all particles are used
		for(int endPart=nParticles;endPart>-1;)
		{
			
			fourVector<double> currentShell;
			//Initial conditions
			currentShell.x=p[index[endPart-1]].x;
			currentShell.y=p[index[endPart-1]].y;
			currentShell.z=p[index[endPart-1]].z;
			currentShell.t=0;//Inner shell radius
			
			int nOld=0;//number of particles before decision is made
			
			///This method isn't that bad
			///But revise the line search loop to do this once, not nDegrees*nParticles
			//Initializes nOld to the smallest shell around particle
			#pragma omp parallel for reduction(+:nOld)
			for(int i=0;i<endPart;i++)
			{
				threeVector<double> d;
				d.x=p[index[i]].x-currentShell.x;
				d.x+=((d.x>s.x/2.0)?-s.x:0);
				d.x+=((d.x<-s.x/2.0)?s.x:0);
				
				d.y=p[index[i]].y-currentShell.y;
				d.y+=((d.y>s.y/2.0)?-s.y:0);
				d.y+=((d.y<-s.y/2.0)?s.y:0);
				
				d.z=p[index[i]].z-currentShell.z;
				d.z+=((d.z>s.z/2.0)?-s.z:0);
				d.z+=((d.z<-s.z/2.0)?s.z:0);
				double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
				
				if(dr<thickness)
					nOld++;
			}
			//std::cout << nOld << '\n';
			//std::cin.get();
			int noChange=0;
			
			//Begin line search, a decision per degree of freedom is made in this loop.
			//Go until number of no changes (decisions) is made within the noChangeRatio.
			//Direction is determined from gradient ascent.
			//The function itself is unknown, but we know the output is based on N particles.
			for(int j=1;double(noChange)/double(j)<noChangeRatio;j++)
			{
				int nNew=0;//number of new particles to decide if it should change
				int nChange=0;//total number of new particles in changes
				
				//try a move in the x direction
				double scale=(2.0*randNum.rand53()-1.0)*dD;
				currentShell.x+=scale;
				///This is a terrible method!
				#pragma omp parallel for reduction(+:nNew)
				for(int i=0;i<endPart;i++)
				{
					threeVector<double> d;
					d.x=p[index[i]].x-currentShell.x;
					d.x+=((d.x>s.x/2.0)?-s.x:0);
					d.x+=((d.x<-s.x/2.0)?s.x:0);
					
					d.y=p[index[i]].y-currentShell.y;
					d.y+=((d.y>s.y/2.0)?-s.y:0);
					d.y+=((d.y<-s.y/2.0)?s.y:0);
					
					d.z=p[index[i]].z-currentShell.z;
					d.z+=((d.z>s.z/2.0)?-s.z:0);
					d.z+=((d.z<-s.z/2.0)?s.z:0);
					double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
					
					//is it within the shell?
					if(dr>currentShell.t && dr<currentShell.t+thickness)
						nNew++;
				}
				if(nNew<=nOld)
				{
					currentShell.x-=scale;
				}
				else
				{
					nChange+=(nNew-nOld);
					nOld=nNew;
				}
				
				//try a move in the y direction
				scale=(2.0*randNum.rand53()-1.0)*dD;
				currentShell.y+=scale;
				nNew=0;
				///This is a terrible method!
				#pragma omp parallel for reduction(+:nNew)
				for(int i=0;i<endPart;i++)
				{
					threeVector<double> d;
					d.x=p[index[i]].x-currentShell.x;
					d.x+=((d.x>s.x/2.0)?-s.x:0);
					d.x+=((d.x<-s.x/2.0)?s.x:0);
					
					d.y=p[index[i]].y-currentShell.y;
					d.y+=((d.y>s.y/2.0)?-s.y:0);
					d.y+=((d.y<-s.y/2.0)?s.y:0);
					
					d.z=p[index[i]].z-currentShell.z;
					d.z+=((d.z>s.z/2.0)?-s.z:0);
					d.z+=((d.z<-s.z/2.0)?s.z:0);
					double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
					
					//is it within the shell?
					if(dr>currentShell.t && dr<currentShell.t+thickness)
						nNew++;
				}
				if(nNew<=nOld)
				{
					currentShell.y-=scale;
				}
				else
				{
					nChange+=(nNew-nOld);
					nOld=nNew;
				}
				
				//try a move in the z direction
				scale=(2.0*randNum.rand53()-1.0)*dD;
				currentShell.z+=scale;
				nNew=0;
				
				///This is a terrible method!
				#pragma omp parallel for reduction(+:nNew)
				for(int i=0;i<endPart;i++)
				{
					threeVector<double> d;
					d.x=p[index[i]].x-currentShell.x;
					d.x+=((d.x>s.x/2.0)?-s.x:0);
					d.x+=((d.x<-s.x/2.0)?s.x:0);
					
					d.y=p[index[i]].y-currentShell.y;
					d.y+=((d.y>s.y/2.0)?-s.y:0);
					d.y+=((d.y<-s.y/2.0)?s.y:0);
					
					d.z=p[index[i]].z-currentShell.z;
					d.z+=((d.z>s.z/2.0)?-s.z:0);
					d.z+=((d.z<-s.z/2.0)?s.z:0);
					double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
					
					//is it within the shell?
					if(dr>currentShell.t && dr<currentShell.t+thickness)
						nNew++;
				}
				if(nNew<=nOld)
				{
					currentShell.z-=scale;
				}
				else
				{
					nChange+=(nNew-nOld);
					nOld=nNew;
				}
				
				//try a move in the r direction
				///r cannot be less than 0 (negative radius)
				for(scale=(2.0*randNum.rand53()-1.0)*dD;scale+currentShell.t<0;)
					scale=(2.0*randNum.rand53()-1.0)*dD;
				currentShell.t+=scale;
				nNew=0;
				
				bool cytoContained=false;
				
				
				///This is a terrible method!
				#pragma omp parallel for reduction(+:nNew)
				for(int i=0;i<endPart;i++)
				{
					cytoContained=(p[index[i]].type==CYTO);
					threeVector<double> d;
					d.x=p[index[i]].x-currentShell.x;
					d.x+=((d.x>s.x/2.0)?-s.x:0);
					d.x+=((d.x<-s.x/2.0)?s.x:0);
					
					d.y=p[index[i]].y-currentShell.y;
					d.y+=((d.y>s.y/2.0)?-s.y:0);
					d.y+=((d.y<-s.y/2.0)?s.y:0);
					
					d.z=p[index[i]].z-currentShell.z;
					d.z+=((d.z>s.z/2.0)?-s.z:0);
					d.z+=((d.z<-s.z/2.0)?s.z:0);
					double dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
					
					//is it within the shell?
					if(dr>currentShell.t && dr<currentShell.t+thickness)
						nNew++;
				}
				if(nNew<=nOld)
				{
					currentShell.t-=scale;
				}
				else
				{
					nChange+=(nNew-nOld);
					nOld=nNew;
				}
					cytoInside=cytoContained;
				
				//part of our stop conditions
				if(nChange==0)
					noChange++;
			}
			
			//Now that a shell has been found, pivot indices to the end
			comCompare.setShell(currentShell,thickness);
			
			///Stable partition from c++ lib maintains order and reshuffles
			endPart=std::stable_partition(&index[0],&index[endPart],comCompare)-index-1;
			
			//std::cout << nOld << '\n';
			shell.push_back(currentShell);//go to next shell
			cytoContent.push_back(cytoInside);
		}
		
		std::vector< std::vector<int> > joinedShells;
		std::vector< std::vector<int> > mergedShells;
		
		//sphere registration (volume of overlap)
		for(int i=0;i<shell.size();i++)
		{
			std::vector<int> merged;
			std::vector<int> joined;
			for(int j=0;j<shell.size() && shell[i].t>0;j++)
			{
				if(i!=j && shell[j].t>0)
				{
					fourVector<double> a=shell[j];
					fourVector<double> b=shell[i];
					fourVector<double> d=shell[i]; 
					d.x-=a.x;
					d.y-=a.y;
					d.z-=a.z;
					d.t=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
					
					double vOverlap=(a.t+b.t-d.t);
					vOverlap*=vOverlap;
					vOverlap*=M_PI*(d.t*d.t+2.0*d.t*a.t-3.0*a.t*a.t+2.0*d.t*b.t+6.0*a.t*b.t-3.0*a.t*a.t);
					vOverlap/=12.0*d.t;
					
					double vCurrent=(4.0/3.0)*M_PI*a.t*a.t*a.t;
					
					///only three conditions should be possible:
					//if it exceeds the maximum spatial registration, it isn't a bleb
					if(vOverlap/vCurrent>regThreshold[2])
					{
						std::cout << "Merged " << j << " into " << i << std::endl;
						merged.push_back(j);
					}
					//if there is registration, record it, it might be a neck
					else if(vOverlap/vCurrent>regThreshold[1])
					{
						std::cout << "Joined " << j << " into " << i << std::endl;
						joined.push_back(j);
					}
					else if(vOverlap/vCurrent>regThreshold[0])
					{
						//nothing for now
					}
				}
			}
			mergedShells.push_back(merged);
			joinedShells.push_back(joined);
			
			if(cytoContent[i])
			{
				std::cout << "Main: " << shell[i].t << std::endl;
			}
			else
			{
				std::cout << "Bleb?: " << shell[i].t << std::endl;
			}
		}
		
		int rootShell=-1;
		
		
		
		time_t end_time=time(NULL);
		
		std::cout << end_time-begin_time << '\n';
		std::cout << shell.size() << '\n';
		std::cout << nBlebs << '\n';
		
		std::cout << 13*shell.size() << "\ntest\n";
		for(int i=0;i<shell.size();i++)
		{
			std::cout << 1 << '\t' << shell[i].x << '\t' << shell[i].y << '\t' << shell[i].z << std::endl;
			std::cout << 2 << '\t' << shell[i].x+shell[i].t << '\t' << shell[i].y << '\t' << shell[i].z << std::endl;
			std::cout << 2 << '\t' << shell[i].x-shell[i].t << '\t' << shell[i].y << '\t' << shell[i].z << std::endl;
			std::cout << 2 << '\t' << shell[i].x << '\t' << shell[i].y+shell[i].t << '\t' << shell[i].z << std::endl;
			std::cout << 2 << '\t' << shell[i].x << '\t' << shell[i].y-shell[i].t << '\t' << shell[i].z << std::endl;
			std::cout << 2 << '\t' << shell[i].x << '\t' << shell[i].y << '\t' << shell[i].z+shell[i].t << std::endl;
			std::cout << 2 << '\t' << shell[i].x << '\t' << shell[i].y << '\t' << shell[i].z-shell[i].t << std::endl;
			std::cout << 3 << '\t' << shell[i].x+shell[i].t+thickness << '\t' << shell[i].y << '\t' << shell[i].z << std::endl;
			std::cout << 3 << '\t' << shell[i].x-shell[i].t-thickness << '\t' << shell[i].y << '\t' << shell[i].z << std::endl;
			std::cout << 3 << '\t' << shell[i].x << '\t' << shell[i].y+shell[i].t+thickness << '\t' << shell[i].z << std::endl;
			std::cout << 3 << '\t' << shell[i].x << '\t' << shell[i].y-shell[i].t-thickness << '\t' << shell[i].z << std::endl;
			std::cout << 3 << '\t' << shell[i].x << '\t' << shell[i].y << '\t' << shell[i].z+shell[i].t+thickness << std::endl;
			std::cout << 3 << '\t' << shell[i].x << '\t' << shell[i].y << '\t' << shell[i].z-shell[i].t-thickness << std::endl;
		}
		
	}
	
	//When sphere doesn't fit, try a cylinder?
	
	//std::string newName("shapes_");
	//newName+=name;
	//newName+=".xyz";
	
	//xyzFormat<double> xyzFile(p,nParticles);
	//xyzFile.open(newName.c_str(), std::ios::out);
	//xyzFile.store();
	//xyzFile.close();
	*/
#ifndef MD_SYSTEM
	xyzFile.close();
#endif
	
	delete index;
	delete exclude;
	delete stack;
	
	return 0;
}

