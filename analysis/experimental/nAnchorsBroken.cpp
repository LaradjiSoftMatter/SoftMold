#include <iostream>
#include <cstdlib>

#define SOLVENT 0

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

using namespace std;

class dist {
	public:
		dist(position<double> *pos, threeVector<double> size)
		{
			this->p=pos;
			this->s=size;
		};
		
		void setCenter(int i)
		{
			center=i;
		};
		
		bool operator () (int i, int j)
		{
			threeVector<double> di,dj;
			di.x=p[i].x-p[center].x;
			di.y=p[i].y-p[center].y;
			di.z=p[i].z-p[center].z;
			
			//Minimum image for i, this makes things a little wierd...
			//things get further away until s/2.0, then they get closer...
			di.x=di.x>s.x/2.0?di.x-s.x:di.x;
			di.y=di.y>s.y/2.0?di.y-s.y:di.y;
			di.z=di.z>s.z/2.0?di.z-s.z:di.z;
			
			di.x=di.x<-s.x/2.0?di.x+s.x:di.x;
			di.y=di.y<-s.y/2.0?di.y+s.y:di.y;
			di.z=di.z<-s.z/2.0?di.z+s.z:di.z;
			
			dj.x=p[j].x-p[center].x;
			dj.y=p[j].y-p[center].y;
			dj.z=p[j].z-p[center].z;
			
			//Minimum image for j
			dj.x=dj.x>s.x/2.0?dj.x-s.x:dj.x;
			dj.y=dj.y>s.y/2.0?dj.y-s.y:dj.y;
			dj.z=dj.z>s.z/2.0?dj.z-s.z:dj.z;
			
			dj.x=dj.x<-s.x/2.0?dj.x+s.x:dj.x;
			dj.y=dj.y<-s.y/2.0?dj.y+s.y:dj.y;
			dj.z=dj.z<-s.z/2.0?dj.z+s.z:dj.z;
			
			//if the distance's are equal, then it will be at the end of the list (same particle)
			if((di.x*di.x+di.y*di.y+di.z*di.z)<(dj.x*dj.x+dj.y*dj.y+dj.z*dj.z))
				return true;
			return false;
		};
		
	private:
		position<double> *p;
		threeVector<double> s;
		int center;
};

#define MD_SYSTEM

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		//so simple, this can't possibly mess it up
		std::cout << "usage: " << argv[0] << "name anchorType headType strandType minTime(tau)\n";
		return 0;
	}
	
	char *name=argv[1];
	int anchorType=atoi(argv[2]);//usually 5
	int headType=atoi(argv[3]);//usually 2, it is assumed that this type starts the list
	int strandType=atoi(argv[4]);//usually 4
	double minTime=atof(argv[5]);
	
#ifdef MD_SYSTEM
	///the variables for the simulation
	Blob<double> System;
	
	///load variables, then initialize them, Script requires some functions from Blob
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
#endif
	position <double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	xyzFormat<double> xyzFile(p, nParticles);
	
	std::string framesFilename("frames_");
	framesFilename+=name;
	framesFilename+=".xyz";
	
	xyzFile.open(framesFilename,std::ios::in);
	
	Cell<double> neighbors(System.getPositions(), System.readNParticles(), System.readCutoff(), System.readSize());
	
	std::vector<int> anchor;//anchor indices
	std::vector<int> strand;//strand indices
	int nBlebs=0;//number of blebs
	std::vector<int> whichBleb;//which bleb an anchor resides
	int *bleb=new int[nParticles];//bleb indices
	std::vector<int> dAnchor;//anchor distances, this is updated everytime a broken anchor is found
	std::vector<bool> broken;//number of broken anchors
	std::vector<double> area;
	int headIndex=-1;
	int strandIndex=-1;
	int headMolecule=-1;
	int strandMolecule=-1;
	int anchorIndex=-1;
	int anchorMolecule=-1;
	//these are just shorter names
	threeVector<double> s=System.readSize();
	molecule<double, fourVector<int> > *m=System.getMolecule();
	
	//dAnchors 
	
#ifndef MD_SYSTEM
	while(xyzFile.load())
#else
	if(System.readInitialTime()>=minTime)
#endif
	{
		
		if(anchor.size()==0)
		{
			//where are the starting indices for these molecules?
			for(int i=0;i<nParticles;i++)
			{
				if(p[i].type==headType && headIndex==-1)
					headIndex=i;
				if(p[i].type==strandType && strandIndex==-1)
					strandIndex=i;
				if(p[i].type==anchorType && anchorIndex==-1)
					anchorIndex=i;
			}
			
			//which molecules are they?
#ifdef MD_SYSTEM
			for(int i=0;i<System.readNMolecules();i++)
			{
				if(m[i].readType()==CHAIN)
				{
					if(m[i].getBonds()[0].s[START]==headIndex)
						headMolecule=i;
					if(m[i].getBonds()[0].s[START]==strandIndex)
						strandMolecule=i;
				}
				if(m[i].readType()==BOND)
				{
					for(int j=0;j<m[i].readNBond();j++)
					{
						if(m[i].getBonds()[j].s[0]==anchorIndex ||
							m[i].getBonds()[j].s[1]==anchorIndex)
							anchorMolecule=i;
					}	
				}
			}
#endif			
			for(int i=0;i<nParticles;i++)
				if(p[i].type==anchorType)
					anchor.push_back(i);
				//else if(p[i].type==strandType)
				//	strand.push_back(i);
		}
#ifdef MD_SYSTEM
		//average length of strands
		double strandLength=0;
		for(int i=0;i<m[strandMolecule].getBonds()[0].s[NCHAINS];i++)
		{
			int first=m[strandMolecule].getBonds()[0].s[START]+i*m[strandMolecule].getBonds()[0].s[CHAINLENGTH];
			int last=m[strandMolecule].getBonds()[0].s[START]+(i+1)*m[strandMolecule].getBonds()[0].s[CHAINLENGTH]-1;
			//now that we have the ends, find the anchors and use those endpoints
			bool firstFound=false;
			bool lastFound=false;
			for(int j=0;j<m[anchorMolecule].readNBond();j++)
			{
				int a=m[anchorMolecule].getBonds()[j].s[0];
				int b=m[anchorMolecule].getBonds()[j].s[1];
				
				//the anchor is the one that isn't the strand, first is already the strand
				if(a==first && !firstFound)
				{
					first=b;
					firstFound=true;
				}
				else if(b==first && !firstFound)
				{
					first=a;
					firstFound=true;
				}
				
				//the anchor is the one that isn't the strand, first is already the strand
				if(a==last && !lastFound)
				{
					last=b;
					lastFound=true;
				}
				else if(b==last && !lastFound)
				{
					last=a;
					lastFound=true;
				}
			}
			
			double dx=p[first].x-p[last].x;
			double dy=p[first].y-p[last].y;
			double dz=p[first].z-p[last].z;
			dx-=(dx>s.x/2.0)?s.x:0;
			dy-=(dy>s.y/2.0)?s.y:0;
			dz-=(dz>s.z/2.0)?s.z:0;
			dx+=(dx<-s.x/2.0)?s.x:0;
			dy+=(dy<-s.y/2.0)?s.y:0;
			dz+=(dz<-s.z/2.0)?s.z:0;
			strandLength+=sqrt(dx*dx+dy*dy+dz*dz);
		}
		
		//this is anchor to anchor length or strand length, although it includes the anchored ends
		strandLength/=double(m[strandMolecule].getBonds()[0].s[NCHAINS]);
		
#endif
		
		broken.assign(anchor.size(), true);
		whichBleb.assign(anchor.size(), -1);
		dAnchor.assign(anchor.size(), 0);
		
		int nBroken=anchor.size();
		neighbors.build();
		
		
		double cutoffSqr=System.readCutoff()*System.readCutoff();
		for(int i=0;i<anchor.size();i++)
		{
			threeVector<double> minImg;
			int iIndex=anchor[i];
			//check all neighbors
			for(int j=neighbors.query(iIndex, minImg);j!=-1;j=neighbors.query(iIndex,minImg))
			{
				//of HEAD type
				if(p[j].type==headType)
				{
					threeVector<double> d;
					d.x=p[iIndex].x-p[j].x+minImg.x;
					d.y=p[iIndex].y-p[j].y+minImg.y;
					d.z=p[iIndex].z-p[j].z+minImg.z;
					
					//if close enough to interact, then it isn't broken
					if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					{
						broken[i]=false;
						nBroken--;
						break;//should break the current loop
					}
				}
			}
		}
		
		std::cout << name << '\t' << nBroken << std::endl;
#ifdef MD_SYSTEM
		//this is a cutoff distance squared, and it should be a little longer than the average strand length
		double strandLengthSqr=strandLength*strandLength*2.25;
		double brokenStrandLength=0;
		area.clear();
		int nBrokenLink=0;
		//Approximation assuming each broken anchor contributes 1/3*6=2 cytoskeleton triangles
		// and edges contribute another 2/3 of a triangle
		for(int i=0;i<anchor.size();i++)
		{
			int iIndex=anchor[i];
			if(broken[i])
			{
				//int j=anchor[i];
				//std::cout << "7\t" << p[j].x << '\t' << p[j].y << '\t' << p[j].z << '\n';
				for(int j=0;j<anchor.size();j++)
				{
					if(!broken[j])
					{
						int jIndex=anchor[j];
						double dx=p[iIndex].x-p[jIndex].x;
						double dy=p[iIndex].y-p[jIndex].y;
						double dz=p[iIndex].z-p[jIndex].z;
						dx-=(dx>s.x/2.0)?s.x:0;
						dy-=(dy>s.y/2.0)?s.y:0;
						dz-=(dz>s.z/2.0)?s.z:0;
						dx+=(dx<-s.x/2.0)?s.x:0;
						dy+=(dy<-s.y/2.0)?s.y:0;
						dz+=(dz<-s.z/2.0)?s.z:0;
						if(dx*dx+dy*dy+dz*dz<strandLengthSqr)
						{
							//edges only, 2/3 of triangle
							area+=sqrt(0.75)/6.0;
						}
					}
				}
				//for everybody, 2 triangle contribution
				area+=sqrt(0.75);
			}
		}
		
		area*=strandLength*strandLength;
		
		std::cout << area << std::endl;
#endif
		/*
		for(int i=0;i<anchor.size();i++)
		{
			//if broken or (not broken and not part of a bleb)
			if(broken[i] || (!broken[i] && bleb[i]!=-1))
			{
				//setup of a list of anchors
				for(int j=0;j<anchor.size();j++)
					dAnchor[j]=anchor[j];
				
				//object to compare distances
				dist objectDist(p,s);
				
				//set center
				objectDist.setCenter(i);
				
				//sort by distance from axis
				sort(&dAnchor[0],&dAnchor[anchor.size()-1],objectDist);
				
				//if broken and not part of a bleb
				if(broken[i] && bleb[i]==-1)
				{
					double sum=0;
					int noNeighbor=0;
					for(int j=0;j<anchor.size() || noNeighbor==0;j++)
					{
						threeVector<double> d;
						
						d.x=p[j].x-p[i].x;
						d.y=p[j].y-p[i].y;
						d.z=p[j].z-p[i].z;
						
						d.x=d.x>s.x/2.0?d.x-s.x:d.x;
						d.y=d.y>s.y/2.0?d.y-s.y:d.y;
						d.z=d.z>s.z/2.0?d.z-s.z:d.z;
						
						d.x=d.x<-s.x/2.0?d.x+s.x:d.x;
						d.y=d.y<-s.y/2.0?d.y+s.y:d.y;
						d.z=d.z<-s.z/2.0?d.z+s.z:d.z;
						
						double dr=(d.x*d.x+d.y*d.y+d.z*d.z);
						
						sum+=dr;
						
						double average=sum/(double)i;
						
						if(dr>average*1.5)
							noNeighbor=1;
					}
					
					if(noNeighbor==1)
					{
						
					}
				}
				
				//if broken and part of a bleb (checking for merging)
				if(broken[i] && bleb[i]!=-1)
				{
					
				}
				
				//if not broken, but was part of a bleb
				if(!broken[i] && bleb[i]!=-1)
				{
					
				}
			}
			
		}
		*/
	}
	
	if(bleb!=NULL)
		delete bleb;
	
	return 0;
}


