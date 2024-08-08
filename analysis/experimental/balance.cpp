/** \brief Testing load balance factor.
 *  This tests the load balance factor to see if there is an optimal load blancing technique.
 */

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

class kSplitting {
	public:
		kSplitting(position<double> *particles, unsigned int nParticles)
		{
			p=particles;
			nP=nParticles;
			k=0;
		};
		~kSplitting(){};
		void chooseK(int kDim)
		{
			k=kDim%3;
		};
		bool operator() (unsigned int i, unsigned int j)
		{
			if(i>nParticles || j>nParticles)
			{
				std::cout << "Particle " << i << " or " << j << " out of range [0," << nP-1 << "]!" << std::endl;
				throw 0;
			}
			return p[i].s[k]>p[j].s[k];
		};
	private:
		position<double> *p;
		unsigned int nP;
		int k;
};

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		std::cout << argv[0] << " name \n";
		return 0;
	}
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	double cutoff=System.readCutoff();
	double cutoffSqr=cutoff*cutoff;
	threeVector<double> size=System.readSize();
	unsigned int nParticles=System.readNParticles();
	position<double> *p=System.getPositions();
	
	double maxSize=size.x>size.y?size.x:size.y;
	maxSize=maxSize>size.z?maxSize:size.z;
	//double patchLength=80;
	for(double patchLength=cutoff;patchLength<maxSize/2.0;patchLength+=1.0)
	{
		//using a basic grid for load balancing, shared memory.
		std::map<unsigned int, unsigned int> cells;
		std::vector< keyVal<unsigned int,unsigned int> > hashIndices;
		fourVector<unsigned int> nCells;
		threeVector<double> cellSize;
		
		nCells.x=floor(size.x/patchLength);
		nCells.y=floor(size.y/patchLength);
		nCells.z=floor(size.z/patchLength);
		nCells.t=nCells.x*nCells.y*nCells.z;
		cellSize.x=size.x/nCells.x;
		cellSize.y=size.y/nCells.y;
		cellSize.z=size.z/nCells.z;
		
		for(unsigned int i=0;i<nParticles;i++)
		{
			unsigned int x=int(p[i].x/cellSize.x)%nCells.x;
			unsigned int y=int(p[i].y/cellSize.y)%nCells.y;
			unsigned int z=int(p[i].z/cellSize.z)%nCells.z;
			
			keyVal<unsigned int, unsigned int> hI;
			
			hI.value=i;
			hI.key=x+y*nCells.x+z*nCells.x*nCells.y;
			hashIndices.push_back(hI);
			
			std::map<unsigned int, unsigned int>::iterator it;
			it=cells.find(hashIndices[i].key);
			//cell linked list
			if(it!=cells.end())
			{
				//push operation
				unsigned int buf=cells[hashIndices[i].key];//use our old cell head as the next list tail
				cells[hashIndices[i].key]=hashIndices[i].value;//new cell list head
				hashIndices[i].value=buf;
			}
			else
			{
				//initialize operation
				cells[hashIndices[i].key]=hashIndices[i].value;
				hashIndices[i].value=-1;//it is now the last iterator
			}
		}
//std::cout << nParticles << std::endl << "test" << std::endl;
		//average is easy to calculate
		double avg=static_cast<double>(nParticles)/static_cast<double>(cells.size());
//int cellNumber=0;
		double stdDev=0;
		double avgOverlap=0;
		//calculate standard deviation and average overlap
		for(std::map<unsigned int, unsigned int>::iterator current=cells.begin();
			current!=cells.end();
			++current)
		{
			double nInCell=0;
			double nNearby=0;
			for(unsigned int m=current->second; m!=-1; m=hashIndices[m].value)
			{
				nInCell++;
//std::cout << cellNumber << '\t' << p[m].x << '\t' << p[m].y << '\t' << p[m].z << std::endl; 
			}
//cellNumber++;
			stdDev+=(nInCell-avg)*(nInCell-avg);
			
			//unhash our keys
			int x=current->first%nCells.x;
			int y=int(current->first/nCells.x)%nCells.y;
			int z=int(current->first/(nCells.x*nCells.y));
			
			//look at nearby keys
			for(int a=-1;a<2;a++)
			{
				for(int b=-1;b<2;b++)
				{
					for(int c=-1;c<2;c++)
					{
						//compute a nearby key
						int nextX=(a+x)%nCells.x;
						int nextY=(b+y)%nCells.y;
						int nextZ=(c+z)%nCells.z;
						
						unsigned int nearbyKey=nextX+(nextY*nCells.x)+(nextZ*nCells.x*nCells.y);
						
						if(cells.find(nearbyKey)!=cells.end() && nearbyKey!=current->first)
						{
							threeVector<double> boundary;
						
							boundary.x=(a<0?static_cast<double>(x)*cellSize.x:0);//left bound
							boundary.x=(a>0?static_cast<double>(nextX)*cellSize.x:boundary.x);//right bound
							boundary.y=(b<0?static_cast<double>(y)*cellSize.y:0);//left bound
							boundary.y=(b>0?static_cast<double>(nextY)*cellSize.y:boundary.y);//right bound
							boundary.z=(c<0?static_cast<double>(z)*cellSize.z:0);//left bound
							boundary.z=(c>0?static_cast<double>(nextZ)*cellSize.z:boundary.z);//right bound
							for(int n=cells[nearbyKey]; n!=-1; n=hashIndices[n].value)
							{
								threeVector<double> d;
								d.x=(a!=0?p[n].x-boundary.x:0);
								d.y=(b!=0?p[n].y-boundary.y:0);
								d.z=(c!=0?p[n].z-boundary.z:0);
//std::cout << 1 << '\t' << p[n].x << '\t' << p[n].y << '\t' << p[n].z << std::endl;
								if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
								{
//std::cout << 1 << '\t' << p[n].x << '\t' << p[n].y << '\t' << p[n].z << std::endl;
									nNearby++;
								}
							}
						}
					}
				}
			}
			avgOverlap+=nNearby;
		}
		stdDev/=static_cast<double>(cells.size());
		stdDev=sqrt(stdDev);
		avgOverlap/=static_cast<double>(cells.size());
		
		std::cerr << "Ghost Patch algorithm splitting:"
		<< "\n\tAverage PPC: " << avg 
		<< "\n\tStandard Deviation: " << stdDev 
		<< "\n\tDensity: " << avg/(cellSize.x*cellSize.y*cellSize.z) 
		<< "\n\tAverage Overlap PPC: " << avgOverlap
		<< std::endl;
		std::cout << cells.size() << '\t' << avg << '\t' << stdDev << '\t' << avg/(cellSize.x*cellSize.y*cellSize.z) << '\t' << avgOverlap << std::endl;
	
	}
	
	kSplitting predicate(p,nParticles);
	std::vector<unsigned int> indices;
	for(unsigned int i=0;i<nParticles;i++)
		indices.push_back(0);
	
	for(int nDivisions=1;nDivisions<nParticles/2;nDivisions*=2)
	{
		
		
		for(unsigned int i=0;i<nParticles;i++)
		{
			unsigned int x=int(p[i].x/cellSize.x)%nCells.x;
			unsigned int y=int(p[i].y/cellSize.y)%nCells.y;
			unsigned int z=int(p[i].z/cellSize.z)%nCells.z;
			
			keyVal<unsigned int, unsigned int> hI;
			
			hI.value=i;
			hI.key=x+y*nCells.x+z*nCells.x*nCells.y;
			hashIndices.push_back(hI);
			
			std::map<unsigned int, unsigned int>::iterator it;
			it=cells.find(hashIndices[i].key);
			//cell linked list
			if(it!=cells.end())
			{
				//push operation
				unsigned int buf=cells[hashIndices[i].key];//use our old cell head as the next list tail
				cells[hashIndices[i].key]=hashIndices[i].value;//new cell list head
				hashIndices[i].value=buf;
			}
			else
			{
				//initialize operation
				cells[hashIndices[i].key]=hashIndices[i].value;
				hashIndices[i].value=-1;//it is now the last iterator
			}
		}
		
		//average is easy to calculate
		double avg=static_cast<double>(nParticles)/static_cast<double>(cells.size());
		
		double stdDev=0;
		double avgOverlap=0;
		//calculate standard deviation and average overlap
		for(std::map<unsigned int, unsigned int>::iterator current=cells.begin();
			current!=cells.end();
			++current)
		{
			double nInCell=0;
			double nNearby=0;
			for(unsigned int m=current->second; m!=-1; m=hashIndices[m].value)
				nInCell++;
			stdDev+=(nInCell-avg)*(nInCell-avg);
			
			//unhash our keys
			int x=current->first%nCells.x;
			int y=int(current->first/nCells.x)%nCells.y;
			int z=int(current->first/(nCells.x*nCells.y));
			
			//look at nearby keys
			for(int a=-1;a<2;a++)
			{
				for(int b=-1;b<2;b++)
				{
					for(int c=-1;c<2;c++)
					{
						//compute a nearby key
						int nextX=(a+x)%nCells.x;
						int nextY=(b+y)%nCells.y;
						int nextZ=(c+z)%nCells.z;
						threeVector<double> boundary;
						
						boundary.x=(a<0?x*cellSize.x:0);//left bound
						boundary.x=(a>0?(x+1)*cellSize.x:0);//right bound
						boundary.y=(b<0?y*cellSize.y:0);//left bound
						boundary.y=(b>0?(y+1)*cellSize.y:0);//right bound
						boundary.z=(c<0?z*cellSize.z:0);//left bound
						boundary.z=(c>0?(z+1)*cellSize.z:0);//right bound
						
						unsigned int nearbyKey=nextX+(nextY*nCells.x)+(nextZ*nCells.x*nCells.y);
						
						if(cells.find(nearbyKey)!=cells.end())
						{
							for(int n=cells[nearbyKey]; n!=-1; n=hashIndices[n].value)
							{
								threeVector<double> d;
								d.x=(boundary.x!=0?p[n].x-boundary.x:0);
								d.y=(boundary.y!=0?p[n].y-boundary.y:0);
								d.z=(boundary.z!=0?p[n].z-boundary.z:0);
								if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
									nNearby++;
							}
						}
					}
				}
			}
			avgOverlap+=nNearby;
		}
		stdDev/=static_cast<double>(pow(2,nDivisions));
		stdDev=sqrt(stdDev);
		avgOverlap/=static_cast<double>(pow(2,nDivisions));
		
		std::cerr << "BSP algorithm splitting:"
		<< "\n\tAverage PPC: " << avg 
		<< "\n\tStandard Deviation: " << stdDev 
		<< "\n\tDensity: " << avg/(cellSize.x*cellSize.y*cellSize.z) 
		<< "\n\tAverage Overlap PPC: " << avgOverlap
		<< std::endl;
		std::cout << (cellSize.x*cellSize.y*cellSize.z) << '\t' << avg << '\t' << stdDev << '\t' << avg/(cellSize.x*cellSize.y*cellSize.z) << '\t' << avgOverlap << std::endl;
	
	}
	
	return 0;
}
