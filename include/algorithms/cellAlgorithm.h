//algorithm for calculating forces and potentials between neighbors
#include "functions.h"
//#ifdef _OPENMP
//#include <omp.h>
//#endif

#ifndef MD_CELL
#define MD_CELL

//This algorithm is capable of doing particles per thread
template <typename T>
class Cell {
	public:
		//regular "all particles"
		Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size);
		//Only indexed particles
		Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, int *index, int nI);
		~Cell();
		
		
		//Query functions:
		int query(int particle, threeVector<T> &minImg);
		int query(int particle);
		
		int queryHalf(int particle, threeVector<T> &minImg);
		int queryHalf(int particle);
		
		int query(position<T> pos, threeVector<T> &minImg);
		int query(position<T> pos);
		
		int queryRange(int particle, T radius);
		int queryRange(int particle, threeVector<T> &minImg, T radius);
		
		int queryRange(position<T> pos, T radius);
		int queryRange(position<T> pos, threeVector<T> &minImg, T radius);
		
		//get the cell index
		threeVector<int> getCell(position<T> pos);//based on position
		threeVector<int> getCell(int particle);//based on index
		
		//This one goes to this cell, returning the head particle
		//Just put in the approximate index, it will wrap the index for you!
		int gotoCell(threeVector<int> cellPos);
		
		//Retrieves next element, -1 when finished with a particular cell
		int nextElement();
		
		//resets currentParticle state for query functions per thread
		void resetState();
		
		//resize functions:
		void resize(threeVector<T> size);
		void changeN(int nI, int *index);
		//T computeCustom();
		//the custom function, don't use this with a compute intensive section, like accelerations
		//virtual T function();
		
		//builds linked lists, do this before calling neighbors
		void build();
		
		//return some properties of the system
		fourVector<int> readNCells()
		{
			return nCells;
		};
		
		threeVector<T> readCellSize()
		{
			return cellSize;
		};
	private:
		threeVector<int> nextBox[13];
		
		position<T> *p;//particle positions
		int nP;//nParticles
		T cutoff;//maximum length of interactions
		threeVector<T> s;//size
		
		//two ways to make this work,
		//create indices then fill Cells while updating linked list normally
		//or create indices then fill Cells while sorting linked list as a key-value list, like a hash
		int *cells;
		bool *flag;//solvent flag
		int *linkedList;
		int *indices;
		
		//#ifdef _OPENMP
		//	threeVector<int> *neighbor;
		//	int *currentCell;//keeps track of current cell being called
		//	int *currentParticle;//keeps track of current particle in linked list
		//	fourVector<int> *inquiryCell;//keeps track of box in which inquiry is queried
		//	int *inquiryParticle;//keeps track of particle in which range is queried
		//	int *nextBoxIndex;
		//#else
			//same as above, but whithout threads
			threeVector<int> neighbor;
			int currentCell;
			int currentParticle;
			position<T> cPos;
			fourVector<int> inquiryCell;
			int inquiryParticle;
			int nextBoxIndex;
		//#endif
		
		fourVector<int> nCells;
		int nCellsOld;
		threeVector<T> cellSize;
		
		//flagged with null and 0
		int *index;
		int nI;
		bool reset;
};

template <typename T>
Cell<T>::Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size)
{
	this->p=particles;
	if(p==NULL)
	{
		std::cout << "p is null in Cell!\n";
		throw 0;
	}
	this->nP=nParticles;
	this->cutoff=cutoff;
	this->s=size;
	
	if(this->cutoff<=0)
	{
		std::cout << "No cutoff in Cell!\n";
		throw 0;
	}
	this->nCells.x=int(this->s.x/this->cutoff);
	this->nCells.y=int(this->s.y/this->cutoff);
	this->nCells.z=int(this->s.z/this->cutoff);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->nCellsOld=nCells.t*2;
	this->cells=new int[nCells.t*2];
	this->flag=new bool[nCells.t*2];
	this->linkedList=new int[nP];
	this->indices=new int[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	//parallel or not, whatever
	//#ifdef _OPENMP
	//	neighbor=new threeVector<int>[omp_get_max_threads()];
	//	currentCell=new int[omp_get_max_threads()];
	//	currentParticle=new int[omp_get_max_threads()];
	//	inquiryCell=new fourVector<int>[omp_get_max_threads()];
	//	inquiryParticle=new int[omp_get_max_threads()];
	//	nextBoxIndex=new int[omp_get_max_threads()];
	//	
	//	#pragma parallel
	//	{
	//		this->resetState();
	//	}
	//#else
		//same as above, but whithout threads
		this->resetState();
	//#endif
	
	this->index=NULL;
	this->nI=0;
	
	nextBox[0].x=-1;
	nextBox[0].y=0;
	nextBox[0].z=0;
	nextBox[1].x=0;
	nextBox[1].y=-1;
	nextBox[1].z=0;
	nextBox[2].x=-1;
	nextBox[2].y=-1;
	nextBox[2].z=0;
	nextBox[3].x=0;
	nextBox[3].y=0;
	nextBox[3].z=-1;
	nextBox[4].x=-1;
	nextBox[4].y=0;
	nextBox[4].z=-1;
	nextBox[5].x=0;
	nextBox[5].y=-1;
	nextBox[5].z=-1;
	nextBox[6].x=-1;
	nextBox[6].y=-1;
	nextBox[6].z=-1;
	nextBox[7].x=1;
	nextBox[7].y=0;
	nextBox[7].z=-1;
	nextBox[8].x=0;
	nextBox[8].y=1;
	nextBox[8].z=-1;
	nextBox[9].x=1;
	nextBox[9].y=-1;
	nextBox[9].z=-1;
	nextBox[10].x=-1;
	nextBox[10].y=1;
	nextBox[10].z=-1;
	nextBox[11].x=1;
	nextBox[11].y=1;
	nextBox[11].z=-1;
	nextBox[12].x=-1;
	nextBox[12].y=1;
	nextBox[12].z=0;
}

template <typename T>
Cell<T>::Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, int *index, int nI)
{
	this->p=particles;
	if(p==NULL)
	{
		std::cout << "p is null in Cell!\n";
		throw 0;
	}
	this->nP=nParticles;
	this->cutoff=cutoff;
	this->s=size;
	
	if(this->cutoff<=0)
	{
		std::cout << "No cutoff in Cell!\n";
		throw 0;
	}
	this->nCells.x=int(this->s.x/this->cutoff);
	this->nCells.y=int(this->s.y/this->cutoff);
	this->nCells.z=int(this->s.z/this->cutoff);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->nCellsOld=nCells.t*2;
	this->cells=new int[nCells.t*2];
	this->flag=new bool[nCells.t*2];
	this->linkedList=new int[nI];
	this->indices=new int[nI];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	//parallel or not, whatever
	//#ifdef _OPENMP
	//	neighbor=new threeVector<int>[omp_get_max_threads()];
	//	currentCell=new int[omp_get_max_threads()];
	//	currentParticle=new int[omp_get_max_threads()];
	//	inquiryCell=new fourVector<int>[omp_get_max_threads()];
	//	inquiryParticle=new int[omp_get_max_threads()];
	//	nextBoxIndex=new int[omp_get_max_threads()];
	//	
	//	#pragma parallel
	//	{
	//		this->resetState();
	//	}
	//#else
		//same as above, but whithout threads
		this->resetState();
	//#endif
	
	this->index=index;
	this->nI=nI;
	
	nextBox[0].x=-1;
	nextBox[0].y=0;
	nextBox[0].z=0;
	nextBox[1].x=0;
	nextBox[1].y=-1;
	nextBox[1].z=0;
	nextBox[2].x=-1;
	nextBox[2].y=-1;
	nextBox[2].z=0;
	nextBox[3].x=0;
	nextBox[3].y=0;
	nextBox[3].z=-1;
	nextBox[4].x=-1;
	nextBox[4].y=0;
	nextBox[4].z=-1;
	nextBox[5].x=0;
	nextBox[5].y=-1;
	nextBox[5].z=-1;
	nextBox[6].x=-1;
	nextBox[6].y=-1;
	nextBox[6].z=-1;
	nextBox[7].x=1;
	nextBox[7].y=0;
	nextBox[7].z=-1;
	nextBox[8].x=0;
	nextBox[8].y=1;
	nextBox[8].z=-1;
	nextBox[9].x=1;
	nextBox[9].y=-1;
	nextBox[9].z=-1;
	nextBox[10].x=-1;
	nextBox[10].y=1;
	nextBox[10].z=-1;
	nextBox[11].x=1;
	nextBox[11].y=1;
	nextBox[11].z=-1;
	nextBox[12].x=-1;
	nextBox[12].y=1;
	nextBox[12].z=0;
}

template <typename T>
Cell<T>::~Cell()
{
	delete this->cells;
	delete this->linkedList;
	delete this->indices;
	delete this->flag;
	//#ifdef _OPENMP
	//	delete this->currentCell;
	//	delete this->currentParticle;
	//	delete this->inquiryCell;
	//	delete this->inquiryParticle;
	//	delete this->neighbor;
	//	delete this->nextBoxIndex;
	//#endif
}

//single threaded only!
template <typename T>
void Cell<T>::build()
{
	//non-sorting version
	int buf,current,x,y,z;
	
	//empties the boxes
	//#pragma omp parallel for
	for(int i=0;i<nCells.t;i++)
	{
		cells[i]=-1;
		//this section doesn't matter until other sections are fixed
		//flag[i]=false;
	}
	
	if(nI==0 && index==NULL)
	{
		//calculate each particle's cell hash
		//#pragma omp parallel for private(x,y,z)
		for(int i=0;i<nP;i++)
		{
			x=int(p[i].x/cellSize.x);
			y=int(p[i].y/cellSize.y);
			z=int(p[i].z/cellSize.z);
			indices[i]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[i]=i;//index of the particle
		}
		
		//fill cells
		for(int i=0;i<nP;i++)
		{
			//push operation
			current=indices[i];
			//This section is broken
			#ifdef SOLVENT_FLAG
				if(p[i].type!=SOLVENT_FLAG)
				{
					flag[current]=true;
					//current box
					x=current%nCells.x;
					y=int(current/nCells.x)%nCells.y;
					z=int(current/(nCells.x*nCells.y));
					for(int k=-1;k<2;k++)
					{
						for(int l=-1;l<2;l++)
						{
							for(int m=-1;m<2;m++)
							{
								int nextX=x+k;
								int nextY=y+l;
								int nextZ=z+m;
					
								int adjacent=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
								flag[adjacent]=true;
							}
						}
					}
				}
			#endif
			//This section doesn't matter
			//#ifndef SOLVENT_FLAG
			//	flag[current]=true;
			//#endif
			buf=cells[current];
			cells[current]=linkedList[i];
			linkedList[i]=buf;
		}
	}
	else
	{
		//calculate each particle's cell hash
		//#pragma omp parallel for private(x,y,z)
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			x=int(p[j].x/cellSize.x);
			y=int(p[j].y/cellSize.y);
			z=int(p[j].z/cellSize.z);
			indices[j]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[i]=i;//index of the particle
		}
		//fill cells
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			//push operation
			current=indices[j];
			//this sections is broken
			#ifdef SOLVENT_FLAG
				if(p[j].type!=SOLVENT_FLAG)
				{
					flag[current]=true;
					//current box
					x=current%nCells.x;
					y=int(current/nCells.x)%nCells.y;
					z=int(current/(nCells.x*nCells.y));
					for(int k=-1;k<2;k++)
					{
						for(int l=-1;l<2;l++)
						{
							for(int m=-1;m<2;m++)
							{
								int nextX=x+k;
								int nextY=y+l;
								int nextZ=z+m;
						
								int adjacent=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
								flag[adjacent]=true;
							}
						}
					}
				}
			#endif
			//this section doesn't matter
			//#ifndef SOLVENT_FLAG
			//	flag[current]=true;
			//#endif
			buf=cells[current];
			cells[current]=linkedList[i];
			linkedList[i]=buf;
		}
	}
	//parallel or not, whatever
	//#ifdef _OPENMP
	//	#pragma parallel
	//	{
	//		this->resetState();
	//	}
	//#else
		//same as above, but whithout threads
		this->resetState();
	//#endif
}

//this can be executed in parallel or not
template <typename T>
void Cell<T>::resetState()
{
	
	//#ifdef _OPENMP
	//	inquiryParticle[omp_get_thread_num()]=-1;
	//#else
		//same as above, but whithout threads
		inquiryParticle=-1;
		cPos.x=-1;
		cPos.y=-1;
		cPos.z=-1;
		reset=true;
	//#endif
}

//this can be executed in parallel or not
template <typename T>
inline int Cell<T>::query(int particle, threeVector<T> &minImg)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
	/*#ifdef _OPENMP
		int i=omp_get_thread_num();
		//reset current cell parameters
		if(particle!=inquiryParticle[i])
		{
			int i=omp_get_thread_num();
			//set inquiryParticle
			inquiryParticle[i]=particle;
			
			//set inquiryCell
			inquiryCell[i].x=int(p[inquiryParticle[i]].x/cellSize.x);
			inquiryCell[i].y=int(p[inquiryParticle[i]].y/cellSize.y);
			inquiryCell[i].z=int(p[inquiryParticle[i]].z/cellSize.z);
			inquiryCell[i].t=inquiryCell[i].x+inquiryCell[i].y*nCells.x+inquiryCell[i].z*nCells.x*nCells.y;
			
			//set neighbor
			neighbor[i].x=-2;
			neighbor[i].y=-1;
			neighbor[i].z=-1;
			
			currentParticle[i]=-1;
		}
		else
		{
			//get next particle in linked list
			currentParticle[i]=linkedList[currentParticle[i]];
		}
		
		//check if it is last particle in cell
		while(currentParticle[i]==-1)
		{
			//go to next cell, this increments through each axis
			neighbor[i].x++;
			if(neighbor[i].x==2)//cache misses are going to be bad with this
			{
				neighbor[i].x=-1;
				neighbor[i].y++;
				if(neighbor[i].y==2)
				{
					neighbor[i].y=-1;
					neighbor[i].z++;
					if(neighbor[i].z==2)
					{
						//last neighbor
						currentParticle[i]=-1;
						return currentParticle[i];
					}
				}
			}
			
			//get current cell
			int x=inquiryCell[i].x+neighbor[i].x;
			int y=inquiryCell[i].y+neighbor[i].y;
			int z=inquiryCell[i].z+neighbor[i].z;
			
			//minimum image, basically returns a size to increment current particle by
			minImg.x=0;
			if(x<0) minImg.x=s.x;
			if(x>=nCells.x) minImg.x=-s.x;
			
			minImg.y=0;
			if(y<0) minImg.y=s.y;
			if(y>=nCells.y) minImg.y=-s.y;
			
			minImg.z=0;
			if(z<0) minImg.z=s.z;
			if(z>=nCells.z) minImg.z=-s.z;
			
			//set currentCell
			currentCell[i]=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle[i]=cells[currentCell[i]];
		}
			
		return currentParticle[i];
	#else
	*/
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int(p[inquiryParticle].x/cellSize.x);
			inquiryCell.y=int(p[inquiryParticle].y/cellSize.y);
			inquiryCell.z=int(p[inquiryParticle].z/cellSize.z);
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			//set neighbor
			neighbor.x=-2;
			neighbor.y=-1;
			neighbor.z=-1;
			
			currentParticle=-1;
		}
		else
		{
			//get next particle in linked list
			currentParticle=linkedList[currentParticle];
		}
		
		//check if it is last particle in cell
		while(currentParticle==-1)
		{
			//go to next cell, this increments through each axis
			neighbor.x++;
			if(neighbor.x==2)//cache misses are going to be bad with this
			{
				neighbor.x=-1;
				neighbor.y++;
				if(neighbor.y==2)
				{
					neighbor.y=-1;
					neighbor.z++;
					if(neighbor.z==2)
					{
						neighbor.z=-1;
						//last neighbor
						currentParticle=-1;
						return currentParticle;
					}
				}
			}
			
			//get current cell
			int x=inquiryCell.x+neighbor.x;
			int y=inquiryCell.y+neighbor.y;
			int z=inquiryCell.z+neighbor.z;
			
			//minimum image, basically returns a size to increment current particle by
			minImg.x=0;
			if(x<0) minImg.x=s.x;
			if(x>=nCells.x) minImg.x=-s.x;
			
			minImg.y=0;
			if(y<0) minImg.y=s.y;
			if(y>=nCells.y) minImg.y=-s.y;
			
			minImg.z=0;
			if(z<0) minImg.z=s.z;
			if(z>=nCells.z) minImg.z=-s.z;
			
			//set currentCell
			currentCell=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle=cells[currentCell];
		}
			
		return currentParticle;
	//#endif
}

//this can be executed in parallel or not, this doesn't use minimum image,
// but it does return particles in cells periodically, just remember to ignore them
template <typename T>
inline int Cell<T>::query(int particle)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
	/*#ifdef _OPENMP
		int i=omp_get_thread_num();
		//reset current cell parameters
		if(particle!=inquiryParticle[i])
		{
			
			//set inquiryParticle
			inquiryParticle[i]=particle;
			
			//set inquiryCell
			inquiryCell[i].x=int(p[inquiryParticle[i]].x/cellSize.x);
			inquiryCell[i].y=int(p[inquiryParticle[i]].y/cellSize.y);
			inquiryCell[i].z=int(p[inquiryParticle[i]].z/cellSize.z);
			inquiryCell[i].t=inquiryCell[i].x+inquiryCell[i].y*nCells.x+inquiryCell[i].z*nCells.x*nCells.y;
			
			//set neighbor
			neighbor[i].x=-2;
			neighbor[i].y=-1;
			neighbor[i].z=-1;
			
			currentParticle[i]=-1;
		}
		else
		{
			//get next particle in linked list
			currentParticle[i]=linkedList[currentParticle[i]];
		}
		
		//check if it is last particle in list
		while(currentParticle[i]==-1)
		{
			//go to next cell, this increments through each axis
			neighbor[i].x++;
			if(neighbor[i].x==2)//cache misses are going to be bad with this
			{
				neighbor[i].x=-1;
				neighbor[i].y++;
				if(neighbor[i].y==2)
				{
					neighbor[i].y=-1;
					neighbor[i].z++;
					if(neighbor[i].z==2)
					{
						//last neighbor
						currentParticle[i]=-1;
						return currentParticle[i];
					}
				}
			}
			
			//get current cell
			int x=inquiryCell[i].x+neighbor[i].x;
			int y=inquiryCell[i].y+neighbor[i].y;
			int z=inquiryCell[i].z+neighbor[i].z;
			
			//set currentCell
			currentCell[i]=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle[i]=cells[currentCell[i]];
		}
		
		return currentParticle[i];
	#else
	*/
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int(p[inquiryParticle].x/cellSize.x);
			inquiryCell.y=int(p[inquiryParticle].y/cellSize.y);
			inquiryCell.z=int(p[inquiryParticle].z/cellSize.z);
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			//set neighbor
			neighbor.x=-2;
			neighbor.y=-1;
			neighbor.z=-1;
			
			currentParticle=-1;
		}
		else
		{
			//get next particle in linked list
			currentParticle=linkedList[currentParticle];
		}
		
		//check if it is last particle in list
		while(currentParticle==-1)
		{
			//go to next cell, this increments through each axis
			neighbor.x++;
			if(neighbor.x==2)//cache misses are going to be bad with this
			{
				neighbor.x=-1;
				neighbor.y++;
				if(neighbor.y==2)
				{
					neighbor.y=-1;
					neighbor.z++;
					if(neighbor.z==2)
					{
						neighbor.z=-1;
						//last neighbor
						currentParticle=-1;
						return currentParticle;
					}
				}
			}
			
			//get current cell
			int x=inquiryCell.x+neighbor.x;
			int y=inquiryCell.y+neighbor.y;
			int z=inquiryCell.z+neighbor.z;
			
			//set currentCell
			currentCell=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle=cells[currentCell];
		}
		
		return currentParticle;
	//#endif
}

//this can be executed in parallel or not
template <typename T>
inline int Cell<T>::queryHalf(int particle, threeVector<T> &minImg)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
	/*#ifdef _OPENMP
		int i=omp_get_thread_num();
		
		//reset current cell parameters
		if(particle!=inquiryParticle[i])
		{
			//set inquiryParticle
			inquiryParticle[i]=particle;
			
			//set inquiryCell
			inquiryCell[i].x=int(p[inquiryParticle[i]].x/cellSize.x);
			inquiryCell[i].y=int(p[inquiryParticle[i]].y/cellSize.y);
			inquiryCell[i].z=int(p[inquiryParticle[i]].z/cellSize.z);
			inquiryCell[i].t=inquiryCell[i].x+inquiryCell[i].y*nCells.x+inquiryCell[i].z*nCells.x*nCells.y;
			
			//set neighbor, gotcha, actually just 0
			neighbor[i].x=0;
			neighbor[i].y=0;
			neighbor[i].z=0;
			
			nextBoxIndex[i]=-1;
			
			currentCell[i]=inquiryCell[i].t;
			
			//minimum image, basically returns a size to increment current particle by
			minImg.x=0;
			minImg.y=0;
			minImg.z=0;
			
			//get current particle, which is one after current particle, unless it is the last one in the list
			currentParticle[i]=inquiryParticle[i];
		}
		
		//branching is bad, m'kay
		//get next particle in linked list
		currentParticle[i]=linkedList[currentParticle[i]];
		
		//check if it is last particle in list
		while(currentParticle[i]==-1 && nextBoxIndex[i]<13)
		{
			//grab next set of coordinates
			nextBoxIndex[i]++;
			neighbor[i]=nextBox[nextBoxIndex[i]];
			
			//get current cell
			int x=inquiryCell[i].x+neighbor[i].x;
			int y=inquiryCell[i].y+neighbor[i].y;
			int z=inquiryCell[i].z+neighbor[i].z;
			
			//minimum image, basically returns a size to increment current particle by
			minImg.x=0;
			if(x<0) minImg.x=s.x;
			if(x>=nCells.x) minImg.x=-s.x;
			
			minImg.y=0;
			if(y<0) minImg.y=s.y;
			if(y>=nCells.y) minImg.y=-s.y;
			
			minImg.z=0;
			if(z<0) minImg.z=s.z;
			if(z>=nCells.z) minImg.z=-s.z;
			
			//set currentCell
			currentCell[i]=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle[i]=cells[currentCell[i]];
		}
		
		return currentParticle[i];
	#else
	*/
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int(p[inquiryParticle].x/cellSize.x);
			inquiryCell.y=int(p[inquiryParticle].y/cellSize.y);
			inquiryCell.z=int(p[inquiryParticle].z/cellSize.z);
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			//set neighbor, gotcha, actually just 0
			neighbor.x=0;
			neighbor.y=0;
			neighbor.z=0;
			
			nextBoxIndex=-1;
			
			currentCell=inquiryCell.t;
			
			//minimum image, basically returns a size to increment current particle by
			minImg.x=0;
			minImg.y=0;
			minImg.z=0;
			
			//get current particle, which is one after current particle, unless it is the last one in the list
			currentParticle=inquiryParticle;
		}
		
		//branching is bad, m'kay
		//get next particle in linked list
		currentParticle=linkedList[currentParticle];
		
		//check if it is last particle in list
		while(currentParticle==-1 && nextBoxIndex<13)
		{
			//grab next set of coordinates
			nextBoxIndex++;
			neighbor=nextBox[nextBoxIndex];
			
			//get current cell
			int x=inquiryCell.x+neighbor.x;
			int y=inquiryCell.y+neighbor.y;
			int z=inquiryCell.z+neighbor.z;
			
			//minimum image, basically returns a size to increment current particle by
			minImg.x=0;
			if(x<0) minImg.x=s.x;
			if(x>=nCells.x) minImg.x=-s.x;
			
			minImg.y=0;
			if(y<0) minImg.y=s.y;
			if(y>=nCells.y) minImg.y=-s.y;
			
			minImg.z=0;
			if(z<0) minImg.z=s.z;
			if(z>=nCells.z) minImg.z=-s.z;
			
			//set currentCell
			currentCell=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle=cells[currentCell];
		}
		
		return currentParticle;
	//#endif
}

//this can be executed in parallel or not, doesn't use minimum image
template <typename T>
inline int Cell<T>::queryHalf(int particle)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
	/*#ifdef _OPENMP
		int i=omp_get_thread_num();
		
		//reset current cell parameters
		if(particle!=inquiryParticle[i])
		{
			//set inquiryParticle
			inquiryParticle[i]=particle;
			
			//set inquiryCell
			inquiryCell[i].x=int(p[inquiryParticle[i]].x/cellSize.x);
			inquiryCell[i].y=int(p[inquiryParticle[i]].y/cellSize.y);
			inquiryCell[i].z=int(p[inquiryParticle[i]].z/cellSize.z);
			inquiryCell[i].t=inquiryCell[i].x+inquiryCell[i].y*nCells.x+inquiryCell[i].z*nCells.x*nCells.y;
			
			//set neighbor, gotcha, actually just 0
			neighbor[i].x=0;
			neighbor[i].y=0;
			neighbor[i].z=0;
			
			nextBoxIndex[i]=-1;
			
			currentCell[i]=inquiryCell[i].t;
			
			//get current particle, which is one after current particle, unless it is the last one in the list
			currentParticle[i]=inquiryParticle[i];
		}
		
		//branching is bad, m'kay
		//get next particle in linked list
		currentParticle[i]=linkedList[currentParticle[i]];
		
		
		//check if it is last particle in list
		while(currentParticle[i]==-1 && nextBoxIndex[i]<13)
		{
			//grab next set of coordinates
			nextBoxIndex[i]++;
			neighbor[i]=nextBox[nextBoxIndex[i]];
			
			//get current cell
			int x=inquiryCell[i].x+neighbor[i].x;
			int y=inquiryCell[i].y+neighbor[i].y;
			int z=inquiryCell[i].z+neighbor[i].z;
			
			//set currentCell
			currentCell[i]=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle[i]=cells[currentCell[i]];
		}
		
		return currentParticle[i];
	#else
	*/
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int(p[inquiryParticle].x/cellSize.x);
			inquiryCell.y=int(p[inquiryParticle].y/cellSize.y);
			inquiryCell.z=int(p[inquiryParticle].z/cellSize.z);
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			//set neighbor, gotcha, actually just 0
			neighbor.x=0;
			neighbor.y=0;
			neighbor.z=0;
			
			nextBoxIndex=-1;
			
			currentCell=inquiryCell.t;
			
			//get current particle, which is one after current particle, unless it is the last one in the list
			currentParticle=inquiryParticle;
		}
		
		//branching is bad, m'kay
		//get next particle in linked list
		currentParticle=linkedList[currentParticle];
		
		
		//check if it is last particle in list
		while(currentParticle==-1 && nextBoxIndex<13)
		{
			//grab next set of coordinates
			nextBoxIndex++;
			neighbor=nextBox[nextBoxIndex];
			
			//get current cell
			int x=inquiryCell.x+neighbor.x;
			int y=inquiryCell.y+neighbor.y;
			int z=inquiryCell.z+neighbor.z;
			
			//set currentCell
			currentCell=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle=cells[currentCell];
		}
		
		return currentParticle;
	//#endif
}

//this can be executed in parallel or not
template <typename T>
inline int Cell<T>::query(position<T> pos, threeVector<T> &minImg)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
		//reset current cell parameters
		if(reset)
		{
			//set inquiryCell
			inquiryCell.x=int(pos.x/cellSize.x);
			inquiryCell.y=int(pos.y/cellSize.y);
			inquiryCell.z=int(pos.z/cellSize.z);
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			if(inquiryCell.x>nCells.x || inquiryCell.x<0 ||
			inquiryCell.y>nCells.y || inquiryCell.y<0 ||
			inquiryCell.z>nCells.z || inquiryCell.z<0)
			{
				std::cout << "Error [cellAlgorithm]:Particle in query is out of range!\n";
				throw 0;
			}
			
			//set neighbor
			neighbor.x=-2;
			neighbor.y=-1;
			neighbor.z=-1;
			
			currentParticle=-1;
			reset=false;
		}
		else
		{
			//get next particle in linked list
			currentParticle=linkedList[currentParticle];
		}
		
		//check if it is last particle in cell
		while(currentParticle==-1 && !reset)
		{
			//go to next cell, this increments through each axis
			neighbor.x++;
			if(neighbor.x==2)
			{
				neighbor.x=-1;
				neighbor.y++;
			}
			
			if(neighbor.y==2)
			{
				neighbor.y=-1;
				neighbor.z++;
			}
			
			if(neighbor.z==2)
			{
				neighbor.z=-1;
				reset=true;
			}
			
			if(!reset)
			{
				//get current cell
				int x=inquiryCell.x+neighbor.x;
				int y=inquiryCell.y+neighbor.y;
				int z=inquiryCell.z+neighbor.z;
				
				//minimum image, basically returns a size to increment current particle by
				minImg.x=0;
				if(x<0) minImg.x=s.x;
				if(x>=nCells.x) minImg.x=-s.x;
				
				minImg.y=0;
				if(y<0) minImg.y=s.y;
				if(y>=nCells.y) minImg.y=-s.y;
				
				minImg.z=0;
				if(z<0) minImg.z=s.z;
				if(z>=nCells.z) minImg.z=-s.z;
				
				//set currentCell
				currentCell=((x+nCells.x)%nCells.x)+
					(((y+nCells.y)%nCells.y)*nCells.x)+
					(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
				//get current particle
				currentParticle=cells[currentCell];
			}
		}
		return currentParticle;
}

//this can be executed in parallel or not, this doesn't use minimum image,
// but it does return particles in cells periodically, just remember to ignore them
template <typename T>
inline int Cell<T>::query(position<T> pos)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
		//reset current cell parameters
		if(pos.x != cPos.x || pos.y != cPos.y || pos.z != cPos.z)
		{
			//set cPos
			cPos=pos;
			
			//set inquiryCell
			inquiryCell.x=int(cPos.x/cellSize.x);
			inquiryCell.y=int(cPos.y/cellSize.y);
			inquiryCell.z=int(cPos.z/cellSize.z);
			
			if(inquiryCell.x>nCells.x || inquiryCell.x<0 ||
			inquiryCell.y>nCells.y || inquiryCell.y<0 ||
			inquiryCell.z>nCells.z || inquiryCell.z<0)
			{
				std::cout << "Error [cellAlgorithm]:Particle in query is out of range!\n";
				throw 0;
			}
			
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			//set neighbor
			neighbor.x=-2;
			neighbor.y=-1;
			neighbor.z=-1;
			
			currentParticle=-1;
		}
		else
		{
			//get next particle in linked list
			currentParticle=linkedList[currentParticle];
		}
		
		//check if it is last particle in list
		while(currentParticle==-1)
		{
			//go to next cell, this increments through each axis
			neighbor.x++;
			if(neighbor.x==2)//cache misses are going to be bad with this
			{
				neighbor.x=-1;
				neighbor.y++;
				if(neighbor.y==2)
				{
					neighbor.y=-1;
					neighbor.z++;
					if(neighbor.z==2)
					{
						//last neighbor
						currentParticle=-1;
						return currentParticle;
					}
				}
			}
			
			//get current cell
			int x=inquiryCell.x+neighbor.x;
			int y=inquiryCell.y+neighbor.y;
			int z=inquiryCell.z+neighbor.z;
			
			//set currentCell
			currentCell=((x+nCells.x)%nCells.x)+
				(((y+nCells.y)%nCells.y)*nCells.x)+
				(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
			//get current particle
			currentParticle=cells[currentCell];
		}
		
		return currentParticle;
	//#endif
}

template <typename T>
inline int Cell<T>::queryRange(int particle, T radius)
{
	//reset current cell parameters
	if(particle!=inquiryParticle)
	{
		//set inquiryParticle
		inquiryParticle=particle;
		
		//set inquiryCell
		inquiryCell.x=int(p[inquiryParticle].x/cellSize.x);
		inquiryCell.y=int(p[inquiryParticle].y/cellSize.y);
		inquiryCell.z=int(p[inquiryParticle].z/cellSize.z);
		inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
		
		//defines our range
		threeVector<int> range;
		range.x=floor(radius/cellSize.x)+1;
		range.y=floor(radius/cellSize.y)+1;
		range.z=floor(radius/cellSize.z)+1;
		
		//set neighbor
		neighbor.x=-range.x-1;
		neighbor.y=-range.y;
		neighbor.z=-range.z;
		
		currentParticle=-1;
	}
	else
	{
		//get next particle in linked list
		currentParticle=linkedList[currentParticle];
	}
	
	//check if it is last particle in list
	while(currentParticle==-1)
	{
		//defines our range
		threeVector<int> range;
		range.x=floor(radius/cellSize.x)+1;
		range.y=floor(radius/cellSize.y)+1;
		range.z=floor(radius/cellSize.z)+1;
		
		//go to next cell, this increments through each axis
		neighbor.x++;
		if(neighbor.x==range.x)//cache misses are going to be bad with this
		{
			neighbor.x=-range.x;
			neighbor.y++;
			if(neighbor.y==range.y)
			{
				neighbor.y=-range.y;
				neighbor.z++;
				if(neighbor.z==range.z)
				{
					neighbor.z=-range.z;
					//last neighbor
					currentParticle=-1;
					return currentParticle;
				}
			}
		}
		
		//get current cell
		int x=inquiryCell.x+neighbor.x;
		int y=inquiryCell.y+neighbor.y;
		int z=inquiryCell.z+neighbor.z;
		
		//set currentCell
		currentCell=((x+nCells.x)%nCells.x)+
		(((y+nCells.y)%nCells.y)*nCells.x)+
		(((z+nCells.z)%nCells.z)*nCells.x*nCells.y);
		
		//get current particle
		currentParticle=cells[currentCell];
	}
	
	return currentParticle;
}

inline int Cell<T>::queryRange(int particle, threeVector<T> &minImg, T radius)
{
	
}

template <typename T>
inline int Cell<T>::queryRange(position<T> pos, T radius)
{
	
}

inline int Cell<T>::queryRange(position<T> pos, threeVector<T> &minImg, T radius)
{
	
}

template <typename T>
int Cell<T>::gotoCell(threeVector<int> cellPos)
{
	//For now we will just assume no one will put anything much greater than nCells
	inquiryCell.x=(cellPos.x+nCells.x)%nCells.x;
	inquiryCell.y=(cellPos.y+nCells.y)%nCells.y;
	inquiryCell.z=(cellPos.z+nCells.z)%nCells.z;
	
	//set currentCell
	currentCell=((inquiryCell.x+nCells.x)%nCells.x)+
		(((inquiryCell.y+nCells.y)%nCells.y)*nCells.x)+
		(((inquiryCell.z+nCells.z)%nCells.z)*nCells.x*nCells.y);
	
	//get current particle
	currentParticle=cells[currentCell];
	
	return currentParticle;
}

template <typename T>
int Cell<T>::nextElement()
{
	if(currentParticle!=-1)
		currentParticle=linkedList[currentParticle];
	return currentParticle;
}

template <typename T>
threeVector<int> Cell<T>::getCell(position<T> pos)
{
	threeVector<int> cellIndex;
	cellIndex.x=int(pos.x/cellSize.x);
	cellIndex.y=int(pos.y/cellSize.y);
	cellIndex.z=int(pos.z/cellSize.z);
	return cellIndex;
}

template <typename T>
threeVector<int> Cell<T>::getCell(int particle)
{
	threeVector<int> cellIndex;
	cellIndex.x=int(p[particle].x/cellSize.x);
	cellIndex.y=int(p[particle].y/cellSize.y);
	cellIndex.z=int(p[particle].z/cellSize.z);
	return cellIndex;
}

//This changes the number of particles by changing the redirection indices
template <typename T>
void Cell<T>::changeN(int nI, int *index)
{
	this->nI=nI;
	this->index=index;
	
	delete linkedList,indices;
	
	linkedList=new int[nI];
	indices=new int[nI];
}

template <typename T>
void Cell<T>::resize(threeVector<T> size)
{
	//rescale number of cells
	nCells.x=size.x/cutoff;
	nCells.y=size.y/cutoff;
	nCells.z=size.z/cutoff;
	nCells.t=nCells.x*nCells.y*nCells.z;
	
	//check if new rescaling is higher than allocated number of boxes, this usually isn't the case
	if(nCellsOld<nCells.t)
	{
		delete cells;
		cells=new int[nCells.t*2];
		nCellsOld=nCells.t*2;
	}
	
	//rescale cell size
	cellSize.x=size.x/nCells.x;
	cellSize.y=size.y/nCells.y;
	cellSize.z=size.z/nCells.z;
	
	//rescale overall size
	s=size;
}


//end of MD_CELL
#endif
