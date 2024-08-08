/**
 * \brief An object that uses the cell algorithm to find nearest neighbors.
 * The algorithm contained in this is easy to use compared to cellOpt.h or
 * cellOptWV.h. Just create the object with the positions you want, and
 * then run the build member function, which builds the head and linked
 * list. The typical format to run this is similar to an all pairs like
 * algorithm.
 */

//algorithm for calculating forces and potentials between neighbors
#include "functions.h"
//#ifdef _OPENMP
//#include <omp.h>
//#endif

#ifndef MD_CELL
#define MD_CELL

#ifndef MAX_CELL_SIZE
#define MAX_CELL_SIZE 512
#endif

/**
 * \brief This is the cell algorithm object
 * There are two main constructors. One is a simple position index type,
 * and the other is a mapped index version. This algorithm is not thread safe.
 * There are a variety of query members. Each one is detailed within the class.
 */
template <typename T>
class Cell {
	public:
		/**
		 * Default constructor, initializes everything to 0 or NULL.
		 */
		Cell()
		{
			this->initialize();
		};
		void initialize()
		{
			//Stuff that is usually input
			this->p=NULL;
			this->nP=0;
			this->cutoff=0;
			this->s=0;
			this->lOffset=0;
			this->pIndex=NULL;
			this->nI=0;
			
			//Stuff created/allocated by this object
			this->nCells.x=0;
			this->nCells.y=0;
			this->nCells.z=0;
			this->nCells.t=0;
			
			this->nCellsAllocated=0;
			this->cellSize=0;
			
			this->cells=NULL;
			this->linkedList=NULL;
			this->indices=NULL;
			this->cellsIndex=NULL;
			this->linkedListIndex=NULL;
			this->flag=NULL;
			
			#ifdef LOW_DENSITY
				this->fullCells=NULL;
				this->nFullCells=0;
			#endif
			
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
			//std::cout << "I have been Called!\n";
			//std::cin.get();
		};
		
		/**
		 * Regular constructor. Just takes a list of positions. The cutoff and size
		 * are required for it to make nearest neighbor lists. Actually calls initializor 
		 * with same arguments.
		 */
		Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size);
		
		/**
		 * Constructor with mapped indices. Just like the regular constructor, but
		 * it also takes some mapped indices as arguments. Actually calls initializor with 
		 * same arguments.
		 */
		Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, int *index, int nI);
		
		/**
		 * Regular constructor. Just takes a list of positions. The cutoff and size
		 * are required for it to make nearest neighbor lists. Actually calls initializor 
		 * with same arguments.
		 */
		Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset);
		
		/**
		 * Constructor with mapped indices. Just like the regular constructor, but
		 * it also takes some mapped indices as arguments. Actually calls initializor with 
		 * same arguments.
		 */
		Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset, int *index, int nI);
		
		/**
		 * Regular initializor. Just takes a list of positions. The cutoff and size
		 * are required for it to make nearest neighbor lists.
		 */
		void initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size);
		
		/**
		 * Initializor with mapped indices. Just like the regular constructor, but
		 * it also takes some mapped indices as arguments.
		 */
		void initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, int *index, int nI);
		
		/**
		 * Regular initializor. Just takes a list of positions. The cutoff and size
		 * are required for it to make nearest neighbor lists.
		 */
		void initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset);
		
		/**
		 * Initializor with mapped indices. Just like the regular constructor, but
		 * it also takes some mapped indices as arguments.
		 */
		void initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset, int *index, int nI);
		
		/**
		 * Destructor. It should work, but if you have an issue with it deleting,
		 * then you probably messed something up.
		 */
		~Cell();
		
		/**
		 * Deletes all memory allocated by the object. Typically, a user won't need this.
		 */
		void deleteAllocation();
		
		/**
		 * Builds linked lists. Make sure to do this prior to calling neighbors.
		 */
		void build();
		
		/**
		 * One way mass query. 
		 * This one is thread safe.
		 */
		template <typename COMPOBJ>
		void oneWayCompare(COMPOBJ &compare);
		
		/**
		 * One way mass query. 
		 * This one is thread safe.
		 */
		template <typename COMPOBJ>
		void oneWayCompareIndex(COMPOBJ &compare);
		
		/**
		 * One way mass query with periodic boundaries. 
		 * This one is thread safe.
		 */
		template <typename COMPOBJ>
		void oneWayComparePeriodic(COMPOBJ compare);
		
		/**
		 * Two way mass query. 
		 * No, it is not thread safe.
		 */
		template <typename COMPOBJ>
		COMPOBJ twoWayCompare(COMPOBJ compare);
		
		/**
		 * Two way mass query with periodic boundaries. 
		 * No, it is not thread safe.
		 */
		template <typename COMPOBJ>
		COMPOBJ twoWayComparePeriodic(COMPOBJ compare);
		
		/**
		 * Two way mass query. 
		 * No, it is not thread safe. Passed index version.
		 */
		template <typename COMPOBJ>
		void twoWayCompareIndex(COMPOBJ &compare);
		
		/**
		 * Two way mass query with periodic boundaries. 
		 * No, it is not thread safe. Passed index version.
		 */
		template <typename COMPOBJ>
		void twoWayComparePeriodicIndex(COMPOBJ &compare);
		
		/**
		 * Query function with particle index and minimum image vector as an argument.
		 */
		int query(int particle, threeVector<T> &minImg);
		
		/**
		 * Query function with particle index as an argument.
		 */
		int query(int particle);
		
		/**
		 * Query only particles nearby with a higher index. Utilizes minimum image vector.
		 */
		int queryHalf(int particle, threeVector<T> &minImg);
		
		/**
		 * Query only particles nearby with a higher index.
		 */
		int queryHalf(int particle);
		
		/**
		 * Query using mapped indices. Utilizes minimum image vector.
		 */
		int queryIndex(int particle, threeVector<T> &minImg);
		
		/**
		 * Query with a position as an argument. Utilizes minimum image vector.
		 */
		int query(position<T> pos, threeVector<T> &minImg);
		
		/**
		 * Query with a position as an argument.
		 */
		int query(position<T> pos);
		
		/**
		 * Grabs a single particle from a particular cell. If you use the same cell, it will access the next one.
		 */
		int queryCell(int cell);
		
		/**
		 * Get a nearest neighbor until they are all exhausted.
		 */
		int getNext();
		
		/**
		 * Resets currentParticle and other states for query functions, per thread.
		 */
		void resetState();
		
		/**
		 * Set a new system size. Updates the cell sizes and such.
		 */
		void resize(threeVector<T> size);
		
		/**
		 * Change the number of positions available to the lists.
		 */
		void changeN(int nI, int *index);
	private:
		threeVector<int> nextBox[13];//Nearby boxes that can be accessed in an always forward fashion. No duplicate boxes, except the current box.
		
		position<T> *p;//Particle positions.
		int nP;//Number of particle positions.
		T cutoff;//Maximum radius of interactions.
		threeVector<T> s;//Size of the system at hand.
		threeVector<T> lOffset;//left hand offset, right hand offset is s+lOffset
		
		//two ways to make this work,
		//create indices then fill Cells while updating linked list normally
		//or create indices then fill Cells while sorting linked list as a key-value list, like a hash
		int *cells;//Our cell head list.
		int *linkedList;//The linked list.
		int *indices;//Indices for particles. E.g. which cell it belongs to. I know, terrible name, due to the index variable
		
		int *cellsIndex;//Cells with mapped index labels for elements.
		int *linkedListIndex;//Linked list with mapped index labels for elements
		
		int *flag;
		
		#ifdef LOW_DENSITY
			int *fullCells;
			int nFullCells;
		#endif
		
		threeVector<int> neighbor;//Neighbor variable.
		int currentCell;//Keeps track of current cell being called.
		int nState;//neighbor state
		int currentParticle;
		position<T> cPos;
		fourVector<int> inquiryCell;
		int inquiryParticle;
		int nextBoxIndex;
		
		fourVector<int> nCells;
		int nCellsAllocated;
		threeVector<T> cellSize;
		
		//flagged with null and 0
		int *pIndex;
		int nI;
};

template <typename T>
Cell<T>::Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size)
{
	this->initialize();
	this->initialize(particles,nParticles,cutoff,size);
}

template <typename T>
Cell<T>::Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, int *index, int nI)
{
	this->initialize();
	this->initialize(particles,nParticles,cutoff,size,index,nI);
}

template <typename T>
Cell<T>::Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset)
{
	this->initialize();
	this->initialize(particles,nParticles,cutoff,size,leftOffset);
}

template <typename T>
Cell<T>::Cell(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset, int *index, int nI)
{
	this->initialize();
	this->initialize(particles,nParticles,cutoff,size,leftOffset,index,nI);
}

template <typename T>
void Cell<T>::initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size)
{
	//this->deleteAllocation();
	this->p=particles;
	if(this->p==NULL)
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
	
	this->nCellsAllocated=nCells.t*2;
	this->cells=new int[nCellsAllocated];
	this->flag=new int[nCellsAllocated];
	this->linkedList=new int[nP];
	this->indices=new int[nP];
	#ifdef LOW_DENSITY
		this->nFullCells=0;
		this->fullCells=new int[nCellsAllocated];
	#endif
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	this->lOffset=0;
}


//The initialization routines should really be condensed, if an option is added, maybe it should be a regular public option
template <typename T>
void Cell<T>::initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, int *index, int nI)
{
	this->initialize(particles,nParticles,cutoff,size);
	this->pIndex=index;
	this->nI=nI;
	this->linkedListIndex=new int[this->nP];
	this->cellsIndex=new int[this->nCellsAllocated];
	this->lOffset=0;
}

template <typename T>
void Cell<T>::initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset)
{
	this->initialize(particles,nParticles,cutoff,size);
	this->lOffset=leftOffset;
}

template <typename T>
void Cell<T>::initialize(position<T> *particles, int nParticles, T cutoff, threeVector<T> size, threeVector<T> leftOffset, int *index, int nI)
{
	this->initialize(particles,nParticles,cutoff,size);
	this->pIndex=index;
	this->nI=nI;
	this->lOffset=leftOffset;
	this->linkedListIndex=new int[this->nP];
	this->cellsIndex=new int[this->nCellsAllocated];
}

template <typename T>
Cell<T>::~Cell()
{
	this->deleteAllocation();
}

template <typename T>
void Cell<T>::deleteAllocation()
{
	if(this->cells!=NULL)
		delete this->cells;
	if(this->linkedList!=NULL)
		delete this->linkedList;
	if(this->indices!=NULL)
		delete this->indices;
	if(this->cellsIndex!=NULL)
		delete this->cellsIndex;
	if(this->linkedListIndex!=NULL)
		delete this->linkedListIndex;
	if(this->flag!=NULL)
		delete this->flag;
	this->cells=NULL;
	this->linkedList=NULL;
	this->indices=NULL;
	this->cellsIndex=NULL;
	this->linkedListIndex=NULL;
	this->flag=NULL;
	
	#ifdef LOW_DENSITY
		if(this->fullCells!=NULL)
			delete this->fullCells;
		this->fullCells=NULL;
	#endif
}

//single threaded only!
template <typename T>
void Cell<T>::build()
{
	if(nI==0 && pIndex==NULL)
	{
		//empties the boxes
		for(int i=0;i<nCells.t;i++)
		{
			cells[i]=-1;
			flag[i]=0;
		}
		
		//calculate each particle's cell hash
		for(int i=0;i<nP;i++)
		{
			int x=int((p[i].x-lOffset.x)/cellSize.x);
			int y=int((p[i].y-lOffset.y)/cellSize.y);
			int z=int((p[i].z-lOffset.z)/cellSize.z);
			indices[i]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[i]=i;//index of the particle
		}
		
		//fill cells
		for(int i=0;i<nP;i++)
		{
			//push operation
			int current=indices[i];
			int buf=cells[current];
			cells[current]=linkedList[i];
			linkedList[i]=buf;
			if(linkedList[i]==1 && i==1)
			{
				std::cout << "wrong!\n";
				std::cout << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\n';
				std::cin.get();
			}
			#ifdef SOLVENT_FLAG
				if(p[i].type!=SOLVENT_FLAG)
					flag[current]=1;
			#else
				flag[current]=1;
			#endif
		}
	}
	else
	{
		//empties the boxes
		for(int i=0;i<nCells.t;i++)
		{
			cells[i]=-1;
			cellsIndex[i]=-1;
			flag[i]=0;
		}
		
		//calculate each particle's cell hash
		for(int i=0;i<nI;i++)
		{
			int j=pIndex[i];
			
			int x=int((p[j].x-lOffset.x)/cellSize.x);
			int y=int((p[j].y-lOffset.y)/cellSize.y);
			int z=int((p[j].z-lOffset.z)/cellSize.z);
			indices[j]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[j]=j;//index of the particle
			linkedListIndex[i]=i;
		}
		
		//fill cells
		for(int i=0;i<nI;i++)
		{
			int j=pIndex[i];
			//push operation
			int current=indices[j];
			int buf=cells[current];
			cells[current]=linkedList[j];
			linkedList[j]=buf;
			
			buf=cellsIndex[current];
			cellsIndex[current]=linkedListIndex[i];
			linkedListIndex[i]=buf;
			
			#ifdef SOLVENT_FLAG
				if(p[j].type!=SOLVENT_FLAG)
					flag[current]=1;
			#else
				flag[current]=1;
			#endif
		}
	}
	
	#ifdef LOW_DENSITY
		nFullCells=0;
		//only full cells
		for(int i=0;i<nCells.t;i++)
		{
			if(flag[i]==1 && cells[i]!=-1)
			{
				//This took forever to hunt down. If you use solvent flag without solvent, it is
				// a 50% performance hit, yikes!
				#ifdef SOLVENT_FLAG
					int x=i%nCells.x;
					int y=int(i/nCells.x)%nCells.y;
					int z=int(i/(nCells.x*nCells.y));
					//Here is the reason, we are flagging all the cells nearby, which is about
					// 13-20 more cells than needed to carry out this calculation.
					for(int j=0;j<27;j++)
					{
						int aa=j%3-1;
						int bb=int(j/3)%3-1;
						int cc=int(j/9)-1;
						int k=((x+aa+nCells.x)%nCells.x)+\
						(((y+bb+nCells.y)%nCells.y)*nCells.x)+\
						(((z+cc+nCells.z)%nCells.z)*nCells.x*nCells.y);
						
						if(flag[k]==0)
						{
							flag[k]=2;
							fullCells[nFullCells]=k;
							nFullCells++;
						}
					}
					fullCells[nFullCells++]=i;
				#endif
				#ifndef SOLVENT_FLAG
					fullCells[nFullCells++]=i;
				#endif
			}
		}
	#endif
	
	this->resetState();
}

//this can be executed in parallel or not
template <typename T>
void Cell<T>::resetState()
{
	inquiryParticle=-1;
	cPos.x=-1;
	cPos.y=-1;
	cPos.z=-1;
	currentCell=0;
	nState=0;
}

template <typename T>
template <typename COMPOBJ>
inline void Cell<T>::oneWayCompareIndex(COMPOBJ &compare)
{
	#ifdef LOW_DENSITY
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the currentbox full?
		if(cells[i]!=-1 && flag[i]!=0)
		{
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<27;j++)
			{
				int a=j%3-1;
				int b=int(j/3)%3-1;
				int c=int(j/9)-1;
				
				//next box
				int nextX=x+a;
				int nextY=y+b;
				int nextZ=z+c;
				
				if(nextX<nCells.x && nextX>=0 && nextY<nCells.y && nextY>=0 && nextZ<nCells.z && nextZ>=0)
				{
					int k=nextX+(nextY*nCells.x)+(nextZ*nCells.x*nCells.y);
					//is the cell empty?
					if(cells[k]!=-1 && flag[k]!=0)
						//compare all of next cell to all of current cell
						for(int m=cells[i]; m!=-1; m=linkedList[m])
							for(int n=cells[k]; n!=-1; n=linkedList[n])
								compare(m,n);
				}
			}
		}
	}
}

//For god sakes, please remember to optimize this, just like the periodic version.
template <typename T>
template <typename COMPOBJ>
inline COMPOBJ Cell<T>::twoWayCompare(COMPOBJ compare)
{
	int nextBoxStatic[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	
	#ifdef LOW_DENSITY
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the currentbox full?
		if(cells[i]!=-1 && flag[i]!=0)
		{
			//compare within a cell
			for(int m=cells[i]; m!=-1; m=linkedList[m])
				for(int n=linkedList[m]; n!=-1; n=linkedList[n])
					compare(p[m],p[n]);
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<13;j++)
			{
				//next box
				int nextX=x+nextBoxStatic[j][0];
				int nextY=y+nextBoxStatic[j][1];
				int nextZ=z+nextBoxStatic[j][2];
				
				if(nextX<=nCells.x-1 && nextX>=0 && nextY<=nCells.y-1 && nextY>=0 && nextZ<=nCells.z-1 && nextZ>=0)
				{
					int k=nextX+(nextY*nCells.x)+(nextZ*nCells.x*nCells.y);
					
					//is the cell empty?
					if(cells[k]!=-1 && flag[k]!=0)
						//compare all of next cell to all of current cell
						for(int m=cells[i]; m!=-1; m=linkedList[m])
							for(int n=cells[k]; n!=-1; n=linkedList[n])
								compare(p[m],p[n]);
				}
			}
		}
	}
	return compare;
}

template <typename T>
template <typename COMPOBJ>
inline void Cell<T>::twoWayCompareIndex(COMPOBJ &compare)
{
	int nextBoxStatic[26][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	
	#ifdef LOW_DENSITY
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the currentbox full?
		if(cells[i]!=-1 && flag[i]!=0)
		{
			//compare within a cell
			for(int m=cells[i]; m!=-1; m=linkedList[m])
				for(int n=linkedList[m]; n!=-1; n=linkedList[n])
					compare(m,n);
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<13;j++)
			{
				//next box
				int nextX=x+nextBoxStatic[j][0];
				int nextY=y+nextBoxStatic[j][1];
				int nextZ=z+nextBoxStatic[j][2];
				
				if(nextX<=nCells.x-1 && nextX>=0 && nextY<=nCells.y-1 && nextY>=0 && nextZ<=nCells.z-1 && nextZ>=0)
				{
					int k=nextX+(nextY*nCells.x)+(nextZ*nCells.x*nCells.y);
					
					//is the cell empty?
					if(cells[k]!=-1 && flag[k]!=0)
						//compare all of next cell to all of current cell
						for(int m=cells[i]; m!=-1; m=linkedList[m])
							for(int n=cells[k]; n!=-1; n=linkedList[n])
								compare(m,n);
				}
			}
		}
	}
}

template <typename T>
template <typename COMPOBJ>
inline COMPOBJ Cell<T>::twoWayComparePeriodic(COMPOBJ compare)
{
	int nextBoxStatic[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	
	#ifdef LOW_DENSITY
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the currentbox full?
		if(cells[i]!=-1 && flag[i]!=0)
		{
			//An optimization
			position<T> tempP[MAX_CELL_SIZE];
			
			//compare within a cell
			
			int ni=0;
			
			for(int m=cells[i]; m!=-1; ni++, m=linkedList[m])
			{
				#ifdef CELL_SIZE_FAILURE
					if(ni>MAX_CELL_SIZE)
					{
						std::cout << "Error(CellOpt): Too many particles in cell!\n";
						throw 0;//I don't know if this even works multithreaded
					}
				#endif
				
				//load an element, this linearizes the operations from here on out
				tempP[ni]=p[m];
				
				//compare current element of current cell to previous elements
				for(int j=0;j<ni;j++)
					compare(tempP[ni],tempP[j]);
			}
			/*
			for(int m=cells[i]; m!=-1; m=linkedList[m])
				for(int n=linkedList[m]; n!=-1; n=linkedList[n])
					compare(p[m],p[n]);
			*/
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<13;j++)
			{
				//these are private to this thread
				threeVector<T> minImg;
				
				//next box
				int nextX=x+nextBoxStatic[j][0];
				int nextY=y+nextBoxStatic[j][1];
				int nextZ=z+nextBoxStatic[j][2];
				
				//minimum image
				minImg.x=0;
				if(nextX<0) minImg.x-=s.x;
				if(nextX>=nCells.x) minImg.x+=s.x;
				
				minImg.y=0;
				if(nextY<0) minImg.y-=s.y;
				if(nextY>=nCells.y) minImg.y+=s.y;
				
				minImg.z=0;
				if(nextZ<0) minImg.z-=s.z;
				if(nextZ>=nCells.z) minImg.z+=s.z;
				
				int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1 && flag[i]!=0)
				{
					//compare all of next cell to all of current cell
					for(int n=cells[k]; n!=-1; n=linkedList[n])
					{
						position<T> pn=p[n];
						pn.x+=minImg.x;
						pn.y+=minImg.y;
						pn.z+=minImg.z;
						
						//for(int m=cells[i]; m!=-1; m=linkedList[m])
						for(int m=0;m<ni;m++)
						{
							compare(tempP[m],pn);
						}
					}
				}
			}
		}
	}
	//std::cout << "Ready to Destroy compare!\n";
	//std::cin.get();
	return compare;
}

template <typename T>
template <typename COMPOBJ>
inline void Cell<T>::twoWayComparePeriodicIndex(COMPOBJ &compare)
{
	int nextBoxStatic[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	
	#ifdef LOW_DENSITY
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the currentbox full?
		if(cells[i]!=-1 && flag[i]!=0)
		{
			
			for(int m=cells[i]; m!=-1; m=linkedList[m])
				for(int n=linkedList[m]; n!=-1; n=linkedList[n])
					compare(m,n);
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<13;j++)
			{
				//these are private to this thread
				threeVector<T> minImg;
				
				//next box
				int nextX=x+nextBoxStatic[j][0];
				int nextY=y+nextBoxStatic[j][1];
				int nextZ=z+nextBoxStatic[j][2];
				
				//minimum image
				minImg.x=0;
				if(nextX<0) minImg.x-=s.x;
				if(nextX>=nCells.x) minImg.x+=s.x;
				
				minImg.y=0;
				if(nextY<0) minImg.y-=s.y;
				if(nextY>=nCells.y) minImg.y+=s.y;
				
				minImg.z=0;
				if(nextZ<0) minImg.z-=s.z;
				if(nextZ>=nCells.z) minImg.z+=s.z;
				
				int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1 && flag[i]!=0)
				{
					//compare all of next cell to all of current cell
					for(int n=cells[k]; n!=-1; n=linkedList[n])
					{
						//for(int m=cells[i]; m!=-1; m=linkedList[m])
						for(int m=cells[i]; m!=-1; m=linkedList[m])
						{
							compare(m,n,minImg);
						}
					}
				}
			}
		}
	}
}

//this can be executed in parallel or not
template <typename T>
int Cell<T>::query(int particle, threeVector<T> &minImg)
{
	//issueance needs to be fast, try to keep branching statements fairly unested
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int((p[inquiryParticle].x-lOffset.x)/cellSize.x);
			inquiryCell.y=int((p[inquiryParticle].y-lOffset.y)/cellSize.y);
			inquiryCell.z=int((p[inquiryParticle].z-lOffset.z)/cellSize.z);
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

//#pragma ECHO "Hello!"

//this can be executed in parallel or not
template <typename T>
int Cell<T>::queryIndex(int particle, threeVector<T> &minImg)
{
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			int i=pIndex[inquiryParticle];
			inquiryCell.x=int((p[i].x-lOffset.x)/cellSize.x);
			inquiryCell.y=int((p[i].y-lOffset.y)/cellSize.y);
			inquiryCell.z=int((p[i].z-lOffset.z)/cellSize.z);
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
			currentParticle=linkedListIndex[currentParticle];
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
			currentParticle=cellsIndex[currentCell];
		}
			
		return currentParticle;
	//#endif
}

//this can be executed in parallel or not, this doesn't use minimum image,
// but it does return particles in cells periodically, just remember to ignore them
template <typename T>
int Cell<T>::query(int particle)
{
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int((p[inquiryParticle].x-lOffset.x)/cellSize.x);
			inquiryCell.y=int((p[inquiryParticle].y-lOffset.y)/cellSize.y);
			inquiryCell.z=int((p[inquiryParticle].z-lOffset.z)/cellSize.z);
			inquiryCell.t=inquiryCell.x+inquiryCell.y*nCells.x+inquiryCell.z*nCells.x*nCells.y;
			
			//set neighbor
			//neighbor.x=-2;
			//neighbor.y=-1;
			//neighbor.z=-1;
			
			nState=0;
			
			currentParticle=-1;
		}
		else
		{
			if(currentParticle==1)
			{
				std::cout << currentParticle << '\t' << linkedList[currentParticle] << '\n';
				std::cin.get();
			}
			//get next particle in linked list
			currentParticle=linkedList[currentParticle];
		}
		
		//check if it is last particle in list
		while(currentParticle==-1)// && nState<27)
		{
			//go to next cell, this increments through each axis
			//neighbor.x++;
			//if(neighbor.x==2)//cache misses are going to be bad with this
			//{
			//	neighbor.x=-1;
			//	neighbor.y++;
			//	if(neighbor.y==2)
			//	{
			//		neighbor.y=-1;
			//		neighbor.z++;
			//		if(neighbor.z==2)
			//		{
			//			//last neighbor
			//			currentParticle=-1;
			//			return currentParticle;
			//		}
			//	}
			//}
			
			nState++;
			
			neighbor.x=nState%3-1;
			neighbor.y=int(nState/3)%3-1;
			neighbor.z=int(nState/9)-1;
			
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

//this can be executed in parallel or not
template <typename T>
int Cell<T>::queryHalf(int particle, threeVector<T> &minImg)
{
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int((p[inquiryParticle].x-lOffset.x)/cellSize.x);
			inquiryCell.y=int((p[inquiryParticle].y-lOffset.y)/cellSize.y);
			inquiryCell.z=int((p[inquiryParticle].z-lOffset.z)/cellSize.z);
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
}

//this can be executed in parallel or not, doesn't use minimum image
template <typename T>
int Cell<T>::queryHalf(int particle)
{
		//reset current cell parameters
		if(particle!=inquiryParticle)
		{
			//set inquiryParticle
			inquiryParticle=particle;
			
			//set inquiryCell
			inquiryCell.x=int((p[inquiryParticle].x-lOffset.x)/cellSize.x);
			inquiryCell.y=int((p[inquiryParticle].y-lOffset.y)/cellSize.y);
			inquiryCell.z=int((p[inquiryParticle].z-lOffset.z)/cellSize.z);
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
}

//this can be executed in parallel or not
template <typename T>
int Cell<T>::query(position<T> pos, threeVector<T> &minImg)
{
		//reset current cell parameters
		if(pos.x != cPos.x && pos.y != cPos.y && pos.z != cPos.z)
		{
			//set cPos
			cPos=pos;
			
			//set inquiryCell
			inquiryCell.x=int(cPos.x/cellSize.x);
			inquiryCell.y=int(cPos.y/cellSize.y);
			inquiryCell.z=int(cPos.z/cellSize.z);
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
}

//this can be executed in parallel or not, this doesn't use minimum image,
// but it does return particles in cells periodically, just remember to ignore them
template <typename T>
int Cell<T>::query(position<T> pos)
{
		//reset current cell parameters
		if(pos.x != cPos.x && pos.y != cPos.y && pos.z != cPos.z)
		{
			//set cPos
			cPos=pos;
			
			//set inquiryCell
			inquiryCell.x=int(cPos.x/cellSize.x);
			inquiryCell.y=int(cPos.y/cellSize.y);
			inquiryCell.z=int(cPos.z/cellSize.z);
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
}

/*
//get the next cell using a hash index
template <typename T>
int Cell<T>::queryCell(int cell)
{
	if(cell!=currentCell)
	{
		currentCell=cell;
		currentParticle=cells[currentCell];
		nextBoxIndex=0;
		inquiryCell.x=currentCell%nCells.x;
		inquiryCell.y=int(currentCell/nCells.x)%nCells.y;
		inquiryCell.z=int(currentCell/(nCells.x*nCells.y))%nCells.z;
	}
	currentParticle=linkedList[currentParticle];
	return currentParticle;
}

//grab an anonymous cell
template <typename T>
int Cell<T>::queryCell()
{
	//get current particle
	currentParticle=linkedList[currentParticle];
	
	//is it the last particle in the list?
	while(currentParticle==-1 && currentCell<nCells.t)
	{
		//grab current particle in cell
		currentParticle=cells[currentCell++];
		
		//reset the box index list
		nextBoxIndex=0;
		
		//set the inquiry cell
		inquiryCell.x=currentCell%nCells.x;
		inquiryCell.y=int(currentCell/nCells.x)%nCells.y;
		inquiryCell.z=int(currentCell/(nCells.x*nCells.y))%nCells.z;
		
		//try this again
		currentParticle=linkedList[currentParticle];
	}
	
	//return next particle if possible
	return currentParticle;
}

//grab an adjacent particle
template <typename T>
int Cell<T>::getNext()
{
	//get next particle
	nextParticle=cells[nextCell];
	
	//is it the last particle in the list?
	while(nextParticle==-1 && nextBoxIndex<13)
	{
		//grab next set of coordinates
		threeVector<int> neighbor=nextBox[nextBoxIndex++];
		
		//get next cell
		neighbor.x+=inquiryCell.x;
		neighbor.y+=inquiryCell.y;
		neighbor.z+=inquiryCell.z;
		
		//set next cell
		nextCell=((neighbor.x+nCells.x)%nCells.x)+
			(((neighbor.y+nCells.y)%nCells.y)*nCells.x)+
			(((neighbor.z+nCells.z)%nCells.z)*nCells.x*nCells.y);
			
		//try this again
		nextParticle=cells[nextCell];
	}
	
	//return next particle if possible
	return nextParticle;
}
*/
//This changes the number of particles by changing the redirection indices
template <typename T>
void Cell<T>::changeN(int nI, int *index)
{
	this->nI=nI;
	this->pIndex=index;
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
	if(nCellsAllocated<nCells.t)
	{
		delete cells;
		cells=new int[nCells.t*2];
		nCellsAllocated=nCells.t*2;
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
