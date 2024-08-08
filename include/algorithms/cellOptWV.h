//Optimal algorithm for calculating forces and potentials between neighbors
#include "functions.h"

#ifndef MD_CELLOPT
#define MD_CELLOPT


//number of maximum cell particles, expand as neccesary
#ifndef MAX_CELL_SIZE
#define MAX_CELL_SIZE 512
#endif

/*
//this one is really generic, P is the positions type, C is the compute object, and LL is the linked list object
template <typename P, typename C, typename LL>
class GenericCell {
	public:
		GenericCell(P *points, LL *linkedlist, C *comp);
		~GenericCell();
		void compute();
		void build();
	private:
		P *p;
		LL *ll;
		C *c;
};
*/

//This algorithm is most efficient for about 43 particles per cell, the efficiency levels off pas that point.
//For less than 43 particles per cell, the threads should be decomposed to a particle per thread rather than
// a cell per thread.
template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
class CellOpt {
	public:
		CellOpt()//Used for initializing null lists
		{
			cells=NULL;
			#ifdef LOW_DENSITY
				fullCells=NULL;
			#endif
			flag=NULL;
			linkedList=NULL;
			indices=NULL;
			#ifdef _OPENMP
				lockedCells=NULL;//used to lock a cell for writing
			#endif
		};
		///DPD
		CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed);
		CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax);
		CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, int nIndices, int *index);
		CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax, int nIndices, int *index);
		~CellOpt();
		
		void initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed);
		void initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax);
		void initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, int nIndices, int *index);
		void initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax, int nIndices, int *index);
		
		void computeForce();//compute forces
		void resize(threeVector<T> size);
		void changeNP(int nParticles);
		//T computeCustom();
		//the custom function, don't use this with a compute intensive section, like accelerations
		//virtual T function();
		
		//builds linked lists, do this before any function
		void build();
	private:
		position<T> *p;//particle positions
		threeVector<T> *a;//accelerations
		threeVector<T> *v;//velocities, usually null for MD, points to velocities for DPD
		T *cF;//particle force constants
		int nP;
		int nT;
		T cutoff;//a value if there is one, 0 otherwise
		threeVector<T> s;
		threeVector<bool> w;
		MTRand **randNum;
		
		//two ways to make this work,
		//create indices then fill Cells while updating linked list normally
		//or create indices then fill Cells while sorting linked list as a key-value list, like a hash
		int *cells;
		#ifdef LOW_DENSITY
			int *fullCells;
			int nFullCells;
		#endif
		bool *flag;//solvent flag
		int *linkedList;
		int *indices;
		#ifdef _OPENMP
			omp_lock_t *lockedCells;//used to lock a cell for writing
		#endif
		
		fourVector<int> nCells;
		int nCellsOld;
		threeVector<T> cellSize;
		
		//flagged with null and 0
		int *index;
		int nI;
};


template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
CellOpt<T,FORCE>::CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed)
{
	this->initialize(particles,velocities,accelerations,Fconstants,nParticles,nTypes,size,wrap,cutoff,seed);
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
CellOpt<T,FORCE>::CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax)
{
	this->initialize(particles,velocities,accelerations,Fconstants,nParticles,nTypes,size,wrap,cutoff,seed,dMax);
}


template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
CellOpt<T,FORCE>::CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, int nIndices, int *index)
{
	this->initialize(particles,velocities,accelerations,Fconstants,nParticles,nTypes,size,wrap,cutoff,seed,nIndices,index);
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
CellOpt<T,FORCE>::CellOpt(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax, int nIndices, int *index)
{
	this->initialize(particles,velocities,accelerations,Fconstants,nParticles,nTypes,size,wrap,cutoff,seed,dMax,nIndices,index);
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed)
{
	this->p=particles;
	this->v=velocities;
	this->a=accelerations;
	this->cF=Fconstants;
	if(p==NULL)
	{
		std::cout << "p is null in CellOpt!\n";
		throw 0;
	}
	if(v==NULL)
	{
		std::cout << "v is null in CellOpt!\n";
		throw 0;
	}
	if(a==NULL)
	{
		std::cout << "a is null in CellOpt!\n";
		throw 0;
	}
	if(cF==NULL)
	{
		std::cout << "cF is null! in CellOpt\n";
		throw 0;
	}
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;

	if(this->cutoff<=0)
	{
		std::cout << "No cutoff in CellOpt!\n";
		throw 0;
	}
	this->nCells.x=int(this->s.x/this->cutoff);
	this->nCells.y=int(this->s.y/this->cutoff);
	this->nCells.z=int(this->s.z/this->cutoff);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t;
	
	this->cells=new int[nCells.t*2];
	this->flag=new bool[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	this->linkedList=new int[nP];
	this->indices=new int[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	#ifdef _OPENMP
		this->lockedCells=NULL;
		//this->lockedCells=(omp_lock_t *)malloc(nCells.t*2*sizeof(omp_lock_t));
		this->lockedCells=new omp_lock_t[nCells.t*2];
		if(lockedCells!=NULL)
			#pragma omp parallel for
			for(int i=0;i<nCells.t;i++)
				omp_init_lock(&lockedCells[i]);
	#endif
	
	#ifdef _OPENMP
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
	
	this->index=NULL;
	this->nI=0;
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax)
{
	this->p=particles;
	this->v=velocities;
	this->a=accelerations;
	this->cF=Fconstants;
	if(p==NULL)
	{
		std::cout << "p is null in CellOpt!\n";
		throw 0;
	}
	if(v==NULL)
	{
		std::cout << "v is null in CellOpt!\n";
		throw 0;
	}
	if(a==NULL)
	{
		std::cout << "a is null in CellOpt!\n";
		throw 0;
	}
	if(cF==NULL)
	{
		std::cout << "cF is null! in CellOpt\n";
		throw 0;
	}
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;

	if(this->cutoff<=0)
		throw "No cutoff!\n";
	this->nCells.x=int(this->s.x/(this->cutoff+dMax));
	this->nCells.y=int(this->s.y/(this->cutoff+dMax));
	this->nCells.z=int(this->s.z/(this->cutoff+dMax));
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t;
	
	this->cells=new int[nCells.t*2];
	this->flag=new bool[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	this->linkedList=new int[nP];
	this->indices=new int[nP];

	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	#ifdef _OPENMP
		//this->lockedCells=(omp_lock_t *)malloc(nCells.t*2*sizeof(omp_lock_t));
		this->lockedCells=new omp_lock_t[nCells.t*2];
		#pragma omp parallel for
		for(int i=0;i<nCells.t;i++)
			omp_init_lock(&lockedCells[i]);
	#endif
	
	#ifdef _OPENMP
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
	
	this->index=NULL;
	this->nI=0;
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, int nIndices, int *index)
{
	this->p=particles;
	this->v=velocities;
	this->a=accelerations;
	this->cF=Fconstants;
	if(p==NULL)
	{
		std::cout << "p is null in CellOpt!\n";
		throw 0;
	}
	if(a==NULL)
	{
		std::cout << "a is null in CellOpt!\n";
		throw 0;
	}
	if(v==NULL)
	{
		std::cout << "v is null in CellOpt!\n";
		throw 0;
	}
	if(cF==NULL)
	{
		std::cout << "cF is null! in CellOpt\n";
		throw 0;
	}
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;

	if(this->cutoff<=0)
	{
		std::cout << "No cutoff in CellOpt!\n";
		throw 0;
	}
	this->nCells.x=int(this->s.x/this->cutoff);
	this->nCells.y=int(this->s.y/this->cutoff);
	this->nCells.z=int(this->s.z/this->cutoff);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t;
	
	this->cells=new int[nCells.t*2];
	this->flag=new bool[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	this->linkedList=new int[nP];
	this->indices=new int[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	#ifdef _OPENMP
		this->lockedCells=NULL;
		//this->lockedCells=(omp_lock_t *)malloc(nCells.t*2*sizeof(omp_lock_t));
		this->lockedCells=new omp_lock_t[nCells.t*2];
		if(lockedCells!=NULL)
			#pragma omp parallel for
			for(int i=0;i<nCells.t;i++)
				omp_init_lock(&lockedCells[i]);
	#endif
	
	#ifdef _OPENMP
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
	
	this->index=index;
	this->nI=nIndices;
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::initialize(position<T> *particles, threeVector<T> *velocities, threeVector<T> *accelerations, T *Fconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int seed, T dMax, int nIndices, int *index)
{
	this->p=particles;
	this->v=velocities;
	this->a=accelerations;
	this->cF=Fconstants;
	if(p==NULL)
	{
		std::cout << "p is null in CellOpt!\n";
		throw 0;
	}
	if(a==NULL)
	{
		std::cout << "a is null in CellOpt!\n";
		throw 0;
	}
	if(v==NULL)
	{
		std::cout << "v is null in CellOpt!\n";
		throw 0;
	}
	if(cF==NULL)
	{
		std::cout << "cF is null! in CellOpt\n";
		throw 0;
	}
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;

	if(this->cutoff<=0)
		throw "No cutoff!\n";
	this->nCells.x=int(this->s.x/(this->cutoff+dMax));
	this->nCells.y=int(this->s.y/(this->cutoff+dMax));
	this->nCells.z=int(this->s.z/(this->cutoff+dMax));
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t;
	
	this->cells=new int[nCells.t*2];
	this->flag=new bool[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	this->linkedList=new int[nP];
	this->indices=new int[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	#ifdef _OPENMP
		//this->lockedCells=(omp_lock_t *)malloc(nCells.t*2*sizeof(omp_lock_t));
		this->lockedCells=new omp_lock_t[nCells.t*2];
		#pragma omp parallel for
		for(int i=0;i<nCells.t;i++)
			omp_init_lock(&lockedCells[i]);
	#endif
	
	#ifdef _OPENMP
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
			
	this->index=index;
	this->nI=nIndices;
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
CellOpt<T,FORCE>::~CellOpt()
{
	if(cells!=NULL)
		delete this->cells;
	#ifdef LOW_DENSITY
		if(fullCells!=NULL)
			delete this->fullCells;
	#endif
	if(linkedList!=NULL)
		delete this->linkedList;
	if(indices!=NULL)
		delete this->indices;
	if(flag!=NULL)
		delete this->flag;
	#ifdef _OPENMP
		if(lockedCells!=NULL)
		{
			#pragma omp parallel for
			for(int i=0;i<nCells.t;i++)
				omp_destroy_lock(&lockedCells[i]);
			delete this->lockedCells;
		}
		//free(this->lockedCells);
	#endif
	#ifdef _OPENMP
		for(int i=0;i<omp_get_max_threads();i++)
			delete randNum[i];
	#else
		delete randNum[0];
	#endif
	delete randNum;
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::build()
{
	int nextBox[13][3]={{-1,0,0},{0,-1,0},{-1,-1,0},{0,0,-1},{-1,0,-1},{0,-1,-1},{-1,-1,-1},{1,0,-1},{0,1,-1},{1,-1,-1},{-1,1,-1},{1,1,-1},{-1,1,0}};
	//non-sorting version
	int buf,current,x,y,z;
	
	//empties the boxes
	//#pragma omp parallel for
	for(int i=0;i<nCells.t;i++)
	{
		cells[i]=-1;
		flag[i]=false;
	}
	
	if(nI==0 && index==NULL)
	{
		//calculate each particle's cell hash
		#pragma omp parallel for private(x,y,z)
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
					
								int adjacent=((nextX+nCells.x)%nCells.x)+\
									(((nextY+nCells.y)%nCells.y)*nCells.x)+\
									(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
								flag[adjacent]=true;
							}
						}
					}
				}
			#endif
			#ifndef SOLVENT_FLAG
				flag[current]=true;
			#endif
			buf=cells[current];
			cells[current]=linkedList[i];
			linkedList[i]=buf;
		}
	}
	else
	{
		//calculate each particle's cell hash
		#pragma omp parallel for private(x,y,z)
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			if(j>=nP || j<0)
			{
				std::cout << "index " << i << " isn't in range (" << j << "?nParticles)\n";
				throw 0;
			}
			x=int(p[j].x/cellSize.x);
			y=int(p[j].y/cellSize.y);
			z=int(p[j].z/cellSize.z);
			indices[j]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[j]=j;//index of the particle
		}
		
		//fill cells
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			//push operation
			current=indices[j];
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
						
								int adjacent=((nextX+nCells.x)%nCells.x)+\
								(((nextY+nCells.y)%nCells.y)*nCells.x)+\
								(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
								flag[adjacent]=true;
							}
						}
					}
				}
			#endif
			#ifndef SOLVENT_FLAG
				flag[current]=true;
			#endif
			buf=cells[current];
			cells[current]=linkedList[j];
			linkedList[j]=buf;
		}
	}
	
	#ifdef LOW_DENSITY
		nFullCells=0;
		//only full cells
		for(int i=0;i<nCells.t;i++)
		{
			if(flag[i] && cells[i]!=-1)
			{
				fullCells[nFullCells]=i;
				nFullCells++;
			}
		}
	#endif
}

template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::computeForce()
{
	int nextBox[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	

	//every cell is decomposed to a thread, not per row or anything else
	//this should have each atom decomposed to a thread for better results on low density systems, possibly
	#ifdef LOW_DENSITY
		#pragma omp parallel for
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		#pragma omp parallel for
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the current cell empty?
		if(flag[i] && cells[i]!=-1)
		{
			//all variables that are private to this thread
			threeVector<T> minImg;
			int x,y,z,j,k,m,n,ni,cindex,nextX,nextY,nextZ;
			//T dx,dy,dz,dr,magnitude;
			threeVector<T> tempA[MAX_CELL_SIZE];
			threeVector<T> tempV[MAX_CELL_SIZE];
			position<T> tempP[MAX_CELL_SIZE];
			threeVector<T> currentA,currentV;
			position<T> currentP;
			
			//load each element of current list while computing accelerations
			for(ni=0, m=cells[i]; m!=-1; ni++, m=linkedList[m])
			{
				#ifdef CELL_SIZE_FAILURE
				if(ni>MAX_CELL_SIZE)
				{
					std::cout << "Cell size: " << ni << std::endl;
					throw 0;//I don't know if this even works multithreaded
				}
				#endif
				
				//load an element of current box and initialize
				tempP[ni]=p[m];
				tempV[ni]=v[m];
				tempA[ni].x=0;
				tempA[ni].y=0;
				tempA[ni].z=0;
				
				//comapre it to previously loaded elements
				for(j=0;j<ni;j++)
				{
					#ifdef _OPENMP

						FORCE(cutoffSquared,nT, tempP[ni], tempP[j], tempV[ni], tempV[j], tempA[ni], tempA[j], cF, 
						      randNum[omp_get_thread_num()]->rand53());
					#else
						FORCE(cutoffSquared,nT, tempP[ni], tempP[j], tempV[ni], tempV[j], tempA[ni], tempA[j], cF,
						      randNum[0]->rand53());
					#endif
				}
			}
			
			//boxes around current
			for(j=0;j<13;j++)
			{
				//current box
				x=i%nCells.x;
				y=int(i/nCells.x)%nCells.y;
				z=int(i/(nCells.x*nCells.y));
				
				//next box
				nextX=x+nextBox[j][0];
				nextY=y+nextBox[j][1];
				nextZ=z+nextBox[j][2];
				
				//minimum image, only works if cell crosses boundary and if wrap, aka w.component, is true
				minImg.x=0;
				if(nextX<0 && w.x) minImg.x-=s.x;
				if(nextX>=nCells.x && w.x) minImg.x+=s.x;
				
				minImg.y=0;
				if(nextY<0 && w.y) minImg.y-=s.y;
				if(nextY>=nCells.y && w.y) minImg.y+=s.y;
				
				minImg.z=0;
				if(nextZ<0 && w.z) minImg.z-=s.z;
				if(nextZ>=nCells.z && w.z) minImg.z+=s.z;
				
				k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1)
				{
					#ifdef _OPENMP
						//lock cell if openmp is used
						omp_set_lock(&lockedCells[k]);
					#endif
					//compare all elements in current cell to every element of next cell
					for(n=cells[k]; n!=-1; n=linkedList[n])
					{
						//this is confusing, currentP is actually an element of the next box
						currentP=p[n];
						currentV=v[n];
						//minimum image, if applicable
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
						currentA.x=0;
						currentA.y=0;
						currentA.z=0;
						//compare element to everything in current box
						for(m=0;m<ni;m++)
						{
							#ifdef _OPENMP
								FORCE(cutoffSquared,nT, tempP[m], currentP, tempV[m], currentV, tempA[m], currentA, cF, 
								      randNum[omp_get_thread_num()]->rand53());
							#else
								FORCE(cutoffSquared,nT, tempP[m], currentP, tempV[m], currentV, tempA[m], currentA, cF, 
								      randNum[0]->rand53());
							#endif
						}
						
						//reduce next acceleration element back to array
						a[n].x+=currentA.x;
						a[n].y+=currentA.y;
						a[n].z+=currentA.z;
					}
					#ifdef _OPENMP
						//remove lock from cell if openmp is used
						omp_unset_lock(&lockedCells[k]);
					#endif
				}
			}
			#ifdef _OPENMP
				//lock current cell if openmp is used
				omp_set_lock(&lockedCells[i]);
			#endif
			//reduce the current accelerations of cell back to array
			for(ni=0, m=cells[i]; m!=-1; ni++, m=linkedList[m])
			{
				a[m].x+=tempA[ni].x;
				a[m].y+=tempA[ni].y;
				a[m].z+=tempA[ni].z;
			}
			#ifdef _OPENMP
				//unlock current cell if openmp is used
				omp_unset_lock(&lockedCells[i]);
			#endif
		}
	}
}

//This changes the number of particles
template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::changeNP(int nParticles)
{
	nP=nParticles;//seriously, this is all this does... Can't reallocate the number of particles, just change where it ends.
}

//this doesn't adjust the sizes of fullCells and lockedCells, so it's broken until further notice
//further notice: fixed
template <typename T, void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &v1, threeVector<T> &v2, threeVector<T> &a1, threeVector<T> &a2,T *fC, T R)>
void CellOpt<T,FORCE>::resize(threeVector<T> size)
{
	int nCellsOld=nCells.t;
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
		#ifdef LOW_DENSITY
			delete fullCells;
			fullCells=new int[nCells.t*2];
		#endif
		#ifdef _OPENMP
			delete lockedCells;
			lockedCells=NULL;
			//this->lockedCells=(omp_lock_t *)malloc(nCells.t*2*sizeof(omp_lock_t));
			lockedCells=new omp_lock_t[nCells.t*2];
			if(lockedCells!=NULL)
				#pragma omp parallel for
				for(int i=0;i<nCells.t;i++)
					omp_init_lock(&lockedCells[i]);
		#endif
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
