/**
 * \brief Optimal algorithm for molecular dynamics calculations of non-bonded forces and potentials between neighbors
 * Utilizes cell lists in the form of linked lists. Uses hash to determine which cell each particle is in.
 * It is optimized for low density, high density, OpenMP, and CPUs with large caches. Actually, I can't verify
 * that last statement; the testing to make sure it utilizes the cache correctly hasn't been done, but the 
 * compiler should at least be directed to optimize it that way. The theory is that a local variable, with a 
 * size MAX_CELL_SIZE, sits on the cache, while the elements of the list sit in the RAM. If I copy needed elements
 * directly to the local memory, then I can speed up memory access with sequential access during calculations and 
 * also maintain thread safety by locking the memory when I copy back. Regardless of the sequential access speedup,
 * the thread local storage definately works.
 * 
 * There are flags for optimizations: 
 *  1. '#define LOW_DENSITY' if you want to use low density optimizations, otherwise it assumes all cells are full.
 *  2. '#define MAX_CELL_SIZE N' if the average number of particles per cell is much lower than N.
 *  3. '#define WARNINGS_ENABLED' if you want to see warnings that affect the system but not the execution. Unpredictable particle movement, etc...
 *  4. '#define ERRORS_ENABLED' if you want to see errors that affect the execution of the program. Seg. Faults, etc...
 */

#include "functions.h"

#ifndef MD_CELLOPT
#define MD_CELLOPT


//number of maximum cell particles, expand as neccesary
#ifndef MAX_CELL_SIZE
#define MAX_CELL_SIZE 512
#endif

/**
 * \brief Optimized cell object for molecular dynamics simulations utilizing a simple pairwise potential.
 * This algorithm is most efficient for about 43 particles per cell, the efficiency gain levels off past that point.
 * For less than 43 particles per cell, the threads should be decomposed to a particle per thread rather than
 *  a cell per thread.
 */
template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), 
	void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
class CellOpt {
	public:
		///Constructor that generates an object that is empty.
		CellOpt()//Used for initializing null lists
		{
			cells=NULL;
			#ifdef LOW_DENSITY
				fullCells=NULL;
			#endif
			flag=NULL;
			hashIndices=NULL;
		};
		///Constructor utilizing a cutoff.
		///Calls initialize with same arguments.
		CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff);
		///Constructor utilizing a cutoff, larger by an addition of dMax. Note that it doesn't check dMax for positive numbers.
		///Calls initialize with same arguments.
		CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax);
		///Constructor utilizing a cutoff and list of indices.
		///Calls initialize with same arguments.
		CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int nIndices, int *index);
		///Constructor utilizing a cutoff, larger by an additon of dMax. Also uses a list of indices. Note that it doesn't check dMax for positive numbers.
		///Calls initialize with same arguments.
		CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax, int nIndices, int *index);
		///Destructor.
		~CellOpt();
		
		///Initializer utilizing a cutoff.
		void initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff);
		///Initializer utilizing a cutoff, larger by an addition of dMax. Note that it doesn't check dMax for positive numbers.
		void initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax);
		///Initializer utilizing a cutoff and list of indices.
		void initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int nIndices, int *index);
		///Initializer utilizing a cutoff, larger by an addition of dMax. Also uses a list of indices. Note that it doesn't check dMax for positive numbers.
		void initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants,
		     int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax, int nIndices, int *index);
		
		///Compute pair-wise force interactions.
		void computeForce();//compute forces
		///Compute pair-wise potential interactions.
		T computePotential();//do everything
		///Compute a difference in potential based on a scaling of size by a factor a.
		T computeDPotential(threeVector<T>);
		///Compute a pair-wise potential for a single particle versus all it's neighbors.
		T computePotential(int f);//do one
		///Compute pair-wise potential interactions and output to a file named 'potential_name.dat'
		T outputPotential(T time, char *name);
		///Resize system for this object. Recalculates memory allocation and constants.
		void resize(threeVector<T> size);
		///Change the number of particles in the object. Does no other allocation. Used for building different lists with the same object.
		void changeNP(int nParticles);
		//T computeCustom();
		//the custom function, don't use this with a compute intensive section, like accelerations
		//virtual T function();
		
		///Builds cell lists. Do this before any compute function.
		void build();
	private:
		position<T> *p;//particle positions
		threeVector<T> *a;//accelerations
		T *cU;//particle potential constants
		T *cF;//particle force constants
		int nP;
		int nT;
		T cutoff;//a value if there is one, 0 otherwise
		threeVector<T> s;
		threeVector<bool> w;
		
		//two ways to make this work,
		//create indices then fill Cells while updating linked list normally
		//or create indices then fill Cells while sorting linked list as a key-value list, like a hash
		int *cells;
		#ifdef LOW_DENSITY
			int *fullCells;
			int nFullCells;
		#endif
		char *flag;//solvent flag
		keyVal<int,int> *hashIndices;//sorted by hash and index, which are key and value respectively
		
		fourVector<int> nCells;
		int nCellsOld;
		threeVector<T> cellSize;
		
		//flagged with null and 0
		int *index;
		int nI;
};

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff)
{
	this->p=particles;
	this->a=accelerations;
	this->cU=Uconstants;
	this->cF=Fconstants;
	#ifdef ERRORS_ENABLED
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
		if(cU==NULL)
		{
			std::cout << "cU is null in CellOpt!\n";
			throw 0;
		}
		if(cF==NULL)
		{
			std::cout << "cF is null! in CellOpt\n";
			throw 0;
		}
	#endif
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;
	#ifdef ERRORS_ENABLED
		if(this->cutoff<=0)
		{
			std::cout << "No cutoff in CellOpt!\n";
			throw 0;
		}
	#endif
	this->nCells.x=int(this->s.x/this->cutoff);
	this->nCells.y=int(this->s.y/this->cutoff);
	this->nCells.z=int(this->s.z/this->cutoff);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t*2;
	
	this->cells=new int[nCells.t*2];
	this->flag=new char[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	hashIndices=new keyVal<int,int>[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->index=NULL;
	this->nI=0;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax)
{
	this->p=particles;
	this->a=accelerations;
	this->cU=Uconstants;
	this->cF=Fconstants;
	#ifdef ERRORS_ENABLED
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
		if(cU==NULL)
		{
			std::cout << "cU is null in CellOpt!\n";
			throw 0;
		}
		if(cF==NULL)
		{
			std::cout << "cF is null! in CellOpt\n";
			throw 0;
		}
	#endif
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;
	#ifdef ERRORS_ENABLED
		if(this->cutoff<=0)
		{
			std::cout << "No cutoff in CellOpt!\n";
			throw 0;
		}
	#endif
	this->nCells.x=int(this->s.x/(this->cutoff+dMax));
	this->nCells.y=int(this->s.y/(this->cutoff+dMax));
	this->nCells.z=int(this->s.z/(this->cutoff+dMax));
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t*2;
	
	this->cells=new int[nCells.t*2];
	this->flag=new char[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	hashIndices=new keyVal<int,int>[nP];

	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->index=NULL;
	this->nI=0;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int nIndices, int *index)
{
	this->p=particles;
	this->a=accelerations;
	this->cU=Uconstants;
	this->cF=Fconstants;
	#ifdef ERRORS_ENABLED
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
		if(cU==NULL)
		{
			std::cout << "cU is null in CellOpt!\n";
			throw 0;
		}
		if(cF==NULL)
		{
			std::cout << "cF is null! in CellOpt\n";
			throw 0;
		}
	#endif
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;
	#ifdef ERRORS_ENABLED
		if(this->cutoff<=0)
		{
			std::cout << "No cutoff in CellOpt!\n";
			throw 0;
		}
	#endif
	this->nCells.x=int(this->s.x/this->cutoff);
	this->nCells.y=int(this->s.y/this->cutoff);
	this->nCells.z=int(this->s.z/this->cutoff);
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t*2;
	
	this->cells=new int[nCells.t*2];
	this->flag=new char[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	hashIndices=new keyVal<int,int>[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->index=index;
	this->nI=nIndices;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::initialize(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax, int nIndices, int *index)
{
	this->p=particles;
	this->a=accelerations;
	this->cU=Uconstants;
	this->cF=Fconstants;
	#ifdef ERRORS_ENABLED
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
		if(cU==NULL)
		{
			std::cout << "cU is null in CellOpt!\n";
			throw 0;
		}
		if(cF==NULL)
		{
			std::cout << "cF is null! in CellOpt\n";
			throw 0;
		}
	#endif
	this->nP=nParticles;
	this->nT=nTypes;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;
	#ifdef ERRORS_ENABLED
		if(this->cutoff<=0)
		{
			std::cout << "No cutoff in CellOpt!\n";
			throw 0;
		}
	#endif
	this->nCells.x=int(this->s.x/(this->cutoff+dMax));
	this->nCells.y=int(this->s.y/(this->cutoff+dMax));
	this->nCells.z=int(this->s.z/(this->cutoff+dMax));
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	this->nCellsOld=nCells.t*2;
	
	this->cells=new int[nCells.t*2];
	this->flag=new char[nCells.t*2];
	#ifdef LOW_DENSITY
		this->fullCells=new int[nCells.t*2];
	#endif
	hashIndices=new keyVal<int,int>[nP];
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->index=index;
	this->nI=nIndices;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
CellOpt<T,POTENTIAL,FORCE>::CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff)
{
	this->initialize(particles,accelerations,Fconstants,Uconstants,nParticles,nTypes,size,wrap,cutoff);
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
CellOpt<T,POTENTIAL,FORCE>::CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax)
{
	this->initialize(particles,accelerations,Fconstants,Uconstants,nParticles,nTypes,size,wrap,cutoff,dMax);
}


template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
CellOpt<T,POTENTIAL,FORCE>::CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, int nIndices, int *index)
{
	this->initialize(particles,accelerations,Fconstants,Uconstants,nParticles,nTypes,size,wrap,cutoff,nIndices,index);
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
CellOpt<T,POTENTIAL,FORCE>::CellOpt(position<T> *particles, threeVector<T> *accelerations, T *Fconstants, T *Uconstants, int nParticles, int nTypes, threeVector<T> size, threeVector<bool> wrap, T cutoff, T dMax, int nIndices, int *index)
{
	this->initialize(particles,accelerations,Fconstants,Uconstants,nParticles,nTypes,size,wrap,cutoff,dMax,nIndices,index);
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
CellOpt<T,POTENTIAL,FORCE>::~CellOpt()
{
	if(cells!=NULL)
		delete this->cells;
	#ifdef LOW_DENSITY
		if(fullCells!=NULL)
			delete this->fullCells;
	#endif
	if(flag!=NULL)
		delete this->flag;
	
	if(hashIndices!=NULL)
		delete hashIndices;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::build()
{
	int nextBox[13][3]={{-1,0,0},{0,-1,0},{-1,-1,0},{0,0,-1},{-1,0,-1},{0,-1,-1},{-1,-1,-1},{1,0,-1},{0,1,-1},{1,-1,-1},{-1,1,-1},{1,1,-1},{-1,1,0}};
	//non-sorting version
	//empties the boxes
	for(int i=0;i<nCells.t;i++)
	{
		cells[i]=-1;
		flag[i]=0;
	}
	
	if(nI==0 && index==NULL)
	{
		
		//calculate each particle's cell hash
		for(int i=0;i<nP;i++)
		{
			int x=int(p[i].x/cellSize.x);
			int y=int(p[i].y/cellSize.y);
			int z=int(p[i].z/cellSize.z);
			
			//correction for older systems, usually when loading
			x=(x>=nCells.x)?x-1:x;
			y=(y>=nCells.y)?y-1:y;
			z=(z>=nCells.z)?z-1:z;
			
			#ifdef ERRORS_ENABLED
				if(x>=nCells.x || y>=nCells.y || z>=nCells.z || x<0 || y<0 || z<0)
				{
					std::cout << "Error(cellOpt): cell placement is on boundary!\n";
					std::cout << nCells.x << '\t' << nCells.y << '\t' << nCells.z << ":\t";
					std::cout << x << '\t' << y << '\t' << z << '\n';
					std::cout << s.x << '\t' << s.y << '\t' << s.z << ":\t";
					std::cout << p[i].type << '\t' << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\n';
					throw 0;
				}
			#endif
			
			hashIndices[i].value=i;
			hashIndices[i].key=x+y*nCells.x+z*nCells.x*nCells.y;
		}
		
		//fill cells via Linked List
		for(int i=0;i<nP;i++)
		{
			//push operation
			int current=hashIndices[i].key;
			#ifdef SOLVENT_FLAG
				if(p[i].type!=SOLVENT_FLAG)
					flag[current]=1;
			#else
				flag[current]=1;
			#endif
			int buf=cells[current];
			cells[current]=hashIndices[i].value;
			hashIndices[i].value=buf;
		}
	}
	else
	{
		//calculate each particle's cell hash
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			
			#ifdef ERRORS_ENABLED
				if(j>=nP || j<0)
				{
					std::cout << "index " << i << " isn't in range (" << j << "?nParticles)\n";
					throw 0;
				}
			#endif
			
			int x=int(p[j].x/cellSize.x);
			int y=int(p[j].y/cellSize.y);
			int z=int(p[j].z/cellSize.z);
			
			//correction for older systems, usually when loading
			x=(x>=nCells.x)?x-1:x;
			y=(y>=nCells.y)?y-1:y;
			z=(z>=nCells.z)?z-1:z;
			
			#ifdef ERRORS_ENABLED
				if(x>=nCells.x || y>=nCells.y || z>=nCells.z || x<0 || y<0 || z<0)
				{
					std::cout << "Error(cellOpt): cell placement is on boundary!\n";
					std::cout << nCells.x << '\t' << nCells.y << '\t' << nCells.z << ":\t";
					std::cout << x << '\t' << y << '\t' << z << '\n';
					std::cout << s.x << '\t' << s.y << '\t' << s.z << ":\t";
					std::cout << p[i].type << '\t' << p[i].x << '\t' << p[i].y << '\t' << p[i].z << '\n';
					throw 0;
				}
			#endif
			
			hashIndices[j].value=j;
			hashIndices[j].key=x+y*nCells.x+z*nCells.x*nCells.y;
		}
		
		//fill cells
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			//push operation
			int current=hashIndices[j].key;
			#ifdef SOLVENT_FLAG
				if(p[j].type!=SOLVENT_FLAG)
					flag[current]=1;
			#else
				flag[current]=1;
			#endif
			int buf=cells[current];
			cells[current]=hashIndices[j].value;
			hashIndices[j].value=buf;
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
					fullCells[nFullCells]=i;
					nFullCells++;
					#ifdef ERRORS_ENABLED
						if(nFullCells>nCells.t)
						{
							std::cout << "Too many fullCells!\n";
							throw 0;
						}
					#endif
				#endif
				#ifndef SOLVENT_FLAG
					fullCells[nFullCells]=i;
					nFullCells++;
				#endif
			}
		}
	#endif
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::computeForce()
{
	int nextBox[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	
	//every cell is decomposed to a thread, not per row or anything else
	//this should have each atom decomposed to a thread for better results on low density systems, possibly
	#ifdef LOW_DENSITY
		for(int ijk=0;ijk<nFullCells;ijk++)
	#else
		for(int i=0;i<nCells.t;i++)
	#endif
	{
		#ifdef LOW_DENSITY
			int i=fullCells[ijk];
		#endif
		//is the current cell empty?
		if(flag[i]!=0 && cells[i]!=-1)
		{
			//all variables that are private to this thread
			threeVector<T> minImg;
			int m,ni;
			//T dx,dy,dz,dr,magnitude;
			threeVector<T> tempA[MAX_CELL_SIZE];
			position<T> tempP[MAX_CELL_SIZE];
			threeVector<T> currentA;
			position<T> currentP;
			
			//load each element of current list while computing accelerations
			for(ni=0, m=cells[i]; m!=-1; ni++, m=hashIndices[m].value)
			{
				#ifdef CELL_SIZE_FAILURE
					if(ni>MAX_CELL_SIZE)
					{
						std::cout << "Error(CellOpt): Too many particles in cell!\n";
						throw 0;//I don't know if this even works multithreaded
					}
				#endif
				
				//load an element of current box and initialize
				tempP[ni]=p[m];
				tempA[ni].x=0;
				tempA[ni].y=0;
				tempA[ni].z=0;
				
				//comapre it to previously loaded elements
				for(int j=0;j<ni;j++)
					FORCE(cutoffSquared,nT, tempP[ni], tempP[j], tempA[ni], tempA[j], cF);
			}
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			//boxes around current
			for(int j=0;j<13;j++)
			{
				//next box
				int nextX=x+nextBox[j][0];
				int nextY=y+nextBox[j][1];
				int nextZ=z+nextBox[j][2];
				
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
				
				int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1)
				{
					for(int n=cells[k]; n!=-1; n=hashIndices[n].value)
					{
						//this is confusing, currentP is actually an element of the next box
						currentP=p[n];
						//minimum image, if applicable
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
						currentA.x=0;
						currentA.y=0;
						currentA.z=0;
						//compare element to everything in current box
						for(m=0;m<ni;m++)
							FORCE(cutoffSquared,nT, tempP[m], currentP, tempA[m], currentA, cF);
						
						//reduce next acceleration element back to array
						a[n].x+=currentA.x;
						a[n].y+=currentA.y;
						a[n].z+=currentA.z;
					}
				}
			}
			//reduce the current accelerations of cell back to array
			for(ni=0, m=cells[i]; m!=-1; ni++, m=hashIndices[m].value)
			{
				a[m].x+=tempA[ni].x;
				a[m].y+=tempA[ni].y;
				a[m].z+=tempA[ni].z;
			}
		}
	}
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
T CellOpt<T,POTENTIAL,FORCE>::computePotential()
{
	int nextBox[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	
	T potential=0;
	//every cell is decomposed to a thread, not per row or anything else
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
			//these are private to this thread
			threeVector<T> minImg;
			int m,ni;
			//T dx,dy,dz,dr,magnitude;
			position<T> tempP[MAX_CELL_SIZE];
			position<T> currentP;
			
			//load current cell
			for(ni=0, m=cells[i]; m!=-1; ni++, m=hashIndices[m].value)
			{
				#ifdef CELL_SIZE_FAILURE
					if(ni>MAX_CELL_SIZE)
					{
						std::cout << "Error(CellOpt): Too many particles in cell!\n";
						throw 0;//I don't know if this even works multithreaded
					}
				#endif
				
				//load an element
				tempP[ni]=p[m];
				
				//compare current element of current cell to previous elements
				for(int j=0;j<ni;j++)
					potential+=POTENTIAL(cutoffSquared, nT, tempP[ni], tempP[j], cU);
			}
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<13;j++)
			{
				//next box
				int nextX=x+nextBox[j][0];
				int nextY=y+nextBox[j][1];
				int nextZ=z+nextBox[j][2];
				
				//minimum image
				minImg.x=0;
				if(nextX<0 && w.x) minImg.x-=s.x;
				if(nextX>=nCells.x && w.x) minImg.x+=s.x;
				
				minImg.y=0;
				if(nextY<0 && w.y) minImg.y-=s.y;
				if(nextY>=nCells.y && w.y) minImg.y+=s.y;
				
				minImg.z=0;
				if(nextZ<0 && w.z) minImg.z-=s.z;
				if(nextZ>=nCells.z && w.z) minImg.z+=s.z;
				
				int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1)
				{
					//compare all of next cell to all of current cell
					for(int n=cells[k]; n!=-1; n=hashIndices[n].value)
					{
						//this is confusing, currentP is actually an element of the next box
						currentP=p[n];
						//minimum image, if applicable
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
						//all of current cell
						for(int m=0;m<ni;m++)
							potential+=POTENTIAL(cutoffSquared, nT, tempP[m], currentP, cU);
					}
				}
			}
		}
	}
	
	return potential;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
T CellOpt<T,POTENTIAL,FORCE>::computeDPotential(threeVector<T> aSize)
{
	int nextBox[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	
	T dPotential=0;
	//every cell is decomposed to a thread, not per row or anything else
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
			//these are private to this thread
			threeVector<T> minImg;
			int m,ni;
			//T dx,dy,dz,dr,magnitude;
			position<T> tempP[MAX_CELL_SIZE];
			position<T> currentP;
			
			//load current cell
			for(ni=0, m=cells[i]; m!=-1; ni++, m=hashIndices[m].value)
			{
				#ifdef CELL_SIZE_FAILURE
					if(ni>MAX_CELL_SIZE)
					{
						std::cout << "Error(CellOpt): Too many particles in cell!\n";
						throw 0;//I don't know if this even works multithreaded
					}
				#endif
				
				//load an element
				tempP[ni]=p[m];
				
				//compare current element of current cell to previous elements
				for(int j=0;j<ni;j++)
				{
					double oldPotential=POTENTIAL(cutoffSquared, nT, tempP[ni], tempP[j], cU);
					position<T> pA=tempP[ni];
					position<T> pB=tempP[j];
					pA.x*=aSize.x;
					pA.y*=aSize.y;
					pA.z*=aSize.z;
					pB.x*=aSize.x;
					pB.y*=aSize.y;
					pB.z*=aSize.z;
					double newPotential=POTENTIAL(cutoffSquared, nT, pA, pB, cU);
					dPotential+=newPotential-oldPotential;
				}
			}
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
			
			for(int j=0;j<13;j++)
			{
				//next box
				int nextX=x+nextBox[j][0];
				int nextY=y+nextBox[j][1];
				int nextZ=z+nextBox[j][2];
				
				//minimum image
				minImg.x=0;
				if(nextX<0 && w.x) minImg.x-=s.x;
				if(nextX>=nCells.x && w.x) minImg.x+=s.x;
				
				minImg.y=0;
				if(nextY<0 && w.y) minImg.y-=s.y;
				if(nextY>=nCells.y && w.y) minImg.y+=s.y;
				
				minImg.z=0;
				if(nextZ<0 && w.z) minImg.z-=s.z;
				if(nextZ>=nCells.z && w.z) minImg.z+=s.z;
				
				int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1)
				{
					//compare all of next cell to all of current cell
					for(int n=cells[k]; n!=-1; n=hashIndices[n].value)
					{
						//this is confusing, currentP is actually an element of the next box
						currentP=p[n];
						//minimum image, if applicable
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
						//all of current cell
						for(int m=0;m<ni;m++)
						{
							double oldPotential=POTENTIAL(cutoffSquared, nT, tempP[m], currentP, cU);
							position<T> pA=tempP[m];
							position<T> pB=currentP;
							pA.x*=aSize.x;
							pA.y*=aSize.y;
							pA.z*=aSize.z;
							pB.x*=aSize.x;
							pB.y*=aSize.y;
							pB.z*=aSize.z;
							double newPotential=POTENTIAL(cutoffSquared, nT, pA, pB, cU);
							dPotential+=newPotential-oldPotential;
						}
					}
				}
			}
		}
	}
	
	return dPotential;
}

template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
T CellOpt<T,POTENTIAL,FORCE>::computePotential(int f)
{
	T cutoffSquared=cutoff*cutoff;
	
	T potential=0;
	fourVector<int> currentCell;
	currentCell.x=p[f].x/cellSize.x;
	currentCell.y=p[f].y/cellSize.y;
	currentCell.z=p[f].z/cellSize.z;
	currentCell.t=(currentCell.x)+(nCells.x*currentCell.y)+(nCells.x*nCells.y*currentCell.z);
	
	//is the currentbox full?
	if(cells[currentCell.t]!=-1 && flag[currentCell.t]!=0)
	{
		//these are private to this thread
		position<T> currentP;
		threeVector<T> minImg;
		
		//check current cell
		for(int m=cells[currentCell.t]; m!=-1; m=hashIndices[m].value)
			if(m!=f)
				potential+=POTENTIAL(cutoffSquared, nT, p[m], p[f], cU);

		for(int a=-1;a<2;a++)
		for(int b=-1;b<2;b++)
		for(int c=-1;c<2;c++)
		{
			//next box
			int nextX=currentCell.x+a;//nextBox[j][0];
			int nextY=currentCell.y+b;//nextBox[j][1];
			int nextZ=currentCell.z+c;//nextBox[j][2];
			
			//minimum image
			minImg.x=0;
			if(nextX<0 && w.x) minImg.x-=s.x;
			if(nextX>=nCells.x && w.x) minImg.x+=s.x;
			
			minImg.y=0;
			if(nextY<0 && w.y) minImg.y-=s.y;
			if(nextY>=nCells.y && w.y) minImg.y+=s.y;
			
			minImg.z=0;
			if(nextZ<0 && w.z) minImg.z-=s.z;
			if(nextZ>=nCells.z && w.z) minImg.z+=s.z;
			
			int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
			//is the cell empty?
			if(cells[k]!=-1 && k!=currentCell.t)
			{
				//compare all of next cell to all of current cell
				for(int n=cells[k]; n!=-1; n=hashIndices[n].value)
				{
					//this is confusing, currentP is actually an element of the next box
					currentP=p[n];
					//minimum image, if applicable
					currentP.x+=minImg.x;
					currentP.y+=minImg.y;
					currentP.z+=minImg.z;
					potential+=POTENTIAL(cutoffSquared, nT, p[f], currentP, cU);
				}
			}
		}
	}
	
	return potential;
}



template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
T CellOpt<T,POTENTIAL,FORCE>::outputPotential(T time, char *name)
{
	int nextBox[13][3]={{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1},{1,1,1},{-1,0,1},{0,-1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,-1,0}};
	T cutoffSquared=cutoff*cutoff;
	
	T potential=0;
	
	//every cell is decomposed to a thread, not per row or anything else
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
			//these are private to this thread
			threeVector<T> minImg;
			int m,ni;
			//T dx,dy,dz,dr,magnitude;
			position<T> tempP[MAX_CELL_SIZE];
			position<T> currentP;
			
			//load current cell
			for(m=cells[i]; m!=-1; m=hashIndices[m].value)
			{
				#ifdef CELL_SIZE_FAILURE
					if(ni>MAX_CELL_SIZE)
					{
						std::cout << "Error(CellOpt): Too many particles in cell!\n";
						throw 0;//I don't know if this even works multithreaded
					}
				#endif
				
				//load an element
				tempP[ni]=p[m];
				
				//compare current element of current cell to previous elements
				for(int j=0;j<ni;j++)
					potential+=POTENTIAL(cutoffSquared, nT, tempP[ni], tempP[j], cU);
			}
			
			//current box
			int x=i%nCells.x;
			int y=int(i/nCells.x)%nCells.y;
			int z=int(i/(nCells.x*nCells.y));
				
			for(int j=0;j<13;j++)
			{
				//next box
				int nextX=x+nextBox[j][0];
				int nextY=y+nextBox[j][1];
				int nextZ=z+nextBox[j][2];
				
				//minimum image
				minImg.x=0;
				if(nextX<0 && w.x) minImg.x-=s.x;
				if(nextX>=nCells.x && w.x) minImg.x+=s.x;
				
				minImg.y=0;
				if(nextY<0 && w.y) minImg.y-=s.y;
				if(nextY>=nCells.y && w.y) minImg.y+=s.y;
				
				minImg.z=0;
				if(nextZ<0 && w.z) minImg.z-=s.z;
				if(nextZ>=nCells.z && w.z) minImg.z+=s.z;
				
				int k=((nextX+nCells.x)%nCells.x)+(((nextY+nCells.y)%nCells.y)*nCells.x)+(((nextZ+nCells.z)%nCells.z)*nCells.x*nCells.y);
				
				//is the cell empty?
				if(cells[k]!=-1)
				{
					//compare all of next cell to all of current cell
					for(int n=cells[k]; n!=-1; n=hashIndices[n].value)
					{
						//this is confusing, currentP is actually a particle from the next cell
						currentP=p[n];
						
						//minimum image, if applicable
						currentP.x+=minImg.x;
						currentP.y+=minImg.y;
						currentP.z+=minImg.z;
						//all of current cell
						for(int m=0;m<ni;m++)
							potential+=POTENTIAL(cutoffSquared, nT, tempP[m], currentP, cU);
					}
				}
			}
		}
	}
	
	std::fstream dataFile;
	std::string buf("potential_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << time << '\t' << potential << '\n';
	dataFile.close();
	
	return potential;
}

//This changes the number of particles
template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::changeNP(int nParticles)
{
	nP=nParticles;//seriously, this is all this does... Can't reallocate the number of particles, just change where it ends.
}

//this doesn't adjust the sizes of fullCells and lockedCells, so it's broken until further notice
//further notice: fixed
template <typename T, T POTENTIAL(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC), void FORCE(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)>
void CellOpt<T,POTENTIAL,FORCE>::resize(threeVector<T> size)
{
	//rescale number of cells
	nCells.x=int(size.x/cutoff);
	nCells.y=int(size.y/cutoff);
	nCells.z=int(size.z/cutoff);
	nCells.t=nCells.x*nCells.y*nCells.z;
	
	//check if new rescaling is higher than allocated number of boxes, this usually isn't the case
	if(nCellsOld<nCells.t)
	{
		delete cells;
		delete flag;
		cells=new int[nCells.t*2];
		flag=new char[nCells.t*2];
		
		#ifdef LOW_DENSITY
			delete fullCells;
			fullCells=new int[nCells.t*2];
		#endif
		nCellsOld=nCells.t*2;
	}
	
	//rescale cell size
	cellSize.x=size.x/nCells.x;
	cellSize.y=size.y/nCells.y;
	cellSize.z=size.z/nCells.z;
	
	//std::cout << "NCells: " << nCells.x << '\t' << nCells.y << '\t' << nCells.z << '\t' << nCells.t << '\n';
	//std::cout << "cellSize: " << cellSize.x << '\t' << cellSize.y << '\t' << cellSize.z << '\n';
	//std::cout << "size: " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	//rescale overall size
	s=size;
}


//end of MD_CELL
#endif
