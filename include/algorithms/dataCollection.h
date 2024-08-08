#include "functions.h"
//this whole file needs to be broken up into smaller parts
//various data collection classes
//These should handle their own output with some parameters,
//e.g. compute returns something, output places it in a file
//"Com" is a center of mass class
//"Radius" is a class for calculatint radius from a center of mass
//"Height" is a class for calculating the height from a plane
//"Diffusion" is a class for calculating the diffusion of particles
//"Volume" is a class for calculating volume properties
//"Surface" is a class for surface properties

#ifndef MD_DATACOLLECTION
#define MD_DATACOLLECTION

//Center of mass class and definitions
template <typename T>
class Com {
	public:
		Com(position<T> *particles, int nParticles, threeVector<T> size);
		Com(position<T> *particles, int nParticles, threeVector<T> size, int *index, int nIndices);
		~Com(){};
		
		threeVector<T> compute();//return a value
		threeVector<T> output(T time, char *name);//return a value and output to a file com_name.dat
		
	private:
		position<T> *p;
		int nP;
		threeVector<T> s;
		//these are flagged by null
		int *index;
		int nI;
};

template <typename T>
Com<T>::Com(position<T> *particles, int nParticles, threeVector<T> size)
{
	this->p=particles;
	this->nP=nParticles;
	this->s=size;
	this->index=NULL;
	this->nI=0;
}

template <typename T>
Com<T>::Com(position<T> *particles, int nParticles, threeVector<T> size, int *index, int nIndices)
{
	this->p=particles;
	this->nP=nParticles;
	this->s=size;
	this->index=index;
	this->nI=nIndices;
}

template <typename T>
threeVector<T> Com<T>::compute()
{
	threeVector<T> comVal;
	comVal.x=0;
	comVal.y=0;
	comVal.z=0;
	T x=0,y=0,z=0;
	if(nI==0 && index==NULL)
	{
		#pragma omp parallel for reduction(+:x) reduction(+:y) reduction(+:z)
		for(int i=0;i<nP;i++)
		{
			x+=p[i].x;
			y+=p[i].y;
			z+=p[i].z;
		}
		comVal.x=x;
		comVal.y=y;
		comVal.z=z;
		comVal.x/=T(nP);
		comVal.y/=T(nP);
		comVal.z/=T(nP);
	}
	else
	{
		#pragma omp parallel for reduction(+:x) reduction(+:y) reduction(+:z)
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			x+=p[j].x;
			y+=p[j].y;
			z+=p[j].z;
		}
		comVal.x=x;
		comVal.y=y;
		comVal.z=z;
		comVal.x/=T(nI);
		comVal.y/=T(nI);
		comVal.z/=T(nI);
	}
	
	return comVal;
}

template <typename T>
threeVector<T> Com<T>::output(T time, char *name)
{
	threeVector<T> comVal;
	comVal.x=0;
	comVal.y=0;
	comVal.z=0;
	T x=0,y=0,z=0;
	if(nI==0 && index==NULL)
	{
		#pragma omp parallel for reduction(+:x) reduction(+:y) reduction(+:z)
		for(int i=0;i<nP;i++)
		{
			x+=p[i].x;
			y+=p[i].y;
			z+=p[i].z;
		}
		comVal.x=x;
		comVal.y=y;
		comVal.z=z;
		comVal.x/=T(nP);
		comVal.y/=T(nP);
		comVal.z/=T(nP);
	}
	else
	{
		#pragma omp parallel for reduction(+:x) reduction(+:y) reduction(+:z)
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			x+=p[j].x;
			y+=p[j].y;
			z+=p[j].z;
		}
		comVal.x=x;
		comVal.y=y;
		comVal.z=z;
		comVal.x/=T(nI);
		comVal.y/=T(nI);
		comVal.z/=T(nI);
	}
	
	std::fstream dataFile;
	std::string buf("com_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << time << '\t' << comVal.x << '\t' << comVal.y << '\t' << comVal.z << std::endl;
	dataFile.close();
	
	return comVal;
}




//Radius class and definitions
template <typename T>
class Radius {
	public:
		Radius(position<T> *particles, int nParticles, threeVector<T> size);
		Radius(position<T> *particles, int nParticles, threeVector<T> size, int *index, int nIndices);
		~Radius(){};
		
		T compute(threeVector<T> comV);//all particles
		T output(threeVector<T> comV, T time, char *name);
		
	private:
		position<T> *p;
		int nP;
		threeVector<T> s;//size of system
		threeVector<T> comVal;
		//these are flagged by null
		int *index;
		int nI;
};

template <typename T>
Radius<T>::Radius(position<T> *particles, int nParticles, threeVector<T> size)
{
	this->p=particles;
	this->nP=nParticles;
	this->s=size;
	this->index=NULL;
	this->nI=0;
}

template <typename T>
Radius<T>::Radius(position<T> *particles, int nParticles, threeVector<T> size, int *index, int nIndices)
{
	this->p=particles;
	this->nP=nParticles;
	this->s=size;
	this->index=index;
	this->nI=nIndices;
}

template <typename T>
T Radius<T>::compute(threeVector<T> comVal)
{
	T radius=0;
	T dx,dy,dz,dr;
	
	if(nI==0 && index==NULL)
	{
		#pragma omp parallel for private(dx,dy,dz,dr) reduction(+:radius)
		for(int i=0;i<nP;i++)
		{
			//this is completely screwed up if anything crosses the boundary...
			dx=p[i].x-comVal.x;
			dy=p[i].y-comVal.y;
			dz=p[i].z-comVal.z;
			dr=sqrt(dx*dx+dy*dy+dz*dz);
			radius+=dr;
		}
		radius/=T(nP);
	}
	else
	{
		#pragma omp parallel for private(dx,dy,dz,dr) reduction(+:radius)
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			//this is completely screwed up if anything crosses the boundary...
			dx=p[j].x-comVal.x;
			dy=p[j].y-comVal.y;
			dz=p[j].z-comVal.z;
			dr=sqrt(dx*dx+dy*dy+dz*dz);
			radius+=dr;
		}
		radius/=T(nI);
	}
	
	return radius;
}

template <typename T>
T Radius<T>::output(threeVector<T> comVal, T time, char *name)
{
	T radius=0;
	T dx,dy,dz,dr;
	
	if(nI==0 && index==NULL)
	{
		#pragma omp parallel for private(dx,dy,dz,dr) reduction(+:radius)
		for(int i=0;i<nP;i++)
		{
			//this is completely screwed up if anything crosses the boundary...
			dx=p[i].x-comVal.x;
			dy=p[i].y-comVal.y;
			dz=p[i].z-comVal.z;
			dr=sqrt(dx*dx+dy*dy+dz*dz);
			radius+=dr;
		}
		radius/=T(nP);
	}
	else
	{
		#pragma omp parallel for private(dx,dy,dz,dr) reduction(+:radius)
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			//this is completely screwed up if anything crosses the boundary...
			dx=p[j].x-comVal.x;
			dy=p[j].y-comVal.y;
			dz=p[j].z-comVal.z;
			dr=sqrt(dx*dx+dy*dy+dz*dz);
			radius+=dr;
		}
		radius/=T(nI);
	}
	
	std::fstream dataFile;
	std::string buf("radius_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << time << '\t' << radius << std::endl;
	dataFile.close();
	
	return radius;
}


template <typename T, typename S>
class Neighbors {
	public:
		//list is position list
		//conformation is list of what conformation is being compared
		//cutoff is the cutoff distance for the position list
		Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
			  threeVector<T> size);
		Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
			  threeVector<T> size, threeVector<bool> wrap);
		Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
			  threeVector<T> size, int nIndices, int *index);
		Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
			  threeVector<T> size, threeVector<bool> wrap, int nIndices, int *index);
		
		~Neighbors();
		
		//builds linked lists, do this before any function
		void build();
		
		//compare is a function that is used to compare conformations 
		int compute(bool compare(S, S, T), T compareVal);//returns number of neighbor domains
		int output(bool compare(S,S), T time, char *name);//returns number of domains, and outputs the number of domains
		
		int *domains;//start indices of domains
		int nDomains;//number of domains
		position<T> *sorted;//sorted list by domain
	private:
		int sL;
		position<T> *l;//abstracted
		S *conf;//really abstracted
		//uses a cutoff list algorithm
		int *cells;
		int *linkedList;
		int *indices;
		int *stack;
		
		threeVector<T> s;
		threeVector<bool> w;
		T cutoff;
		fourVector<T>nCells;
		threeVector<T>cellSize;
		
		//flagged with null and 0
		int *index;
		int nI;
};

template <typename T, typename S>
Neighbors<T,S>::Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
		     threeVector<T> size)
{
	this->l=list;
	this->sL=sizeList;
	this->conf=conformation;
	this->cutoff=cutoff;
	this->s=size;

	if(this->cutoff<=0)
		throw "No cutoff!\n";
	this->nCells.x=this->s.x/this->cutoff;
	this->nCells.y=this->s.y/this->cutoff;
	this->nCells.z=this->s.z/this->cutoff;
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->w.x=false;
	this->w.y=false;
	this->w.z=false;
	
	this->cells=new int[nCells.t];
	this->linkedList=new int[sL];
	this->indices=new int[sL];
	this->stack=new int[sL];
	this->domains=new int[sL];
	this->sorted=new position<T>[sL];
	
	this->index=NULL;
	this->nI=0;
}

template <typename T, typename S>
Neighbors<T,S>::Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
		     threeVector<T> size, threeVector<bool> wrap)
{
	this->l=list;
	this->sL=sizeList;
	this->conf=conformation;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;
	
	if(this->cutoff<=0)
		throw "No cutoff!\n";
	this->nCells.x=this->s.x/this->cutoff;
	this->nCells.y=this->s.y/this->cutoff;
	this->nCells.z=this->s.z/this->cutoff;
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->cells=new int[nCells.t];
	this->linkedList=new int[sL];
	this->indices=new int[sL];
	this->stack=new int[sL];
	this->domains=new int[sL];
	this->sorted=new position<T>[sL];
	
	this->index=NULL;
	this->nI=0;
}

template <typename T, typename S>
Neighbors<T,S>::Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
		     threeVector<T> size, int nIndices, int *index)
{
	this->l=list;
	this->sL=sizeList;
	this->conf=conformation;
	this->cutoff=cutoff;
	this->s=size;
	
	if(this->cutoff<=0)
		throw "No cutoff!\n";
	this->nCells.x=this->s.x/this->cutoff;
	this->nCells.y=this->s.y/this->cutoff;
	this->nCells.z=this->s.z/this->cutoff;
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->w.x=false;
	this->w.y=false;
	this->w.z=false;
	
	this->cells=new int[nCells.t];
	this->linkedList=new int[sL];
	this->indices=new int[sL];
	this->stack=new int[sL];
	this->domains=new int[sL];
	this->sorted=new position<T>[sL];
	
	this->index=index;
	this->nI=nIndices;
}

template <typename T, typename S>
Neighbors<T,S>::Neighbors(position<T> *list, int sizeList, S *conformation, T cutoff,
		     threeVector<T> size, threeVector<bool> wrap, int nIndices, int *index)
{
	this->l=list;
	this->sL=sizeList;
	this->conf=conformation;
	this->cutoff=cutoff;
	this->s=size;
	this->w=wrap;
	
	if(this->cutoff<=0)
		throw "No cutoff!\n";
	this->nCells.x=this->s.x/this->cutoff;
	this->nCells.y=this->s.y/this->cutoff;
	this->nCells.z=this->s.z/this->cutoff;
	this->nCells.t=nCells.x*nCells.y*nCells.z;
	
	this->cellSize.x=this->s.x/nCells.x;
	this->cellSize.y=this->s.y/nCells.y;
	this->cellSize.z=this->s.z/nCells.z;
	
	this->cells=new int[nCells.t];
	this->linkedList=new int[sL];
	this->indices=new int[sL];
	this->stack=new int[sL];
	this->domains=new int[sL];
	this->sorted=new position<T>[sL];
	
	this->index=index;
	this->nI=nIndices;
}

template <typename T, typename S>
Neighbors<T,S>::~Neighbors()
{
	delete cells;
	delete linkedList;
	delete indices;
	delete stack;
	delete domains;
	delete sorted;
}

template <typename T, typename S>
void Neighbors<T,S>::build()
{
	//non-sorting version
	int buf,current,x,y,z;
	
	//empties the boxes
	//#pragma omp parallel for
	for(int i=0;i<nCells.t;i++)
	{
		cells[i]=-1;
	}

	if(nI==0 && index==NULL)
	{
		//calculate each point's cell hash
		for(int i=0;i<sL;i++)
		{
			sorted[i]=l[i];
			sorted[i].type=-1;
			x=int(sorted[i].x/cellSize.x);
			y=int(sorted[i].y/cellSize.y);
			z=int(sorted[i].z/cellSize.z);
			indices[i]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[i]=i;//index of the particle
		}

		//fill cells
		for(int i=0;i<sL;i++)
		{
			//push operation
			current=indices[i];
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
			sorted[j]=l[j];
			sorted[j].type=-1;
			x=int(sorted[j].x/cellSize.x);
			y=int(sorted[j].y/cellSize.y);
			z=int(sorted[j].z/cellSize.z);
			indices[j]=x+y*nCells.x+z*nCells.x*nCells.y;//index of the cell
			linkedList[j]=j;//index of the particle
		}

		//fill cells
		for(int i=0;i<nI;i++)
		{
			int j=index[i];
			//push operation
			current=indices[j];
			buf=cells[current];
			cells[current]=linkedList[j];
			linkedList[j]=buf;
		}
	}
	
	//we will just traverse all boxes for now
	/*
	for(int i=0,nFullCells=0;i<nCells.t;i++)
	{
		if(cells[i]!=-1)
		{
			fullCells[nFullCells]=cells[i];
			nFullCells++;
		}
	}
	*/
}

template <typename T, typename S>
int Neighbors<T,S>::compute(bool compare(S, S, T), T compareVal)
{
	nDomains=0;
	T cutoffSquared=cutoff*cutoff;
	//go through each point
	for(int i=0;i<sL;i++)
	{
		if(sorted[i].type==-1)
		{
			//all variables that are private to this cell
			int x,y,z;
			T dx,dy,dz,dr;
			
			//flag it
			sorted[i].type=nDomains;
			
			//group nearby points within cutoff
			stack[0]=i;
			for(int j=0;j>-1;)
			{
				//get current point and cell it is in
				int current=stack[j];
				
				x=int(sorted[current].x/cellSize.x);
				y=int(sorted[current].y/cellSize.y);
				z=int(sorted[current].z/cellSize.z);
				
				int cCell=x+y*nCells.x+z*nCells.x*nCells.y;
				
				//pop current point
				j--;//don't put this at the end of the loop!
				
				//look around neigboring cells
				for(int tx=-1;tx<2;tx++)
				{
					for(int ty=-1;ty<2;ty++)
					{
						for(int tz=-1;tz<2;tz++)
						{
							//next cell index
							int k=((tx+x+nCells.x)%nCells.x)+
							(((ty+y+nCells.y)%nCells.y)*nCells.x)+
							(((tz+z+nCells.z)%nCells.z)*nCells.x*nCells.y);
							
							if(cells[k]!=-1)
							{
								for(int m=cells[k];m!=-1;m=linkedList[m])
								{
									dx=l[m].x-l[current].x;
									dy=l[m].y-l[current].y;
									dz=l[m].z-l[current].z;
									dr=dx*dx+dy*dy+dz*dz;
									if(dr<cutoffSquared && compare(conf[m],conf[current]))
									{
										j++;
										if(j>sL)
										{
											std::cout << "Stack overflow in neighbor computation!\n";
											throw 0;
										}
										stack[j]=m;
										
										//flag it
										sorted[m].type=nDomains;
									}
								}
								
							}
						}
					}
				}
			}
			nDomains++;
		}
	}
	
	return nDomains;
}

template <typename T>
class Kinetic {
	public:
		Kinetic(){};
		Kinetic(threeVector<T> *velocity, int nParticles);
		void initialize(threeVector<T> *velocity, int nParticles);
		T compute();
		T output(T time, char *name);
	private:
		threeVector<T> *v;
		int nP;
};

template <typename T>
Kinetic<T>::Kinetic(threeVector<T> *velocity, int nParticles)
{
	this->initialize(velocity, nParticles);
}

template <typename T>
void Kinetic<T>::initialize(threeVector<T> *velocity, int nParticles)
{
	v=velocity;
	nP=nParticles;
}

template <typename T>
T Kinetic<T>::compute()
{
	T energy=0;
	#pragma omp parallel for reduction(+:energy)
	for(int i=0;i<nP;i++)
		energy+=((v[i].x*v[i].x+v[i].y*v[i].y+v[i].z*v[i].z)/2.0);
	return energy;
}

template <typename T>
T Kinetic<T>::output(T time, char *name)
{
	T energy=compute();//this really should be done for all classes similar to this one
	
	std::fstream dataFile;
	std::string buf("kinetic_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << time << '\t' << energy << std::endl;
	dataFile.close();
	
	return energy;
}

/*
//recursive template?
//T needs to provide:
// flagged(int index)
//

//simple recursion that doesn't use cache memory.
//T requires:
// operator ()
template <typename T, T compute(T,T)>
void breadthRecursion(T *input, int current, int length)
{
	std::vector<int> stack;
	stack.push_back(compute(input[current], input[0]));
	while(!stack.empty())
	{
		while(!stack.back().end())
		{
			stack.push_back(stack.back().compute());
		}
	}
}

template <typename T>
void assignLocality(T input, std::vector<T> &output, int nInput)
{
	std::vector<int> stack;
	///This might be the best stack implementation I've made. All others need to use this.
	for(int i=0;i<nInput;i++)
	{
		//if it hasn't been flagged
		if(!input.flagged(i))
		{
			std::vector<int> bleb;
			//push it on stack, flag it, and begin neighbor search
			stack.clear();
			stack.push_back(i);
			input.flag(i);
			while(!stack.empty())
			{
				//pop it off stack
				int current=stack.back();
				stack.pop_back();
				output.push_back(current);
				
				//find nearest neighbors
				threeVector<double> minImg(0);
				for(int k=neighbors.queryIndex(current,minImg);k!=-1;k=neighbors.queryIndex(current,minImg))
				{
					//only check range if it isn't excluded
					if(!input.flagged(k))
					{
						int iIndex=index[current];
						int kIndex=index[k];
						threeVector<double> d;
						d.x=p[iIndex].x-p[kIndex].x+minImg.x;
						d.y=p[iIndex].y-p[kIndex].y+minImg.y;
						d.z=p[iIndex].z-p[kIndex].z+minImg.z;
						
						double dr=d.x*d.x+d.y*d.y+d.z*d.z;
						
						//is it within the range?
						if(dr<cutoff*cutoff)
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
	
}
*/
/*
///Bond length and bend length, use in conjunction with bondedMolecules bond functions
template <typename T>
inline T bondLength(position<T> &p1, position<T> &p2, T *uC, threeVector<T> s)
{
	double dx=p1.x-p2.x;
	double dy=p1.y-p2.y;
	double dz=p1.z-p2.z;
	if(dx>=s.x/2.0) dx-=s.x;
	if(dx<=-s.x/2.0) dx+=s.x;
	if(dy>=s.y/2.0) dy-=s.y;
	if(dy<=-s.y/2.0) dy+=s.y;
	if(dz>=s.z/2.0) dz-=s.z;
	if(dz<=-s.z/2.0) dz+=s.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//This isn't logical, but whatever
template <typename T>
inline T bendLength(position<T> &p1,position<T> &p2,position<T> &p3,
		   threeVector<T> &a1,threeVector<T> &a2,threeVector<T> &a3,T *uC, threeVector<T> s)
{
	return 0;
}
*/
//end of MD_DATACOLLECTION
#endif