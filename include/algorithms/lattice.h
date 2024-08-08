#include "functions.h"
#include "bondedMolecules.h"

#ifndef MD_LATTICE
#define MD_LATTICE


//T is a simple data type restriction (int, float, double...).
//U is a class with a check function that takes a position type and outputs a position type.
// The output is either type=-1 for a failure or type=something for success and save.
template <typename T, typename U>
class Lattice {
	public:
		//constraint is an object that 
		Lattice(position<T> * &p, molecules<T> * &m);//normal version, just outputs the lattice
		Lattice(position<T> * &p, molecules<T> * &m, U *constraint);//alternate version
		~Lattice();
		void fcc(threeVector<T> size, threeVector<int> nNodes);//makes an fcc lattice
		void bcc(threeVector<T> size, threeVector<int> nNodes);//makes a bcc lattice
		void cubic(threeVector<T> size, threeVector<int> nNodes);
		int sizeP(){return nP;};
		int sizeM(){return nM;};
	private:
		position<T> **p;
		molecules<T> **m;
		U *constraint;
		bool allocated;
		
		//variables required by this class for internal use
};

template <typename T, typename U>
Lattice<T,U>::Lattice(position<T> * &p, U *constraint)
{
	this->p=&p;
	this->m=&m;
	this->constraint=constraint;
	bool allocated=false;
}

//the destructor should be defined, but usually used to deallocate memory
template <typename T, typename U>
Lattice<T,U>::~Lattice()
{
	if(allocated)
	{
		delete (*p);
		delete (*m);
	}
}

//for rebuilding run time data structures related to execution
//this isn't always needed, but is a good idea if you don't need to build all the time
//e.g. when one data structure lasts 20 time steps before needing to update
template <typename T, typename U>
void Lattice<T,U>::fcc(threeVector<T> size, threeVector<int> nNodes)
{

}

//This should just compute and output a value, but you could pass a pointer and find other values
template <typename T, typename U>
void Lattice<T,U>::bcc(threeVector<T> size, threeVector<int> nNodes)
{

}

void Lattice<T,U>::cubic(threeVector<T> size, threeVector<int> nNodes)
{
	if(allocated)
	{
		delete (*p);
		delete (*m);
	}
	threeVector<T> length;
	length.x=size.x/(T)nNodes.x;
	length.y=size.y/(T)nNodes.y;
	length.z=size.z/(T)nNodes.z;
	
	nP=nNodes.x*nNodes.y*nNodes.z;
	nM=nNodes.x*nNodes.y*nNodes.z*9;
	(*p)=new position<T>[nP];
	(*m)=new molecule<T>[nM];
	int currentNP=0;
	position<T> point;
	for(int i=0;i<nNodes.x;i++)
	{
		for(int j=0;j<nNodes.y;j++)
		{
			for(int k=0;k<nNodes.z;k++)
			{
				point.type=0;
				point.x=(T)i*length.x;
				point.y=(T)j*length.y;
				point.z=(T)k*length.z;
				position<T> test=constraint.check(point);
				if(test.type!=-1)
					(*p)[currentNP++]=test;
			}
		}
	}
	nP=currentNP;
	
	currentNM=0;
	
	//BOND types
	for(int i=0;i<nP;i++)
	{
		for(int j=i+1;j<nP;j++)
		{
			threeVector<T> ds;
			ds.x=p[i].x-p[j].x;
			ds.y=p[i].y-p[j].y;
			ds.z=p[i].z-p[j].z;
			//I'm assuming nothing
			if(sqrt(ds.x*ds.x+ds.y*ds.y+ds.z*ds.z)<constraint.bondLength())
			{
				(*m)[currentNM].bond=new genericVector<int,4>[1];
				(*m)[currentNM].type=BOND;
				(*m)[currentNM].bond[0]=i;
				(*m)[currentNM].bond[1]=j;
				currentNM++;
			}
		}
	}
	
	nM=currentNM;
	
	//BEND types
	for(int i=0;i<nM;i++)
	{
		for(j=i+1;j<nM;j++)
		{
			if((*m)[i].bond[1]=(*m)[j].bond[0])
			{
				(*m)[currentNM].bond=new genericVector<int,4>[1];
				(*m)[currentNM].type=BEND;
				(*m)[currentNM].bond[0]=(*m)[i].bond[0];
				(*m)[currentNM].bond[1]=(*m)[i].bond[1];
				(*m)[currentNM].bond[2]=(*m)[i].bond[1];
			}
		}
	}
			
}

#endif