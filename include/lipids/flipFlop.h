/*
 * \brief Algorithm for determining flip flop of lipids. Also does multiple components.
 *
 *
 */

#ifndef FLIPFLOP_MD
#define FLIPFLOP_MD

#include "../algorithms/molecules.h"

template <typename T>
class flipFlop {
	public:
		flipFlop();
		flipFlop(position<T> *particles, int nParticles, threeVector<double> size, 
		         T cutoff, molecule<T, fourVector<int> > * mol, int nMolecules
		         int *vecBeginTypes, int *vecEndTypes, int nVec);
		
		void initialize(position<T> *particles, int nParticles, threeVector<double> size, 
		         T cutoff, molecule<T, fourVector<int> > * mol, int nMolecules
		         int *vecBeginTypes, int *vecEndTypes, int nVec);
		
		~flipFlop();
		
		void build();
		void compute();
		
	private:
		//constants
		T rc;//cutoff
		position<T> *p;//particles
		int nP;//number of particles
		threeVector<T> s;//size
		molecule<T, fourVector<int> > *m;//molecules
		int nM;//number of molecules
		int *vBT;//vector begin types
		int *vET;//vector end types
		int nV;//number of vector types
		
		//values
		//inner surface is a surface where all vectors converge
		//outer surface is a surface where all vectors diverge
		//for open surfaces inner is net negative and outer is net positive
		bool *flip;
		
		int nS;
		
		int *surface;//index which surface a vector rests on
};

template <typename T>
flipFlop<T>::flipFlop()
{
	p=NULL;
	m=NULL;
	vBT=NULL;
	vET=NULL;
	nV=0;
	nM=0;
	nP=0;
	s.x=0;
	s.y=0;
	s.z=0;
	rc=0;
	flip=NULL;
	surface=NULL;
	nS=0;
}

template <typename T>
flipFlop<T>::flipFlop(position<T> *particles, int nParticles, threeVector<double> size, 
		         T cutoff, molecule<T, fourVector<int> > * mol, int nMolecules
		         int *vecBeginTypes, int *vecEndTypes, int nVec)
{
	this->initialize(particles, nParticles, size, cutoff, mol, nMolecules, vecBeginTypes, vecEndTypes, nVec);
}

template <typename T>
void flipFlop<T>::initialize(position<T> *particles, int nParticles, threeVector<double> size, 
		         T cutoff, molecule<T, fourVector<int> > * mol, int nMolecules
		         int *vecBeginTypes, int *vecEndTypes, int nVec);
{
	this->p=particles;
	this->nP=nParticles;
	this->s=size;
	this->rc=cutoff;
	this->vBT=vecBeginTypes;
	this->vET=vecEndTypes;
	this->nV=nVec;
	this->nM=nMolecules;
	this->m=mol;
	flip=new bool[nP];
	surface=new int[nP];


}

template <typename T>
flipFlop<T>::~flipFlop()
{
	if(surface!=NULL)
		delete surface;
	if(flip!=NULL)
		delete flip;
}



#endif
