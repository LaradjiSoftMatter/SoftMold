//This file contains data measurements for lipid membranes
//split this file up
#include "../algorithms/functions.h"

#ifndef LIPIDS
#define LIPIDS

template <typename T>
class Domains {
	public:
		Domains(position<T> *particles, T *FConstants, T *UConstants, T *FBondConstants, T *UBondConstants,
			T *FBendConstants, T *UBendConstants, T *FDihedralConstants, T *UDihedralConstants, 
			T *FTorsionConstants, T *UTorsionConstants, molecule<T> *mol, int nParticles, int nTypes,
			int nMolecules, threeVector<T> size);
		
		~Domains();
		
		T compute();//average size in system units
		T output(T time, char *name);//returns average size, outputs average size, nDomains
		
		int *domains;//start indices of mCom domains
		int nDomains;//number of domains
		position<T> *mCom;//center of mass for all molecules
	private:
		position<T> *p;//particle positions
		threeVector<T> *a;//accelerations
		molecule<T> *m;
		
		//constants
		T *cU;
		T *cF;
		T *cBondU;//bond potential constants
		T *cBondF;//bond force constants
		T *cBendU;//bend potential constants
		T *cBendF;//bend force constants
		T *cDihedralU;//dihedral potential constants
		T *cDihedralF;//dihedral force constants
		T *cTorsionU;//torsion potential constants
		T *cTorsionF;//torsion force constants

		//number of each constant
		int nP;//is this even needed?
		int nT;
		threeVector<T> s;//size of system
		int nM;//number of molecules
};

template <typename T>
Domains<T>::Domains(position<T> *particles, T *FConstants, T *UConstants, T *FBondConstants, T *UBondConstants,
		 T *FBendConstants, T *UBendConstants, T *FDihedralConstants, T *UDihedralConstants, 
		 T *FTorsionConstants, T *UTorsionConstants, molecule<T> *mol, int nParticles, int nTypes,
		 int nMolecules, threeVector<T> size)
{
	
}

template <typename T>
class Order {
	public:
		Order(position<T> *particles, molecule<T> mol, int nParticles, int nTypes, threeVector<T> size);
		
		~Order();
		
		//compute an order paramter against some normal
		T compute(threeVector<T> normal){return 0;};
		T output(T time, char *name, threeVector<T> normal){return 0;};
		
		//compute a local order paramter within some radius
		T compute(T radius){return 0;};
		T output(T time, char *name, T radius){return 0;};
		
		//compute the straightness of the chains as an order parameter
		T compute(){return 0;};//returns order parameter
		T output(){return 0;};//returns order parameter, outputs average size, nDomains, order parameter, 
		
		//compute the order parameter as a function of inside, outside, and average of a sphere
		threeVector<T> compute(threeVector<T> com, T radius);
		threeVector<T> output(T time, char *name, threeVector<T> com, T radius);
		
		//this should be done, just not now
		int *domains;//start indices of mCom domains
		int nDomains;//number of domains
		threeVector<T> *mCom;//center of mass for all molecules
		threeVector<T> *orientation;
	private:
		position<T> *p;//particle positions
		molecule<T> m;//one molecule list

		//number of each constant
		int nP;//is this even needed?
		int nT;
		threeVector<T> s;//size of system
};

template <typename T>
Order<T>::Order(position<T> *particles, molecule<T> mol, int nParticles, int nTypes, threeVector<T> size)
{
	this->p=particles;
	this->m=mol;
	if(m.type!=CHAIN || m.nBonded!=1)
	{
		std::cout << "Molecule not type CHAIN!\n";
		throw 0;
	}
	if(m.nBonded!=1)
	{
		std::cout << "More than one group in molecule!\n";
		throw 0;
	}
	this->nP=nParticles;
	this->nT=nTypes;
	this->s=size;
	this->mCom=new threeVector<T>[m.bonded[0].group[NCHAINS]];
	this->orientation=new threeVector<T>[m.bonded[0].group[NCHAINS]];
}

template <typename T>
Order<T>::~Order()
{
	
}

template <typename T>
threeVector<T> Order<T>::compute(threeVector<T> com, T radius)
{
	threeVector<T> order;
	order.x=0;//inner pointing
	order.y=0;//outer
	order.z=0;//average
	threeVector<int> nDir;
	nDir.x=0;
	nDir.y=0;
	nDir.z=m.bonded[0].group[NCHAINS];
	
	for(int i=0;i<m.bonded[0].group[NCHAINS];i++)
	{
		threeVector<T> vec;
		T currentOrder=0;
		int head=m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i;
		int tail=m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*(i+1)-1;
		bool swaped=false;
		
		if(p[head].type>p[tail].type)
		{
			int buf=head;
			head=tail;
			tail=buf;
			swaped=true;
		}
		
		orientation[i].x=p[head].x-p[tail].x;
		orientation[i].y=p[head].y-p[tail].y;
		orientation[i].z=p[head].z-p[tail].z;
		
		orientation[i]=unitVector<T>(orientation[i]);//this is our "normal" vector
		
		mCom[i].x=0;
		mCom[i].y=0;
		mCom[i].z=0;
		
		for(int j=m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i;
		    j<m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*(i+1);
		    j++)
		{
			mCom[i].x+=p[j].x;
			mCom[i].y+=p[j].y;
			mCom[i].z+=p[j].z;
			
			if(swaped && j<m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*(i+1)-1)
			{
				vec.x=p[j+1].x-p[j].x;
				vec.y=p[j+1].y-p[j].y;
				vec.z=p[j+1].z-p[j].z;
				vec=unitVector<T>(vec);
				currentOrder+=dotProduct<T>(vec,orientation[i]);
			}
			
			if(!swaped && j<m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*(i+1)-1)
			{
				vec.x=p[j].x-p[j+1].x;
				vec.y=p[j].y-p[j+1].y;
				vec.z=p[j].z-p[j+1].z;
				vec=unitVector<T>(vec);
				currentOrder+=dotProduct<T>(vec,orientation[i]);
			}
		}
		
		currentOrder/=(T)m.bonded[0].group[LENGTH]-1;
		
		threeVector<T> first,second;
		first.x=p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i].x-p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+1].x;
		first.y=p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i].y-p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+1].y;
		first.z=p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i].z-p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+1].z;
		
		second.x=p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+1].x-p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+2].x;
		second.y=p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+1].y-p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+2].y;
		second.z=p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+1].z-p[m.bonded[0].group[START]+m.bonded[0].group[LENGTH]*i+2].z;
		
		first=unitVector<T>(first);
		second=unitVector<T>(second);
		
		currentOrder=dotProduct<T>(first,second);
		currentOrder*=currentOrder;
		
		mCom[i].x/=(T)m.bonded[0].group[LENGTH];
		mCom[i].y/=(T)m.bonded[0].group[LENGTH];
		mCom[i].z/=(T)m.bonded[0].group[LENGTH];
		
		//mCom[i]=unitVector<T>(mCom[i]);
		
		vec.x=mCom[i].x-com.x;
		vec.y=mCom[i].y-com.y;
		vec.z=mCom[i].z-com.z;
		
		T direction=dotProduct<T>(orientation[i], vec);
		
		if(direction<=0)
		{
			nDir.x++;
			order.x+=currentOrder;
		}
		else
		{
			nDir.y++;
			order.y+=currentOrder;
		}
	}
	
	order.z=order.x+order.y;
	order.x/=(T)nDir.x;
	order.y/=(T)nDir.y;
	order.z/=(T)nDir.z;
	
	return order;
}

template <typename T>
threeVector<T> Order<T>::output(T time, char *name, threeVector<T> com, T radius)
{
	threeVector<T> order=compute(com,radius);
	
	std::fstream dataFile;
	std::string buf("order_");
	buf+=name;
	buf+=".dat";
	dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
	dataFile << time << '\t' << order.x << '\t' << order.y << '\t' << order.z << std::endl;//potentially more values here, maybe even more files
	dataFile.close();
	
	return order;
}
/*
template <typename T,T space(position <T> ) >
class BilayerDistribution
{
	public:
		BilayerDistribution(position<T> *particles, molecule<T> mol, int nParticles, int nTypes, threeVector<T> size, );
		~BilayerDistribution();
		void build();
	private:
		
}
*/
#endif
