//algorithms for integrating motions of a particle system
//rename this verlet.h or something
#include "functions.h"

#ifndef MD_VERLET
#define MD_VERLET

template <typename T>
class Verlet {
	public:
		Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		       int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap);
		Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		       int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap, int *index, int nIndices);
		Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		       int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap, position<T> *absolutePosition);
		Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		       int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap, position<T> *absolutePosition,
		       int *index, int nIndices);
		~Verlet();
		void first();//first part of integration
		void second();//after acceleration updates, second part
		void resize(threeVector<T> size){s=size;};
		
	private:
		//position, acceleration, velocity
		position<T> *p;
		threeVector<T> *a;
		threeVector<T> *v;
#ifdef VERLET_HALF_CARRY		
		threeVector<T> *vHalf;
		threeVector<T> *aHalf;
#endif		
		int nP;//number of particles
		threeVector<T> s;//size
		threeVector<bool> w;//for wrapping simulation in a box
		T dt;//time step
		//for adjusting when integration occurs, currently it is just at dt/2.0, implement in future
		T dtFirst;
		T dtSecond;
		
		bool diffusion;
		position<T> *aP;//an array dedicated to tracking the diffusion of p particles, absolute position
		
		int *index;
		int nI;
};

template <typename T>
Verlet<T>::Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		  int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap)
{
	this->p=particles;
	this->a=accelerations;
	this->v=velocities;
	this->nP=nParticles;
	this->s=size;
	this->dt=timeStep;
	this->w=wrap;
	this->aP=NULL;
	this->diffusion=false;
	this->index=NULL;
	this->nI=0;
#ifdef VERLET_HALF_CARRY	
	vHalf=new threeVector<T>[nParticles];
	aHalf=new threeVector<T>[nParticles];
#endif
}

template <typename T>
Verlet<T>::Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		  int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap, int *index, int nIndices)
{
	this->p=particles;
	this->a=accelerations;
	this->v=velocities;
	this->nP=nParticles;
	this->s=size;
	this->dt=timeStep;
	this->w=wrap;
	this->aP=NULL;
	this->diffusion=false;
	this->index=index;
	this->nI=nIndices;
#ifdef VERLET_HALF_CARRY	
	vHalf=new threeVector<T>[nParticles];
	aHalf=new threeVector<T>[nParticles];
#endif
}

template <typename T>
Verlet<T>::Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		  int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap, position<T> *absolutePosition)
{
	this->p=particles;
	this->a=accelerations;
	this->v=velocities;
	this->nP=nParticles;
	this->s=size;
	this->dt=timeStep;
	this->w=wrap;
	this->aP=absolutePosition;
	this->diffusion=true;
	this->index=NULL;
	this->nI=0;
#ifdef VERLET_HALF_CARRY	
	vHalf=new threeVector<T>[nParticles];
	aHalf=new threeVector<T>[nParticles];
#endif
}

template <typename T>
Verlet<T>::Verlet(position<T> *particles, threeVector<T> *accelerations, threeVector<T> *velocities,
		  int nParticles, threeVector<T> size, T timeStep, threeVector<bool> wrap, position<T> *absolutePosition,
		  int *index, int nIndices)
{
	this->p=particles;
	this->a=accelerations;
	this->v=velocities;
	this->nP=nParticles;
	this->s=size;
	this->dt=timeStep;
	this->w=wrap;
	this->aP=absolutePosition;
	this->diffusion=true;
	this->index=index;
	this->nI=nIndices;
#ifdef VERLET_HALF_CARRY	
	vHalf=new threeVector<T>[nParticles];
	aHalf=new threeVector<T>[nParticles];
#endif
}

template <typename T>
Verlet<T>::~Verlet()
{
#ifdef VERLET_HALF_CARRY
	delete vHalf,aHalf;
#endif
}

template <typename T>
void Verlet<T>::first()
{
	T halfDt=0.5*dt;
#ifdef VERLET_HALF_CARRY
	T halfDtSquared=0.5*dt*dt;
	if(nI==0 && index==NULL)
	{
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			//increment position
			if(p[i].type!=0)
			{
				p[i].x+=v[i].x*dt+halfDtSquared*a[i].x;
				p[i].y+=v[i].y*dt+halfDtSquared*a[i].y;
				p[i].z+=v[i].z*dt+halfDtSquared*a[i].z;
			}
		}
		
		//check if diffusion is active...
		if(aP!=NULL)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			//increment position
			if(aP[i].type!=0)
			{
				aP[i].x+=v[i].x*dt+halfDtSquared*a[i].x;
				aP[i].y+=v[i].y*dt+halfDtSquared*a[i].y;
				aP[i].z+=v[i].z*dt+halfDtSquared*a[i].z;
			}
		}
		
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			//increment velocity
			if(p[i].type!=0)
			{
				vHalf[i]=v[i];
				aHalf[i]=a[i];
				v[i].x+=(a[i].x*halfDt);
				v[i].y+=(a[i].y*halfDt);
				v[i].z+=(a[i].z*halfDt);
			}
		}
		
		//check if wrapping of component is active
		if(w.x)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].x>=s.x) p[i].x-=s.x;
			if(p[i].x<0) p[i].x+=s.x;
		}
		
		if(w.y)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].y>=s.y) p[i].y-=s.y;
			if(p[i].y<0) p[i].y+=s.y;
		}
		
		if(w.z)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].z>=s.z) p[i].z-=s.z;
			if(p[i].z<0) p[i].z+=s.z;
		}
	}
	else
	{
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			//increment position
			if(p[i].type!=0)
			{
				p[i].x+=v[i].x*dt+halfDtSquared*a[i].x;
				p[i].y+=v[i].y*dt+halfDtSquared*a[i].y;
				p[i].z+=v[i].z*dt+halfDtSquared*a[i].z;
			}
		}
		
		//check if diffusion is active...
		if(aP!=NULL)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			//increment position
			if(aP[i].type!=0)
			{
				aP[i].x+=v[i].x*dt+halfDtSquared*a[i].x;
				aP[i].y+=v[i].y*dt+halfDtSquared*a[i].y;
				aP[i].z+=v[i].z*dt+halfDtSquared*a[i].z;
			}
		}
		
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			//increment velocity
			if(p[i].type!=0)
			{
				vHalf[i]=v[i];
				aHalf[i]=a[i];
				v[i].x+=(a[i].x*halfDt);
				v[i].y+=(a[i].y*halfDt);
				v[i].z+=(a[i].z*halfDt);
			}
		}
		
		//check if wrapping of component is active
		if(w.x)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].x>=s.x) p[i].x-=s.x;
			if(p[i].x<0) p[i].x+=s.x;
		}
		
		if(w.y)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].y>=s.y) p[i].y-=s.y;
			if(p[i].y<0) p[i].y+=s.y;
		}
		
		if(w.z)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].z>=s.z) p[i].z-=s.z;
			if(p[i].z<0) p[i].z+=s.z;
		}
	}	
#else
//No verlet half carry
	if(nI==0 && index==NULL)
	{
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			//increment velocity
			if(p[i].type!=0)
			{
				v[i].x+=(a[i].x*halfDt);
				v[i].y+=(a[i].y*halfDt);
				v[i].z+=(a[i].z*halfDt);
			}
		}
		
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			//increment position
			if(p[i].type!=0)
			{
				p[i].x+=v[i].x*dt;
				p[i].y+=v[i].y*dt;
				p[i].z+=v[i].z*dt;
			}
		}
		
		//check if diffusion is active...
		if(aP!=NULL)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			//increment position
			if(aP[i].type!=0)
			{
				aP[i].x+=v[i].x*dt;
				aP[i].y+=v[i].y*dt;
				aP[i].z+=v[i].z*dt;
			}
		}
		
		
		//check if wrapping of component is active
		if(w.x)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
		
			if(p[i].x>s.x) p[i].x-=s.x;
			if(p[i].x<0) p[i].x+=s.x;
		}
		
		if(w.y)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].y>s.y) p[i].y-=s.y;
			if(p[i].y<0) p[i].y+=s.y;
		}
		
		if(w.z)
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].z>s.z) p[i].z-=s.z;
			if(p[i].z<0) p[i].z+=s.z;
		}
	}
	else
	{
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			//increment velocity
			if(p[i].type!=0)
			{
				v[i].x+=(a[i].x*halfDt);
				v[i].y+=(a[i].y*halfDt);
				v[i].z+=(a[i].z*halfDt);
			}
		}
		
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			//increment position
			if(p[i].type!=0)
			{
				p[i].x+=v[i].x*dt;
				p[i].y+=v[i].y*dt;
				p[i].z+=v[i].z*dt;
			}
		}
		
		//check if diffusion is active...
		if(aP!=NULL)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			//increment position
			if(aP[i].type!=0)
			{
				aP[i].x+=v[i].x*dt;
				aP[i].y+=v[i].y*dt;
				aP[i].z+=v[i].z*dt;
			}
		}
		
		//check if wrapping of component is active
		if(w.x)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].x>s.x) p[i].x-=s.x;
			if(p[i].x<0) p[i].x+=s.x;
		}
		
		if(w.y)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].y>s.y) p[i].y-=s.y;
			if(p[i].y<0) p[i].y+=s.y;
		}
		
		if(w.z)
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].z>s.z) p[i].z-=s.z;
			if(p[i].z<0) p[i].z+=s.z;
		}
	}
#endif
}

template <typename T>
void Verlet<T>::second()
{
	T halfDt=dt*0.5;//this is not what this is... it will be replaced
#ifdef VERLET_HALF_CARRY
	if(index==NULL && nI==0)
	{
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].type!=0)
			{
				v[i].x=vHalf[i].x+(a[i].x+aHalf[i].x)*halfDt;
				v[i].y=vHalf[i].y+(a[i].y+aHalf[i].y)*halfDt;
				v[i].z=vHalf[i].z+(a[i].z+aHalf[i].z)*halfDt;
			}
		}
	}
	else
	{
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].type!=0)
			{
				v[i].x=vHalf[i].x+(a[i].x+aHalf[i].x)*halfDt;
				v[i].y=vHalf[i].y+(a[i].y+aHalf[i].y)*halfDt;
				v[i].z=vHalf[i].z+(a[i].z+aHalf[i].z)*halfDt;
			}
		}
	}
#else
//No verlet half carry
	if(index==NULL && nI==0)
	{
		#pragma omp parallel for
		for(int i=0;i<nP;i++)
		{
			if(p[i].type!=0)
			{
				v[i].x+=(a[i].x*halfDt);
				v[i].y+=(a[i].y*halfDt);
				v[i].z+=(a[i].z*halfDt);
			}
		}
	}
	else
	{
		#pragma omp parallel for
		for(int j=0;j<nI;j++)
		{
			int i=index[j];
			if(p[i].type!=0)
			{
				v[i].x+=(a[i].x*halfDt);
				v[i].y+=(a[i].y*halfDt);
				v[i].z+=(a[i].z*halfDt);
			}
		}
	}
#endif
}

//end of MD_VERLET
#endif
