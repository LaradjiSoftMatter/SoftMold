//algorithms for thermostat functions
//rename langevin.h or something

#include "functions.h"

#ifndef MD_LANGEVIN
#define MD_LANGEVIN

template <typename T>
class Langevin {
	public:
		Langevin()
		{
			this->a=NULL;
			this->v=NULL;
			this->p=NULL;
			this->nP=0;
			this->g=0;
			this->gT=NULL;
			this->dT=0;
			this->index=NULL;
			this->nI=0;
			this->randNum=NULL;
			this->sT=NULL;
		};
		Langevin(threeVector<T> *acceleration, threeVector<T> *velocity, int nParticles, T gamma,
			 T timeStep, int seed)
		{
			this->initialize(acceleration, velocity, nParticles, gamma, timeStep, seed);
		};
		Langevin(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions,
			 int nParticles, T gamma, T timeStep, int seed)
		{
			this->initialize(acceleration, velocity, positions, nParticles, gamma, timeStep, seed);
		};
		Langevin(threeVector<T> *acceleration, threeVector<T> *velocity, int nParticles, T gamma,
			 T timeStep, int seed, int *index, int nIndices)
		{
			this->initialize(acceleration, velocity, nParticles, gamma, timeStep, seed, index, nIndices);
		};
		Langevin(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, int nParticles, T* gammaType,
			 T timeStep, int seed, int nTypes)
		{
			this->initialize(acceleration, velocity, positions, nParticles, gammaType, timeStep, seed);
		};
		Langevin(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, int nParticles, T* gammaType,
			 T timeStep, int seed, int nTypes, int *index, int nIndices)
		{
			this->initialize(acceleration, velocity, positions, nParticles, gammaType, timeStep, seed, index, nIndices);
		};
		
		void initialize(threeVector<T> *acceleration, threeVector<T> *velocity, int nParticles, T gamma,
			 T timeStep, int seed);
		void initialize(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, 
			 int nParticles, T gamma, T timeStep, int seed);
		void initialize(threeVector<T> *acceleration, threeVector<T> *velocity, int nParticles, T gamma,
			 T timeStep, int seed, int *index, int nIndices);
		void initialize(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, int nParticles, T* gammaType,
			 T timeStep, int seed, int nTypes);
		void initialize(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, int nParticles, T* gammaType,
			 T timeStep, int seed, int nTypes, int *index, int nIndices);
		~Langevin();
		void compute(T temperature);//ramp your own temperature
	private:
		//acceleration, velocity, positions
		threeVector<T> *a;
		threeVector<T> *v;
		position<T> *p;
		
		int nP;//number of particles
		T g;//gamma
		T *gT;//gamma by type
		T *sT;//sigma by type
		int nT;//number of types
		T dT;//time step
		MTRand **randNum;
		//these are flagged with NULL
		int *index;
		int nI;
		int nThreads;
};

template <typename T>
void Langevin<T>::initialize(threeVector<T> *acceleration, threeVector<T> *velocity, int nParticles, T gamma,
		   T timeStep, int seed)
{
	this->a=acceleration;
	this->v=velocity;
	this->p=NULL;
	this->nP=nParticles;
	this->g=gamma;
	this->gT=NULL;
	this->dT=timeStep;
	this->index=NULL;
	this->nI=0;
	this->sT=NULL;
	this->nT=0;
	#ifdef _OPENMP
		nThreads=omp_get_max_threads();
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		nThreads=1;
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
}

template <typename T>
void Langevin<T>::initialize(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions,
		 int nParticles, T gamma, T timeStep, int seed)
{
	this->a=acceleration;
	this->v=velocity;
	this->p=positions;
	this->nP=nParticles;
	this->g=gamma;
	this->gT=NULL;
	this->dT=timeStep;
	this->index=NULL;
	this->nI=0;
	this->sT=NULL;
	this->nT=0;
	#ifdef _OPENMP
		nThreads=omp_get_max_threads();
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		nThreads=1;
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
}

template <typename T>
void Langevin<T>::initialize(threeVector<T> *acceleration, threeVector<T> *velocity, int nParticles, T gamma,
		   T timeStep, int seed, int *index, int nIndices)
{
	this->a=acceleration;
	this->v=velocity;
	this->p=NULL;
	this->nP=nParticles;
	this->gT=NULL;
	this->g=gamma;
	this->dT=timeStep;
	this->index=index;
	this->nI=nIndices;
	this->sT=NULL;
	this->nT=0;
	#ifdef _OPENMP
		nThreads=omp_get_max_threads();
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		nThreads=1;
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
}

template <typename T>
void Langevin<T>::initialize(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, int nParticles, T* gammaType,
		   T timeStep, int seed, int nTypes)
{
	this->a=acceleration;
	this->v=velocity;
	this->p=positions;
	this->nP=nParticles;
	this->gT=gammaType;
	this->g=0;
	this->dT=timeStep;
	this->index=NULL;
	this->nI=0;
	this->sT=NULL;
	this->nT=nTypes;
	#ifdef _OPENMP
		nThreads=omp_get_max_threads();
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		nThreads=1;
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
}

template <typename T>
void Langevin<T>::initialize(threeVector<T> *acceleration, threeVector<T> *velocity, position<T> *positions, int nParticles, T* gammaType,
		   T timeStep, int seed, int nTypes, int *index, int nIndices)
{
	this->a=acceleration;
	this->v=velocity;
	this->p=positions;
	this->nP=nParticles;
	this->gT=gammaType;
	this->g=0;
	this->dT=timeStep;
	this->index=index;
	this->nI=nIndices;
	this->sT=NULL;
	this->nT=nTypes;
	#ifdef _OPENMP
		nThreads=omp_get_max_threads();
		this->randNum=new MTRand *[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++)
			randNum[i]=new MTRand(seed+i);
	#else
		nThreads=1;
		this->randNum=new MTRand *[1];
		randNum[0]=new MTRand(seed);
	#endif
}

template <typename T>
Langevin<T>::~Langevin()
{
	#ifdef _OPENMP
		for(int i=0;i<omp_get_max_threads();i++)
			delete randNum[i];
	#else
		delete randNum[0];
	#endif
	if(sT!=NULL)
		delete sT;
	delete randNum;
}

template <typename T>
void Langevin<T>::compute(T temperature)
{
	T sigma=sqrt((6.0*temperature*g)/dT);
	if(gT!=NULL && sT==NULL)
	{
		sT=new T[nT];
		for(int i=0;i<nT;i++)
			sT[i]=sqrt((6.0*temperature*gT[i])/dT);
	}
	else
	{
		sigma=sqrt((6.0*temperature*g)/dT);
	}
	if(nI==0 && index==NULL)
	{
		if(gT!=NULL)
		{
			#pragma omp parallel for
			for(int i=0;i<nP;i++)
			{
				T psi;
				int type=0;//p[i].type;
				//if a segmentation fault occurs, this section is probably not
				//responsible, check the allPairs header for related problems
				//(out of bounds particles etc...)
				#ifdef _OPENMP
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].x+=(-gT[type]*v[i].x+sT[type]*psi);
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].y+=(-gT[type]*v[i].y+sT[type]*psi);
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].z+=(-gT[type]*v[i].z+sT[type]*psi);
				
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].x+=(-g*v[i].x+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].y+=(-g*v[i].y+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].z+=(-g*v[i].z+sigma*psi);
				#else
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].x+=(-gT[type]*v[i].x+sT[type]*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].y+=(-gT[type]*v[i].y+sT[type]*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].z+=(-gT[type]*v[i].z+sT[type]*psi);
				#endif
			}
		}
		else
		{
			#pragma omp parallel for
			for(int i=0;i<nP;i++)
			{
				T psi;
				//if a segmentation fault occurs, this section is probably not
				//responsible, check the allPairs header for related problems
				//(out of bounds particles etc...)
				#ifdef _OPENMP
					/*threeVector<T> d;
					d.x=p[i].x-10.0;
					d.y=p[i].y-10.0;
					if(d.x*d.x+d.y*d.y<32.0)
					{
						T gAlt=g;
						T tAlt=2.4;
						T sigmaAlt=sqrt((6.0*tAlt*gAlt)/dT);
						psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
						a[i].x+=(-gAlt*v[i].x+sigmaAlt*psi);
						psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
						a[i].y+=(-gAlt*v[i].y+sigmaAlt*psi);
						psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
						a[i].z+=(-gAlt*v[i].z+sigmaAlt*psi);
					}
					else*/
					//{
						psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
						a[i].x+=(-g*v[i].x+sigma*psi);
						psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
						a[i].y+=(-g*v[i].y+sigma*psi);
						psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
						a[i].z+=(-g*v[i].z+sigma*psi);
					//}
					
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].x+=(-g*v[i].x+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].y+=(-g*v[i].y+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].z+=(-g*v[i].z+sigma*psi);
				#else
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].x+=(-g*v[i].x+sigma*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].y+=(-g*v[i].y+sigma*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].z+=(-g*v[i].z+sigma*psi);
				#endif
			}
		}
	}
	else
	{
		if(gT!=NULL)
		{
			#pragma omp parallel for
			for(int j=0;j<nI;j++)
			{
				int i=index[j];
				T psi;
				int type=p[i].type;
				//if a segmentation fault occurs, this section is probably not
				//responsible, check the allPairs header for related problems
				//(out of bounds particles etc...)
				#ifdef _OPENMP
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].x+=(-gT[type]*v[i].x+sT[type]*psi);
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].y+=(-gT[type]*v[i].y+sT[type]*psi);
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].z+=(-gT[type]*v[i].z+sT[type]*psi);
				
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].x+=(-g*v[i].x+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].y+=(-g*v[i].y+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].z+=(-g*v[i].z+sigma*psi);
				#else
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].x+=(-gT[type]*v[i].x+sT[type]*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].y+=(-gT[type]*v[i].y+sT[type]*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].z+=(-gT[type]*v[i].z+sT[type]*psi);
				#endif
			}
		}
		else
		{
			#pragma omp parallel for
			for(int j=0;j<nI;j++)
			{
				int i=index[j];
				T psi;
				//if a segmentation fault occurs, this section is probably not
				//responsible, check the allPairs header for related problems
				//(out of bounds particles etc...)
				#ifdef _OPENMP
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].x+=(-g*v[i].x+sigma*psi);
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].y+=(-g*v[i].y+sigma*psi);
					psi=2.0*randNum[omp_get_thread_num()]->rand53()-1.0;
					a[i].z+=(-g*v[i].z+sigma*psi);
				
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].x+=(-g*v[i].x+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].y+=(-g*v[i].y+sigma*psi);
					//psi=2.0*randNum[i%nThreads]->rand53()-1.0;
					//a[i].z+=(-g*v[i].z+sigma*psi);
				#else
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].x+=(-g*v[i].x+sigma*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].y+=(-g*v[i].y+sigma*psi);
					psi=2.0*randNum[0]->rand53()-1.0;
					a[i].z+=(-g*v[i].z+sigma*psi);
				#endif
			}
		}
	}
}

//end of MD_LANGEVIN
#endif
