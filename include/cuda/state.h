#include <cstring>
#include "dataTypes.h"
#ifndef MPD_STATE
#define MPD_STATE

namespace mpd {
	
	template <typename T>
	struct state {
		using value_type=T;
		state(position<T> *p, threeVector<T> *v, threeVector<T> *a, T *NBFc, T *NBPc,
		      int nP, int nT, T dT, T g, T t, threeVector<T> s, T *m):
			p_h(p),v_h(v),a_h(a), NBFconstants_h(NBFc),NBPconstants_h(NBPc),
			nParticles(nP),nTypes(nT),deltaT(dT),gamma(g),
			temperature(t),size(s),mass_h(m),
			p_d(NULL),v_d(NULL),a_d(NULL),NBFconstants_d(NULL),NBPconstants_d(NULL),
			mass_d(NULL)
		{
			
		}
		
		~state()
		{
			if(p_d!=NULL) CUDA_API_Warnings(cudaFree(p_d));
			if(v_d!=NULL) CUDA_API_Warnings(cudaFree(v_d));
			if(a_d!=NULL) CUDA_API_Warnings(cudaFree(a_d));
			if(NBFconstants_d!=NULL) CUDA_API_Warnings(cudaFree(NBFconstants_d));
			if(NBPconstants_d!=NULL) CUDA_API_Warnings(cudaFree(NBPconstants_d));
			if(mass_d!=NULL) CUDA_API_Warnings(cudaFree(mass_d));
		}
		
		void resize(threeVector<T> s)
		{
			size=s;
		}
		
		void toDevice()
		{
			if(p_d==NULL) 
			CUDA_API_Errors( cudaMalloc((void **)&p_d, nParticles*sizeof(position<T>)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) p_d, (char*)p_h, 
				nParticles*sizeof(position<T>), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			if(v_d==NULL) 
			CUDA_API_Errors( cudaMalloc((void **)&v_d, nParticles*sizeof(threeVector<T>)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) v_d, (char*)v_h, 
				nParticles*sizeof(threeVector<T>), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			if(a_d==NULL) 
			CUDA_API_Errors( cudaMalloc((void **)&a_d, nParticles*sizeof(threeVector<T>)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) a_d, (char*)a_h, 
				nParticles*sizeof(threeVector<T>), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			if(NBFconstants_d==NULL) 
			CUDA_API_Errors( cudaMalloc((void **)&NBFconstants_d, 6*nTypes*nTypes*sizeof(T)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) NBFconstants_d, (char*)NBFconstants_h, 
				6*nTypes*nTypes*sizeof(T), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			if(NBPconstants_d==NULL) 
			CUDA_API_Errors( cudaMalloc((void **)&NBPconstants_d, 6*nTypes*nTypes*sizeof(T)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) NBPconstants_d, (char*)NBPconstants_h, 
				6*nTypes*nTypes*sizeof(T), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			if(mass_d==NULL) 
			CUDA_API_Errors( cudaMalloc((void **)&mass_d, nParticles*sizeof(T)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) mass_d, (char*)mass_h, 
				nParticles*sizeof(T), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		void toHost()
		{
			CUDA_API_Errors(cudaMemcpy((char *) p_h, (char*)p_d, 
				nParticles*sizeof(position<T>), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) v_h, (char*)v_d, 
				nParticles*sizeof(threeVector<T>), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) a_h, (char*)a_d, 
				nParticles*sizeof(threeVector<T>), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) NBFconstants_h, (char*)NBFconstants_d, 
				6*nTypes*nTypes*sizeof(T), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) NBPconstants_h, (char*)NBPconstants_d, 
				6*nTypes*nTypes*sizeof(T), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) mass_h, (char*)mass_d, 
				nParticles*sizeof(T), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		//pass to the host/device functions, this doesn't allocate memory
		struct copyState {
			using value_type=T;
			copyState(position<T> *pos, threeVector<T> *vel, threeVector<T> *acc, T *NBFc, T *NBPc,
				  int nP, int nT, T dT, T g, T t, threeVector<T> s, T *m):
				p(pos),v(vel),a(acc),NBFconstants(NBFc),NBPconstants(NBPc),
				nParticles(nP),nTypes(nT),deltaT(dT),gamma(g),
				temperature(t),size(s),mass(m)
			{}
			
			//Device pointers
			threeVector<T> *v,*a;
			position<T> *p;
			T *NBFconstants;
			T *NBPconstants;
			T *mass;
			
			int nParticles;
			int nTypes;
			threeVector<T> size;
			T deltaT;
			T gamma;
			T temperature;
		};
		
		copyState deviceState()
		{
			return copyState(p_d,v_d,a_d,NBFconstants_d,NBPconstants_d,nParticles,nTypes,deltaT,gamma,
					 temperature,size,mass_d);
		}
		
		copyState hostState()
		{
			return copyState(p_h,v_h,a_h,NBFconstants_h,NBPconstants_h,nParticles,nTypes,deltaT,gamma,
					 temperature,size,mass_h);
		}
		
		//Host pointers
		threeVector<T> *v_h,*a_h;
		position<T> *p_h;
		T *NBFconstants_h;
		T *NBPconstants_h;
		T *mass_h;
		uint *gridParticleIndex_h;//just allocate
		uint *gridParticleHash_h;//just allocate
		
		//Device pointers
		threeVector<T> *v_d,*a_d;
		position<T> *p_d;
		T *NBFconstants_d;
		T *NBPconstants_d;
		T *mass_d;
		uint *gridParticleIndex_d;//just allocate
		uint *gridParticleHash_d;//just allocate
		
		//Various simulation parameters
		int nParticles;
		int nTypes;
		threeVector<T> size;
		T deltaT;
		T gamma;
		T temperature;
	};
}

#endif
