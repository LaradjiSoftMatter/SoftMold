#include <vector>
#ifndef MPD_DATACOLLECTION
#define MPD_DATACOLLECTION

namespace mpd {
	
	template <typename T>
	struct dataCollection {
		using value_type=T;
		dataCollection(uint nP):nParticles(nP),
			potentialEnergy_h(nP),potentialEnergy_d(NULL),
			kineticEnergy_h(nP),kineticEnergy_d(NULL)
		{
			CUDA_API_Errors(cudaMalloc((void **)&potentialEnergy_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMalloc((void **)&kineticEnergy_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		}
		
		~dataCollection()
		{
			if(potentialEnergy_d!=NULL) CUDA_API_Warnings(cudaFree(potentialEnergy_d));
			if(potentialEnergy_d!=NULL) CUDA_API_Warnings(cudaFree(kineticEnergy_d));
		}
		
		void toDevice()
		{
			if(potentialEnergy_d==NULL) 
			CUDA_API_Errors(cudaMalloc((void **)&potentialEnergy_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) potentialEnergy_d, (char*)potentialEnergy_h.data(), 
				nParticles*sizeof(T), cudaMemcpyHostToDevice));
		CUDA_API_Errors(cudaDeviceSynchronize());
			if(kineticEnergy_d==NULL) 
			CUDA_API_Errors(cudaMalloc((void **)&kineticEnergy_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) kineticEnergy_d, (char*)kineticEnergy_h.data(), 
				nParticles*sizeof(T), cudaMemcpyHostToDevice));
		CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		void toHost()
		{
			CUDA_API_Errors(cudaMemcpy((char *) potentialEnergy_h.data(), (char*)potentialEnergy_d, 
				nParticles*sizeof(T), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaMemcpy((char *) kineticEnergy_h.data(), (char*)kineticEnergy_d, 
				nParticles*sizeof(T), cudaMemcpyDeviceToHost));
		}
		
		struct copyState;
		
		copyState deviceState()
		{
			return copyState(nParticles,potentialEnergy_d,kineticEnergy_d);
		}
		
		copyState hostState()
		{
			return copyState(nParticles,potentialEnergy_h.data(),kineticEnergy_h.data());
		}
		
		//pass to the host/device functions, this doesn't allocate memory
		struct copyState {
			using value_type=T;
			copyState(uint nP, T *pE, T *kE):nParticles(nP),potentialEnergy(pE),kineticEnergy(kE)
			{}
			
			//Device pointers
			T *potentialEnergy;
			T *kineticEnergy;
			
			uint nParticles;
		};
		
		//Host pointers
		std::vector<T> potentialEnergy_h;
		std::vector<T> kineticEnergy_h;
		
		//Device pointers
		T *potentialEnergy_d;
		T *kineticEnergy_d;
		
		//Various simulation parameters
		uint nParticles;
	};
}

#endif
