#include <vector>
#ifndef MPD_DATACOLLECTION
#define MPD_DATACOLLECTION

namespace mpd {
	
	template <typename T>
	struct dataCollection {
		using value_type=T;
		dataCollection(uint nP):nParticles(nP),
			potentialEnergy_h(nP),potentialEnergy_d(NULL),
			kineticEnergy_h(nP),kineticEnergy_d(NULL),
			beadPotential_h(nP),beadPotential_d(NULL)
		{
			CUDA_API_Errors(cudaMalloc((void **)&potentialEnergy_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMalloc((void **)&kineticEnergy_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMalloc((void **)&beadPotential_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		}
		
		~dataCollection()
		{
			if(potentialEnergy_d!=NULL) CUDA_API_Warnings(cudaFree(potentialEnergy_d));
			if(kineticEnergy_d!=NULL) CUDA_API_Warnings(cudaFree(kineticEnergy_d));
			if(beadPotential_d!=NULL) CUDA_API_Warnings(cudaFree(beadPotential_d));
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
			if(beadPotential_d==NULL) 
			CUDA_API_Errors(cudaMalloc((void **)&beadPotential_d, nParticles*sizeof(T)) );
		CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) beadPotential_d, (char*)beadPotential_h.data(), 
				nParticles*sizeof(T), cudaMemcpyHostToDevice));
		CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		void toHost()
		{
			CUDA_API_Errors(cudaMemcpy((char *) potentialEnergy_h.data(), (char*)potentialEnergy_d, 
				nParticles*sizeof(T), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaMemcpy((char *) kineticEnergy_h.data(), (char*)kineticEnergy_d, 
				nParticles*sizeof(T), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaMemcpy((char *) beadPotential_h.data(), (char*)beadPotential_d, 
				nParticles*sizeof(T), cudaMemcpyDeviceToHost));
		}
		
		struct copyState;
		
		copyState deviceState()
		{
			return copyState(nParticles,potentialEnergy_d,kineticEnergy_d,beadPotential_d);
		}
		
		copyState hostState()
		{
			return copyState(nParticles,potentialEnergy_h.data(),kineticEnergy_h.data(),beadPotential_h.data());
		}
		
		//pass to the host/device functions, this doesn't allocate memory
		struct copyState {
			using value_type=T;
			copyState(uint nP, T *pE, T *kE, T *bP):nParticles(nP),potentialEnergy(pE),kineticEnergy(kE),beadPotential(bP)
			{}
			
			//Device pointers
			T *potentialEnergy;
			T *kineticEnergy;
			T *beadPotential;
			
			uint nParticles;
		};
		
		//Host pointers
		std::vector<T> potentialEnergy_h;
		std::vector<T> kineticEnergy_h;
		std::vector<T> beadPotential_h;
		
		//Device pointers
		T *potentialEnergy_d;
		T *kineticEnergy_d;
		T *beadPotential_d;
		
		//Various simulation parameters
		uint nParticles;
	};
}

#endif
