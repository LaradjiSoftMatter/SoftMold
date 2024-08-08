#include "dataTypes.h"
#include "other.h"
#ifndef MPD_APPLYMASS
#define MPD_APPLYMASS

namespace mpd {
	
	template <typename STATE>
	__global__ 
	void applyMass_kernel(STATE input)
	{
		using T=typename STATE::value_type;
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=input.nParticles) return;
		input.a[i]/=input.mass[i];
	}
	
	template <typename STATE>
	void applyMass_device(STATE input)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=1024;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		applyMass_kernel<<<numBlocks,numThreads>>>(input);
		CUDA_Kernel_Errors();
	}
	
	template <typename STATE>
	void applyMass_host(STATE input)
	{
		using T=typename STATE::value_type;
		#pragma omp parallel for
		for(uint i=0;i<input.nParticles;i++)
			input.a[i]/=input.mass[i];
	}
}

#endif
