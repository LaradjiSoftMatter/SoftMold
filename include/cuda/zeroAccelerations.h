#include "dataTypes.h"
#include "other.h"
#ifndef MPD_ZEROACCELERATIONS
#define MPD_ZEROACCELERATIONS

namespace mpd {
	
	template <typename STATE>
	__global__ 
	void zeroAccelerations_kernel(STATE input)
	{
		using T=typename STATE::value_type;
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=input.nParticles) return;
		input.a[i]=threeVector<T>(0.0,0.0,0.0);
	}
	
	template <typename STATE>
	void zeroAccelerations_device(STATE input)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=1024;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		zeroAccelerations_kernel<<<numBlocks,numThreads>>>(input);
		CUDA_Kernel_Errors();
		checkCudaErrors(cudaDeviceSynchronize());
	}
	
	template <typename STATE>
	void zeroAccelerations_host(STATE input)
	{
		using T=typename STATE::value_type;
		#pragma omp parallel for
		for(uint i=0;i<input.nParticles;i++)
			input.a[i]=threeVector<T>(0.0,0.0,0.0);
	}
}

#endif
