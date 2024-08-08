#include "dataTypes.h"
#include "other.h"
#ifndef MPD_RESCALE
#define MPD_RESCALE

namespace mpd {
	
	template <typename STATE, typename T>
	__global__ 
	void rescale_kernel(STATE input, T scale)
	{
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=input.nParticles) return;
		input.p[i].x*=scale.x;
		input.p[i].y*=scale.y;
		input.p[i].z*=scale.z;
	}
	
	template <typename STATE, typename T>
	void rescale_device(STATE input, T scale)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=1024;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		rescale_kernel<<<numBlocks,numThreads>>>(input,scale);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename STATE, typename T>
	void rescale_host(STATE input, T scale)
	{
		#pragma omp parallel for
		for(uint i=0;i<input.nParticles;i++)
		{
			input.p[i].x*=scale.x;
			input.p[i].y*=scale.y;
			input.p[i].z*=scale.z;
		}
	}
}

#endif
