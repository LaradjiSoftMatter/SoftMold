#include <cuda.h>
#include "errors.h"
#ifndef MPD_ZEROPOTENTIAL
#define MPD_ZEROPOTENTIAL

namespace mpd {
	
	template <typename DATACOLLECTION>
	void zeroPotential_host(DATACOLLECTION data)
	{
		#pragma omp parallel for
		for(uint i=0;i<data.nParticles;i++)
			data.potentialEnergy[i]=0.0;
	}
	template <typename DATACOLLECTION>
	__global__
	void zeroPotential_kernel(DATACOLLECTION data)
	{
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=data.nParticles) return;
		data.potentialEnergy[i]=0.0;
	}
	
	template <typename DATACOLLECTION>
	void zeroPotential_device(DATACOLLECTION data)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(data.nParticles, blockSize, numBlocks, numThreads);
		
		zeroPotential_kernel<<<numBlocks,numThreads>>>(data);
		CUDA_Kernel_Errors();
		CUDA_API_Errors((cudaDeviceSynchronize());
	}
}

#endif
