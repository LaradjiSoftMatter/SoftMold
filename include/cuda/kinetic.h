
#include <cuda.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include "dataTypes.h"
#include "other.h"
#ifndef MPD_KINETIC
#define MPD_KINETIC

namespace mpd {
	
	template <typename T>
	__host__ __device__
	T kinetic(const threeVector<T> &v)
	{
		return 0.5*(v.x*v.x+v.y*v.y+v.z*v.z);
	}
	
	template <typename STATE, typename DATACOLLECTION>
	__global__ 
	void kinetic_kernel(STATE input, DATACOLLECTION data)
	{
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=input.nParticles) return;
		data.kineticEnergy[i]+=kinetic(input.v[i]);
	}
	
	template <typename STATE,typename DATACOLLECTION>
	void kinetic_device(STATE input, DATACOLLECTION data)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=512;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		kinetic_kernel<<<numBlocks,numThreads>>>(input,data);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename STATE,typename DATACOLLECTION>
	void kinetic_host(STATE input, DATACOLLECTION data)
	{
		#pragma omp parallel for
		for(uint i=0;i<input.nParticles;i++)
			data.kineticEnergy[i]+=kinetic(input.v[i]);
	}
}

#endif
