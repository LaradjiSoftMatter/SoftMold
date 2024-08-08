#include <cuda.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#ifndef MPD_REDUCEPOTENTIAL
#define MPD_REDUCEPOTENTIAL

namespace mpd {
	
	template <typename DATACOLLECTION>
	typename DATACOLLECTION::value_type reducePotential_host(DATACOLLECTION data)
	{
		typename DATACOLLECTION::value_type potentialEnergy=0.0;
		#pragma omp parallel for reduction(+:potentialEnergy)
		for(uint i=0;i<data.nParticles;i++)
			potentialEnergy+=data.potentialEnergy[i];
		return potentialEnergy;
	}
	
	template <typename DATACOLLECTION>
	typename DATACOLLECTION::value_type reducePotential_device(DATACOLLECTION data)
	{
		typename DATACOLLECTION::value_type potentialEnergy=0.0;
		potentialEnergy=thrust::reduce(thrust::device,data.potentialEnergy,
					data.potentialEnergy+data.nParticles,0.0);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		return potentialEnergy;
	}
}

#endif
