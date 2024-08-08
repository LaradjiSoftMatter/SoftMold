#include <cuda.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#ifndef MPD_REDUCERANGE
#define MPD_REDUCERANGE

namespace mpd {
	
	template <typename T>
	T reduceRange_host(T *data, int n)
	{
		T total(0);
		#pragma omp parallel for reduction(+:total)
		for(uint i=0;i<n;i++)
			total+=data[i];
		return total;
	}
	
	template <typename T>
	T reduceRange_device(T *data, int n)
	{
		T total=thrust::reduce(thrust::device,data,data+n,T(0));
		CUDA_API_Errors(cudaDeviceSynchronize());
		return total;
	}
}

#endif
