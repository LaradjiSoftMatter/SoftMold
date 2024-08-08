#include <cuda.h>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#ifndef MPD_ZERORANGE
#define MPD_ZERORANGE

namespace mpd {
	
	template <typename T>
	void zeroRange_host(T *data, int n)
	{
		#pragma omp parallel for
		for(uint i=0;i<n;i++)
			data[i]=T(0.0);
	}
	
	template <typename T>
	void zeroRange_device(T *data, int n)
	{
		thrust::fill(thrust::device_ptr<T>(data),
			     thrust::device_ptr<T>(data+n),
			     T(0));
	}
}

#endif
