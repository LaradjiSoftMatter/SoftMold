#include "dataTypes.h"
#include "other.h"
#ifndef MPD_VERLET
#define MPD_VERLET

namespace mpd {
	
	//Should be (combining half step velocity update):
	//	1. Initialize accelerations (if restart, update velocities
	//	2. Update velocities "verletVelocity"
	//	3. Update positions "verletPosition"
	//	4. Update accelerations
	//	5. Update velocities "verletVelocity"
	//	6. Goto 2
	
	template <typename T>
	__host__ __device__
	position<T> verletPosition(position<T> p, const threeVector<T> &v, const T &dT, 
				   const threeVector<T> &s)
	{
		p.x+=v.x*dT;
		p.y+=v.y*dT;
		p.z+=v.z*dT;
		
		p=bounds(p,s);
		return p;
	}
	
	template <typename T>
	__host__ __device__
	threeVector<T> verletVelocity(threeVector<T> v, const threeVector<T> &a, const T &dT, 
				      const threeVector<T> &s)
	{
		v.x+=a.x*dT*0.5;
		v.y+=a.y*dT*0.5;
		v.z+=a.z*dT*0.5;
		return v;
	}
	
	template <typename STATE>
	void verletFirst_host(STATE input)
	{
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			input.v[i]=verletVelocity(input.v[i],input.a[i],input.deltaT,input.size);
			input.p[i]=verletPosition(input.p[i],input.v[i],input.deltaT,input.size);
		}
	}
	
	template <typename STATE>
	void verletSecond_host(STATE input)
	{
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
			input.v[i]=verletVelocity(input.v[i],input.a[i],input.deltaT,input.size);
	}

	template <typename STATE>
	__global__ 
	void verletFirst_kernel(STATE input)
	{
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=input.nParticles) return;
		input.v[i]=verletVelocity(input.v[i],input.a[i],input.deltaT,input.size);
		input.p[i]=verletPosition(input.p[i],input.v[i],input.deltaT,input.size);
	}
	
	template <typename STATE>
	void verletFirst_device(STATE input)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=1024;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		verletFirst_kernel<<<numBlocks,numThreads>>>(input);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename STATE>
	__global__ 
	void verletSecond_kernel(STATE input)
	{
		uint i=blockIdx.x* blockDim.x + threadIdx.x;
		if (i>=input.nParticles) return;
		input.v[i]=verletVelocity(input.v[i],input.a[i],input.deltaT,input.size);
	}
	
	template <typename STATE>
	void verletSecond_device(STATE input)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=1024;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		verletSecond_kernel<<<numBlocks,numThreads>>>(input);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
}
#endif
