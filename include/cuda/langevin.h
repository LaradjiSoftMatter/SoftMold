#include <omp.h>
#include <cuda.h>
#include <random>

#include "dataTypes.h"
#include "randomAdaptor.h"
#include "other.h"
#ifndef MPD_LANGEVIN
#define MPD_LANGEVIN

namespace mpd {
	
	template <typename T>
	__host__ __device__
	threeVector<T> langevin(const threeVector<T> &v, const T &dT, const T &t, 
				const T &g, const threeVector<T> randVector)
	{
		threeVector<T> a;
		T sigma=sqrt((6.0*t*g)/dT);
		T psi=2.0*randVector.x-1.0;
		a.x=-g*v.x+sigma*psi;
		psi=2.0*randVector.y-1.0;
		a.y=-g*v.y+sigma*psi;
		psi=2.0*randVector.z-1.0;
		a.z=-g*v.z+sigma*psi;
		return a;
	}
	
	template <typename T>
	__host__ __device__
	T langevinDim(const threeVector<T> &v, const T &dT, const T &t, 
				const T &g, T randVector, int dim)
	{
		T a(0.0f);
		T sigma=sqrt((6.0*t*g)/dT);
		T psi=2.0*randVector-1.0;
		a=-g*v.s[dim]+sigma*psi;
		return a;
	}
	
	template <typename STATE, typename R>
	__global__ 
	void langevin_kernel(STATE input, R randNum)
	{
		using T=typename STATE::value_type;
		uint offset=gridDim.x*blockDim.x;
		//It appears that this requires the full random state to be updated.
		for(uint i=blockIdx.x* blockDim.x + threadIdx.x;i<input.nParticles||i<offset;i+=offset)
		{
		threeVector<T> randVector(0.0f,0.0f,0.0f);
		randVector.x=curand_uniform(&randNum.randState[blockIdx.x]);
		randVector.y=curand_uniform(&randNum.randState[blockIdx.x]);
		randVector.z=curand_uniform(&randNum.randState[blockIdx.x]);
		if(i<input.nParticles)
		input.a[i]+=langevin(input.v[i],input.deltaT,input.temperature,input.gamma,randVector);
		}
	}
	
	template <typename STATE, typename R>
	__global__ 
	void langevinDim_kernel(STATE input, R randNum, int dim)
	{
		using T=typename STATE::value_type;
		uint offset=gridDim.x*blockDim.x;
		for(uint i=blockIdx.x* blockDim.x + threadIdx.x;i<input.nParticles;i+=offset)
		{
		T randVector;
		randVector=curand_uniform(&randNum.randState[blockIdx.x]);
		input.a[i].s[dim]+=langevinDim(input.v[i],input.deltaT,input.temperature,input.gamma,randVector,dim);
		}
	}
	
	template <typename STATE, typename R>
	void langevin_device(STATE input, R randNum)
	{
		langevin_kernel<<<randNum.nBlocks,randNum.nThreads>>>(input,randNum);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename STATE, typename R>
	void langevin_host(STATE input, R randNum)
	{
		using T=typename STATE::value_type;
		#pragma omp parallel for
		for(uint i=0;i<input.nParticles;i++)
		{
			int thread=omp_get_thread_num();
			std::uniform_real_distribution<double> dist(0.0,1.0);
			threeVector<T> randVector;
			randVector.x=dist(randNum.randState[thread]);
			randVector.y=dist(randNum.randState[thread]);
			randVector.z=dist(randNum.randState[thread]);
			input.a[i]+=langevin(input.v[i],input.deltaT,input.temperature,
					     input.gamma,randVector);
		}
	}
}

#endif
