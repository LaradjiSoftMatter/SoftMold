
#include <cuda.h>
#include <curand_kernel.h>
#include <curand_mtgp32_host.h>
#include <curand_mtgp32dc_p_11213.h>
#include <random>
#include "other.h"
#ifndef MPD_RANDOM
#define MPD_RANDOM


namespace mpd {
	
	struct randomAdaptor {
		//using randType=std::mt19937;
		using randType=curandStateMtgp32;
		using randParams=mtgp32_kernel_params;
		randomAdaptor(unsigned long long seed)
		{
			_nBlocks=64;
			_nThreads=256;
			cudaMalloc((void **)&_randKernelParams,sizeof(randParams));
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			cudaMalloc((void **)&_randState,_nBlocks*sizeof(randType));
			CUDA_API_Errors(cudaDeviceSynchronize());
			cudaDeviceSynchronize();
			
			CURAND_Errors(curandMakeMTGP32Constants(mtgp32dc_params_fast_11213,_randKernelParams));
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			CURAND_Errors(curandMakeMTGP32KernelState(_randState,
					mtgp32dc_params_fast_11213,_randKernelParams,_nBlocks,seed));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		~randomAdaptor()
		{
			if(_randState!=NULL) CUDA_API_Warnings(cudaFree(_randState));
			if(_randKernelParams!=NULL) CUDA_API_Warnings(cudaFree(_randKernelParams));
		}
		
		struct copyState;
		
		copyState deviceState()
		{
			return copyState(_randState,_randKernelParams,_nBlocks,_nThreads);
		}
		
		//Force this to work with <random> later...
		//I'm way to tired to do this.
		//copyState host()
		//{
		//	return copyState(_randState,_randKernelParams,_nBlocks,_nThreads);
		//}
		
		struct copyState {
			copyState(randType *r,randParams *rk,int nB, int nT):
				  randState(r),randKernelParams(rk),nBlocks(nB),nThreads(nT)
			{}
			randType *randState;
			randParams *randKernelParams;
			int nBlocks,nThreads;
		};
		
		private:
			randType *_randState;
			randParams *_randKernelParams;
			int _nBlocks;
			int _nThreads;
	};
	
}

#endif
