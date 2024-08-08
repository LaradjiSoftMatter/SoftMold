#include <cuda.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include "potentials.h"
#include "dataTypes.h"
#include "other.h"
#include "errors.h"
#ifndef MPD_BALLINTERACTIONS
#define MPD_BALLINTERACTIONS

namespace mpd {
	
	template <typename BALLLIST, typename STATE>
	void ballForces_host(BALLLIST ballList, STATE input)
	{
		using T=typename STATE::value_type;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=i;
				ballList.elements[j].nextOffset!=ballList.endIndex();
				j=ballList.elements[j].nextOffset)
			{
				int *index=ballList.elements[j].index;
				threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
				d=minImg(d,input.size);
				
				threeVector<T> aTemp=harmonicHalfF(d, ballList.elements[j].offset,
					      ballList.elements[j].constants);
				ballList.aTemp[j]-=aTemp;
				aTotal+=aTemp;
			}
			input.a[i]+=aTotal;
		}
	}
	
	template <typename BALLLIST, typename STATE, typename DATACOLLECTION>
	void ballPotential_host(BALLLIST ballList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			for(int j=i;
				ballList.elements[j].nextOffset!=ballList.endIndex();
				j=ballList.elements[j].nextOffset)
			{
				int *index=ballList.elements[j].index;
				threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
				d=minImg(d,input.size);
				dc.potentialEnergy[i]+=harmonicHalfP(d, ballList.elements[j].offset,
								 ballList.elements[j].constants)/2.0f;
			}
		}
	}
	
	template <typename BALLLIST, typename STATE>
	__global__
	void ballForces_kernel(BALLLIST ballList, STATE input)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		//threeVector<T> aTotal(0.0,0.0,0.0);
		for(int j=i;
		    ballList.elements[j].nextOffset!=ballList.endIndex();
		    j=ballList.elements[j].nextOffset)
		{
			int *index=ballList.elements[j].index;
			threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
			d=minImg(d,input.size);
			
			threeVector<T> aTemp=harmonicHalfF(d, ballList.elements[j].offset,
			      ballList.elements[j].constants);
			ballList.aTemp[index[1]]+=aTemp;
			input.a[index[1]]-=aTemp;//this has to be correct
		}
	}
	
	template <typename T>
	__global__
	void singleBall_kernel(T *a, T input, int j)
	{
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		if (i > 0) return;
		a[j]+=input;
	}
	
	template <typename BALLLIST, typename STATE>
	void ballForces_device(BALLLIST ballList, STATE input)
	{
		if(ballList.size()>0)
		{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		using T=threeVector<typename STATE::value_type>;
		thrust::fill(thrust::device_ptr<T>(ballList.aTemp), 
			thrust::device_ptr<T>(ballList.aTemp+input.nParticles), T(0));
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		ballForces_kernel<<<numBlocks,numThreads>>>(ballList,input);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		T aTemp=thrust::reduce(thrust::device,thrust::device_ptr<T>(ballList.aTemp), 
			thrust::device_ptr<T>(ballList.aTemp+input.nParticles),T(0),
			tvPlus<typename STATE::value_type>());
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		singleBall_kernel<<<1,1>>>(input.a,aTemp,ballList.center);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		}
	}
	
	template <typename BALLLIST, typename STATE, typename DATACOLLECTION>
	__global__
	void ballPotential_kernel(BALLLIST ballList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		for(int j=i;
		    ballList.elements[j].nextOffset!=ballList.endIndex();
		    j=ballList.elements[j].nextOffset)
		{
			int *index=ballList.elements[j].index;
			threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
			d=minImg(d,input.size);
			dc.potentialEnergy[i]+=harmonicHalfP(d, ballList.elements[j].offset,
						 ballList.elements[j].constants);
		}
	}
	
	template <typename BALLLIST, typename STATE, typename DATACOLLECTION>
	void ballPotential_device(BALLLIST ballList, STATE input, DATACOLLECTION dc)
	{
		if(ballList.size()>0)
		{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		ballPotential_kernel<<<numBlocks,numThreads>>>(ballList,input,dc);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		}
	}
	
	template <typename BALLLIST, typename STATE, typename T>
	__global__
	void ballDPotential_kernel(BALLLIST ballList, STATE input, T *dPotential, threeVector<T> scale)
	{
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		for(int j=i;
		    ballList.elements[j].nextOffset!=ballList.endIndex();
		    j=ballList.elements[j].nextOffset)
		{
			int *index=ballList.elements[j].index;
			threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
			d=minImg(d,input.size);
			T potentialA=harmonicHalfP(d, ballList.elements[j].offset,
					       ballList.elements[j].constants);
			d.x*=scale.x;
			d.y*=scale.y;
			d.z*=scale.z;
			T potentialB=harmonicHalfP(d, ballList.elements[j].offset,
					       ballList.elements[j].constants);
			dPotential[i]+=(potentialA-potentialB);
		}
	}
	
	template <typename BALLLIST, typename STATE, typename BAROSTAT, typename SCALE>
	void ballDPotential_device(BALLLIST ballList, STATE input, BAROSTAT bState, SCALE scale)
	{
		if(ballList.size()>0)
		{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		ballDPotential_kernel<<<numBlocks,numThreads>>>(ballList,input,bState.dPotential,scale);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		}
	}
}

#endif

