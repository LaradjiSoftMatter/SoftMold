#include <cuda.h>
#include "potentials.h"
#include "dataTypes.h"
#include "other.h"
#include "errors.h"
#ifndef MPD_BONDINGINTERACTIONS
#define MPD_BONDINGINTERACTIONS

namespace mpd {

	template <typename BONDLIST, typename STATE>
	void bondForces_host(BONDLIST bondList, STATE input)
	{
		using T=typename STATE::value_type;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=i;
				bondList.elements[j].nextOffset!=bondList.endIndex();
				j=bondList.elements[j].nextOffset)
			{
				int *index=bondList.elements[j].index;
				threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
				d=minImg(d,input.size);
				
				aTotal+=harmonicF(d, bondList.elements[j].offset,
					      bondList.elements[j].constants);
			}
			input.a[i]+=aTotal;
		}
	}
	
	template <typename BONDLIST, typename STATE, typename DATACOLLECTION>
	void bondPotential_host(BONDLIST bondList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			for(int j=i;
				bondList.elements[j].nextOffset!=bondList.endIndex();
				j=bondList.elements[j].nextOffset)
			{
				int *index=bondList.elements[j].index;
				threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
				d=minImg(d,input.size);
				dc.potentialEnergy[i]+=harmonicP(d, bondList.elements[j].offset,
								 bondList.elements[j].constants)/2.0f;
			}
		}
	}
	
	template <typename BONDLIST, typename STATE>
	__global__
	void bondForces_kernel(BONDLIST bondList, STATE input)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		threeVector<T> aTotal(0.0,0.0,0.0);
		for(int j=i;
		    bondList.elements[j].nextOffset!=bondList.endIndex();
		    j=bondList.elements[j].nextOffset)
		{
			int *index=bondList.elements[j].index;
			threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
			d=minImg(d,input.size);
			
			aTotal+=harmonicF(d, bondList.elements[j].offset,
					  bondList.elements[j].constants);
		}
		input.a[i]+=aTotal;
	}
	
	template <typename BONDLIST, typename STATE>
	void bondForces_device(BONDLIST bondList, STATE input)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		bondForces_kernel<<<numBlocks,numThreads>>>(bondList,input);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename BONDLIST, typename STATE, typename DATACOLLECTION>
	__global__
	void bondPotential_kernel(BONDLIST bondList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		for(int j=i;
		    bondList.elements[j].nextOffset!=bondList.endIndex();
		    j=bondList.elements[j].nextOffset)
		{
			int *index=bondList.elements[j].index;
			threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
			d=minImg(d,input.size);
			dc.potentialEnergy[i]+=harmonicP(d, bondList.elements[j].offset,
							 bondList.elements[j].constants)/2.0f;
		}
	}
	
	template <typename BONDLIST, typename STATE, typename DATACOLLECTION>
	void bondPotential_device(BONDLIST bondList, STATE input, DATACOLLECTION dc)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		bondPotential_kernel<<<numBlocks,numThreads>>>(bondList,input,dc);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename BONDLIST, typename STATE, typename T>
	__global__
	void bondDPotential_kernel(BONDLIST bondList, STATE input, T *dPotential, threeVector<T> scale)
	{
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		for(int j=i;
		    bondList.elements[j].nextOffset!=bondList.endIndex();
		    j=bondList.elements[j].nextOffset)
		{
			int *index=bondList.elements[j].index;
			threeVector<T> d=difference(input.p[index[0]],input.p[index[1]]);
			d=minImg(d,input.size);
			T potentialA=harmonicP(d, bondList.elements[j].offset,
					       bondList.elements[j].constants)/2.0f;
			d.x*=scale.x;
			d.y*=scale.y;
			d.z*=scale.z;
			T potentialB=harmonicP(d, bondList.elements[j].offset,
					       bondList.elements[j].constants)/2.0f;
			dPotential[i]+=(potentialA-potentialB);
		}
	}
	
	template <typename BONDLIST, typename STATE, typename BAROSTAT, typename SCALE>
	void bondDPotential_device(BONDLIST bondList, STATE input, BAROSTAT bState, SCALE scale)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		bondDPotential_kernel<<<numBlocks,numThreads>>>(bondList,input,bState.dPotential,scale);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
}

#endif

