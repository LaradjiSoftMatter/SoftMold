
#include <cuda.h>
#include "potentials.h"
#include "dataTypes.h"
#include "other.h"
#include "errors.h"
#ifndef MPD_BENDINGINTERACTIONS
#define MPD_BENDINGINTERACTIONS

namespace mpd {

	template <typename BENDLIST, typename STATE>
	void bendForces_host(BENDLIST bendList, STATE input)
	{
		using T=typename STATE::value_type;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=i;
				bendList.elements[j].nextOffset!=bendList.endIndex();
				j=bendList.elements[j].nextOffset)
			{
				int *index=bendList.elements[j].index;
				threeVector<T> da=difference(input.p[index[0]],input.p[index[1]]);
				da=minImg(da,input.size);
				threeVector<T> db=difference(input.p[index[1]],input.p[index[2]]);
				db=minImg(db,input.size);
				
				aTotal+=bendF(da, db, bendList.elements[j].offset,
					      bendList.elements[j].constants);
			}
			input.a[i]+=aTotal;
		}
	}
	
	template <typename BENDLIST, typename STATE, typename DATACOLLECTION>
	void bendPotential_host(BENDLIST bendList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=i;
				bendList.elements[j].nextOffset!=bendList.endIndex();
				j=bendList.elements[j].nextOffset)
			{
				int *index=bendList.elements[j].index;
				threeVector<T> da=difference(input.p[index[0]],input.p[index[1]]);
				da=minImg(da,input.size);
				threeVector<T> db=difference(input.p[index[1]],input.p[index[2]]);
				db=minImg(db,input.size);
				
				dc.potentialEnergy[i]+=bendP(da, db, bendList.elements[j].offset, 
							     bendList.elements[j].constants)/3.0;
				
			}
		}
	}
	
	template <typename BENDLIST, typename STATE>
	__global__
	void bendingForces_kernel(BENDLIST bendList, STATE input)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		threeVector<T> aTotal(0.0,0.0,0.0);
		for(int j=i;
			bendList.elements[j].nextOffset!=bendList.endIndex();
			j=bendList.elements[j].nextOffset)
		{
			int *index=bendList.elements[j].index;
			threeVector<T> da=difference(input.p[index[0]],input.p[index[1]]);
			da=minImg(da,input.size);
			threeVector<T> db=difference(input.p[index[1]],input.p[index[2]]);
			db=minImg(db,input.size);
			
			aTotal+=bendF(da, db, bendList.elements[j].offset,
				      bendList.elements[j].constants);
		}
		input.a[i]+=aTotal;
	}
	
	template <typename BENDLIST, typename STATE>
	void bendForces_device(BENDLIST bendList, STATE input)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		bendingForces_kernel<<<numBlocks,numThreads>>>(bendList,input);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename BENDLIST, typename STATE, typename DATACOLLECTION>
	__global__
	void bendingPotential_kernel(BENDLIST bendList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		for(int j=i;
			bendList.elements[j].nextOffset!=bendList.endIndex();
			j=bendList.elements[j].nextOffset)
		{
			int *index=bendList.elements[j].index;
			threeVector<T> da=difference(input.p[index[0]],input.p[index[1]]);
			da=minImg(da,input.size);
			threeVector<T> db=difference(input.p[index[1]],input.p[index[2]]);
			db=minImg(db,input.size);
			
			dc.potentialEnergy[i]+=bendP(da, db, bendList.elements[j].offset,
						     bendList.elements[j].constants)/3.0;
		}
	}
	
	template <typename BENDLIST, typename STATE, typename DATACOLLECTION>
	void bendPotential_device(BENDLIST bendList, STATE input, DATACOLLECTION dc)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		bendingPotential_kernel<<<numBlocks,numThreads>>>(bendList,input,dc);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename BENDLIST, typename STATE, typename T>
	__global__
	void bendingDPotential_kernel(BENDLIST bendList, STATE input, T *dPotential, threeVector<T> scale)
	{
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		for(int j=i;
			bendList.elements[j].nextOffset!=bendList.endIndex();
			j=bendList.elements[j].nextOffset)
		{
			int *index=bendList.elements[j].index;
			threeVector<T> da=difference(input.p[index[0]],input.p[index[1]]);
			da=minImg(da,input.size);
			threeVector<T> db=difference(input.p[index[1]],input.p[index[2]]);
			db=minImg(db,input.size);
			
			T potentialA=bendP(da, db, bendList.elements[j].offset,
					   bendList.elements[j].constants)/3.0;
			da.x*=scale.x;
			da.y*=scale.y;
			da.z*=scale.z;
			db.x*=scale.x;
			db.y*=scale.y;
			db.z*=scale.z;
			T potentialB=bendP(da, db, bendList.elements[j].offset,
					   bendList.elements[j].constants)/3.0;
			dPotential[i]+=(potentialA-potentialB);
		}
	}
	
	template <typename BENDLIST, typename STATE, typename BAROSTAT, typename SCALE>
	void bendDPotential_device(BENDLIST bendList, STATE input, BAROSTAT bState, SCALE scale)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		bendingDPotential_kernel<<<numBlocks,numThreads>>>(bendList,input,bState.dPotential,scale);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
}

#endif

