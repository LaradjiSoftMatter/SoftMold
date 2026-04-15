#include <cuda.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include "cell.h"
#include "../potentials/laradjiSpangler.h"
#include "dataTypes.h"
#include "zeroAccelerations.h"
#include "zeroPotential.h"
#include "other.h"
#include "errors.h"
#ifndef MPD_BEADINTERACTIONS
#define MPD_BEADINTERACTIONS

//This is the fourth contstant, 
// this should probably be obtained from MD.h.
// Refactor later
//#define BEADRADIUS 4

// This should perform:
// 1. bead to other particles interactions
// 2. other particles to bead interactions
// 3. bead to bead interactions

namespace mpd {

	template <typename BEADLIST, typename STATE>
	void beadForces_host(BEADLIST beadList, STATE input)
	{
		using T=typename STATE::value_type;
		
		//T rminD=c[0];
		//T cutoffSqr=input.cutoff*input.cutoff;
		T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
		int *index=beadList.elements;
		int *nbMask=beadList.nbMask;
		//everyone to each bead
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=0;j<beadList.last;j++)
			{
				int k=index[j];
				if(j!=k)
				{
					threeVector<T> d=difference(input.p[i],input.p[k]);
					d=minImg(d,input.size);
					if(nbMask[i]==0)//normal to bead
					{
						threeVector<T> a=laradjiSpanglerF(d, cutoffSqr, 
							beadList.constants, input.p[i].type, 
							input.p[k].type, input.nTypes);
						aTotal+=a;
					}
					else//bead to bead
					{
						threeVector<T> a=laradjiSpanglerBBF(d, cutoffSqr, 
							beadList.constants, input.p[i].type, 
							input.p[k].type, input.nTypes);
						aTotal+=a;
					}
				}
			}
			input.a[i]+=aTotal;
		}
		
		//bead to everyone, this is too heavy here, 
		//but shouldn't be called since it is just 
		//a template for the kernel below
		#pragma omp parallel for
		for(int i=0;i<beadList.last;i++)
		{
			int k=index[i];
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=0;j<input.nParticles;j++)
			{
				if(k!=j)
				{
					threeVector<T> d=difference(input.p[k],input.p[j]);
					d=minImg(d,input.size);
					threeVector<T> a=laradjiSpanglerF(d, cutoffSqr, beadList.constants, input.p[k].type, 
						input.p[j].type, input.nTypes);
					aTotal+=a;
				}
			}
			input.a[k]+=aTotal;
		}
	}
	
	template <typename BEADLIST, typename STATE, typename DATACOLLECTION>
	void beadPotential_host(BEADLIST beadList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		
		//T rminD=c[0];
		//T cutoffSqr=input.cutoff*input.cutoff;
		T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
		int *index=beadList.elements;
		int *nbMask=beadList.nbMask;
		
		#pragma omp parallel for
		for(int i=0;i<input.nParticles;i++)
		{
			threeVector<T> aTotal(0.0,0.0,0.0);
			for(int j=0;j<beadList.last;j++)
			{
				int k=index[j];
				if(j!=k)
				{
					threeVector<T> d=difference(input.p[i],input.p[k]);
					d=minImg(d,input.size);
					if(nbMask[i]==0)//normal to bead
					{
						T pot=laradjiSpanglerP(d, cutoffSqr, 
							beadList.constants, input.p[i].type, 
							input.p[k].type, input.nTypes);
						dc.potentialEnergy[i]+=pot;
						dc.beadPotential[i]+=pot;
					}
					else//bead to bead
					{
						T pot=laradjiSpanglerBBP(d, cutoffSqr, 
							beadList.constants, input.p[i].type, 
							input.p[k].type, input.nTypes);
						dc.potentialEnergy[i]+=pot;
						dc.beadPotential[i]+=pot;
					}
				}
			}
		}
		
		//bead to everyone, this is too heavy here, 
		//but shouldn't be called since it is just 
		//a template for the kernel below
		#pragma omp parallel for
		for(int i=0;i<beadList.last;i++)
		{
			int k=index[i];
			for(int j=0;j<input.nParticles;j++)
			{
				if(k!=j)
				{
					threeVector<T> d=difference(input.p[k],input.p[j]);
					d=minImg(d,input.size);
					T pot=laradjiSpanglerP(d, cutoffSqr, beadList.constants, input.p[k].type, 
						input.p[j].type, input.nTypes);
					dc.potentialEnergy[i]+=pot;
					dc.beadPotential[i]+=pot;
				}
			}
		}
	}
	
	/**
	 * This is for every particle to beads and handles
	 * bead-bead forces, too.
	 *
	 **/
	template <typename BEADLIST, typename STATE>
	__global__
	void beadForcesA_kernel(BEADLIST beadList, STATE input)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
		int *index=beadList.elements;
		int *nbMask=beadList.nbMask;
		
		threeVector<T> aTotal(0.0,0.0,0.0);
		for(int j=0;j<beadList.last;j++)
		{
			int k=index[j];
			if(i!=k)
			{
				threeVector<T> d=difference(input.p[i],input.p[k]);
				d=minImg(d,input.size);
				if(nbMask[i]==0)//normal to bead
				{
					threeVector<T> a=laradjiSpanglerF(d, cutoffSqr, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
					aTotal+=a;
				}
				else//bead to bead
				{
					T cutoffSqr2=beadList.constants[4]+beadList.constants[0];//2*BEADRADIUS+rc
					cutoffSqr2*=cutoffSqr2;
					threeVector<T> a=laradjiSpanglerBBF(d, cutoffSqr2, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
					aTotal+=a;
				}
			}
		}
		input.a[i]+=aTotal;
	}
	

	// collide a particle against all other particles in a given cell
	template <typename T>
	__host__ __device__
	threeVector<T> collideCellBeadForce(threeVector<int> gridPos,
				   threeVector<int> gridSize,
				   threeVector<T> wrapOffset,
				   int *index,
				   position<T> pos,
				   position<T> *oldPos,
				   uint *cellBegin,
				   uint *cellEnd,
				   T *constants,
				   T cutoffSqr,
				   int nTypes)
	{
		uint gridHash = calcGridHash(gridPos,gridSize);
		
		// get start of bucket for this cell
		uint startIndex = cellBegin[gridHash];
		
		threeVector<T> force(0.0f,0.0f,0.0f);// = make_float3(0.0f);
		if (startIndex != 0xffffffff) // cell is not empty
		{
			// iterate over particles in this cell
			uint endIndex = cellEnd[gridHash];
			
			for (uint j=startIndex; j<endIndex; j++)
			{
				position<T> pos2 = oldPos[j];
				threeVector<T> d;
				d.x=pos.x-pos2.x-wrapOffset.x;
				d.y=pos.y-pos2.y-wrapOffset.y;
				d.z=pos.z-pos2.z-wrapOffset.z;
				// collide two spheres if not the same
				if(d.x*d.x+d.y*d.y+d.z*d.z!=0.0 && pos.type!=pos2.type)
					force+=laradjiSpanglerF(d, cutoffSqr, constants, pos.type, 
								pos2.type, nTypes);
			}
		}
		
		return force;
	}
	
	/**
	 * This is for beads to nearby particles.
	 *
	 **/
	template <typename BEADLIST, typename STATE, typename CELL>
	__global__
	void beadForcesB_kernel(BEADLIST beadList, STATE input, CELL cData)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= beadList.last) return;
		
		T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
		int *index=beadList.elements;
		int j=index[i];
		// read particle data from unsorted array
		position<T> pos = input.p[j];
		
		// get address in grid
		threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
		
		// get max size search radius
		threeVector<int> nDist;
		nDist.x=1;//int(cutoffSqr/cData.cellSize.x)+1;
		nDist.y=1;//int(cutoffSqr/cData.cellSize.y)+1;
		nDist.z=1;//int(cutoffSqr/cData.cellSize.z)+1;
		
		threeVector<T> aTotal(0.0,0.0,0.0);
		for (int z=-nDist.z; z<=nDist.z; z++)
		{
		for (int y=-nDist.y; y<=nDist.y; y++)
		{
		for (int x=-nDist.x; x<=nDist.x; x++)
		{
			threeVector<T> wrapOffset(0.0f,0.0f,0.0f);
			threeVector<int> neiPos = gridPos;
			neiPos.x+=x;
			neiPos.y+=y;
			neiPos.z+=z;
			wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
			
			aTotal+=collideCellBeadForce(neiPos,cData.gridSize,wrapOffset,index,
						     pos,cData.p,cData.cellBegin,cData.cellEnd,
						     beadList.constants,cutoffSqr,input.nTypes);
		}
		}
		}
		
		input.a[j]+=aTotal;
	}
	
	/**
	 * This is for every particle to beads and handles
	 * bead-bead forces, too.
	 *
	 **/
	template <typename BEADLIST, typename STATE, typename DATACOLLECTION>
	__global__
	void beadPotential_kernel(BEADLIST beadList, STATE input, DATACOLLECTION dc)
	{
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		int *index=beadList.elements;
		int *nbMask=beadList.nbMask;
		
		T potential=0;
		for(int j=0;j<beadList.last;j++)
		{
			int k=index[j];
			if(i!=k)
			{
				threeVector<T> d=difference(input.p[i],input.p[k]);
				d=minImg(d,input.size);
				if(nbMask[i]==0)//normal to bead
				{
					T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
					potential+=laradjiSpanglerP(d, cutoffSqr, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
				}
				else//bead to bead
				{
					T cutoffSqr2=beadList.constants[4]+beadList.constants[0];//2*BEADRADIUS+rc
					cutoffSqr2*=cutoffSqr2;
					potential+=laradjiSpanglerBBP(d, cutoffSqr2, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
				}
			}
		}
		dc.potentialEnergy[i]+=potential;
		dc.beadPotential[i]+=potential;
	}

	/* //These next two shouldn't be needed since we would double count the potential energy
	// collide a particle against all other particles in a given cell
	template <typename T>
	__host__ __device__
	threeVector<T> collideCellBeadPotential(threeVector<int> gridPos,
				   threeVector<int> gridSize,
				   threeVector<T> wrapOffset,
				   uint index,
				   position<T> pos,
				   position<T> *oldPos,
				   uint *cellBegin,
				   uint *cellEnd,
				   T *constants,
				   T cutoffSqr,
				   uint nTypes)
	{
		uint gridHash = calcGridHash(gridPos,gridSize);
		
		// get start of bucket for this cell
		uint startIndex = cellBegin[gridHash];
		
		T potential=0;
		if (startIndex != 0xffffffff) // cell is not empty
		{
			// iterate over particles in this cell
			uint endIndex = cellEnd[gridHash];
			
			for (uint j=startIndex; j<endIndex; j++)
			{
				position<T> pos2 = oldPos[j];
				threeVector<T> d;
				d.x=pos.x-pos2.x-wrapOffset.x;
				d.y=pos.y-pos2.y-wrapOffset.y;
				d.z=pos.z-pos2.z-wrapOffset.z;
				// collide two spheres if not the same
				if(d.x*d.x+d.y*d.y+d.z*d.z!=0.0)
					potential+=laradjiSpanglerP(d, cutoffSqr, constants, pos.type, 
								pos2.type, nTypes);
			}
		}
		
		return potential;
	}
	*/
	/**
	 * This is for beads to nearby particles.
	 *
	 **/
	/*
	template <typename BEADLIST, typename STATE, typename CELL, typename DATACOLLECTION>
	__global__
	void beadPotentialB_kernel(BEADLIST beadList, STATE input, CELL cData, DATACOLLECTION dc)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		using T=typename STATE::value_type;
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= beadList.last) return;
		
		T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
		int *index=beadList.elements;
		// read particle data from unsorted array
		position<T> pos = input.p[index[i]];
		
		// get address in grid
		threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
		
		// get max size search radius
		threeVector<int> nDist;
		nDist.x=int(cutoffSqr/cData.cellSize.x)+1;
		nDist.y=int(cutoffSqr/cData.cellSize.y)+1;
		nDist.z=int(cutoffSqr/cData.cellSize.z)+1;
		
		T potential=0;
		for (int z=-nDist.z; z<=nDist.z; z++)
		{
		for (int y=-nDist.y; y<=nDist.y; y++)
		{
		for (int x=-nDist.x; x<=nDist.x; x++)
		{
			threeVector<T> wrapOffset(0.0f,0.0f,0.0f);
			threeVector<int> neiPos = gridPos;
			neiPos.x+=x;
			neiPos.y+=y;
			neiPos.z+=z;
			wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
			
			potential+=collideCellBeadPotential(neiPos,cData.gridSize,wrapOffset,index,
						     pos,cData.p,cData.cellBegin,cData.cellEnd,
						     constants,cutoffSqr,nTypes);
		}
		}
		}
		
		dc.potentialEnergy[i]+=potential;
	}
	*/
	template <typename BEADLIST, typename STATE, typename CELL>
	void beadForces_device(BEADLIST beadList, STATE input, CELL cData, const int &timeStep)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		uint smemSize=sizeof(uint)*(numThreads+1);
		//if(cData.timeStep!=timeStep || timeStep<=0)
		{
			emptyCells_device(cData);
			CUDA_API_Errors(cudaDeviceSynchronize());
			calcHash_kernel<<<numBlocks,numThreads>>>(cData,input.p);
			CUDA_Kernel_Errors();
			CUDA_API_Errors(cudaDeviceSynchronize());
			sortParticles_device(cData);
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			reorderDataAndFindCellStart_kernel<<<numBlocks,numThreads,smemSize>>>(cData,input.p,input.v);
			CUDA_Kernel_Errors();
			CUDA_API_Errors(cudaDeviceSynchronize());
			cData.timeStep=timeStep;
		}
		
		beadForcesA_kernel<<<numBlocks,numThreads>>>(beadList,input);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		numBlocks=0;
		numThreads=0;
		blockSize=128;
		computeGridSize(beadList.last, blockSize, numBlocks, numThreads);
		beadForcesB_kernel<<<numBlocks,numThreads>>>(beadList,input,cData);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		//This shouldn't be needed since we are ignoring the reordering for beadForcesB_kernel
		//vRedistributionByParticle_kernel<<<numBlocks,numThreads>>>(cData.gridParticleIndex, 
		//							   cData.a, 
		//							   input.a, 
		//							   input.nParticles);
		//CUDA_Kernel_Errors();
		//CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	
	template <typename BEADLIST, typename STATE, typename DATACOLLECTION>
	void beadPotential_device(BEADLIST beadList, STATE input, DATACOLLECTION dc)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		beadPotential_kernel<<<numBlocks,numThreads>>>(beadList,input,dc);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename BEADLIST, typename STATE, typename T>
	__global__
	void beadDPotential_kernel(BEADLIST beadList, STATE input, T *dPotential, threeVector<T> scale)
	{
		uint i = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (i >= input.nParticles) return;
		
		int *index=beadList.elements;
		int *nbMask=beadList.nbMask;
		
		T potential=0;
		for(int j=0;j<beadList.last;j++)
		{
			int k=index[j];
			if(i!=k)
			{
				threeVector<T> d=difference(input.p[i],input.p[k]);
				d=minImg(d,input.size);
				if(nbMask[i]==0)//normal to bead
				{
					T cutoffSqr=beadList.constants[0]*beadList.constants[0];//BEADRADIUS+rc
					T potentialA=laradjiSpanglerP(d, cutoffSqr, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
					d.x*=scale.x;
					d.y*=scale.y;
					d.z*=scale.z;
					T potentialB=laradjiSpanglerP(d, cutoffSqr, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
					dPotential[i]+=(potentialA-potentialB);
				}
				else//bead to bead
				{
					T cutoffSqr2=beadList.constants[4]+beadList.constants[0];//2*BEADRADIUS+rc
					cutoffSqr2*=cutoffSqr2;
					T potentialA=laradjiSpanglerBBP(d, cutoffSqr2, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
					d.x*=scale.x;
					d.y*=scale.y;
					d.z*=scale.z;
					T potentialB=laradjiSpanglerBBP(d, cutoffSqr2, 
						beadList.constants, input.p[i].type, 
						input.p[k].type, input.nTypes);
					dPotential[i]+=(potentialA-potentialB);
				}
				
			}
		}
	}
	
	template <typename BEADLIST, typename STATE, typename BAROSTAT, typename SCALE>
	void beadDPotential_device(BEADLIST beadList, STATE input, BAROSTAT bState, SCALE scale)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(input.nParticles, blockSize, numBlocks, numThreads);
		
		beadDPotential_kernel<<<numBlocks,numThreads>>>(beadList,input,bState.dPotential,scale);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}

}
#endif
