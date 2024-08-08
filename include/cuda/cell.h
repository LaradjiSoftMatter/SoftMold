
#include <cuda.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>

#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>
#include <thrust/fill.h>
#include <thrust/device_ptr.h>
#include <cooperative_groups.h>

namespace cg=cooperative_groups;

#include "dataTypes.h"
#include "potentials.h"
#include "other.h"

#ifndef MPD_CELL
#define MPD_CELL

//This algorithm should be the following:
//1. cell<float> cData(...)
//2. emptyCells(cData.deviceCell())
//3. calcHash(cData.deviceCell(),p)
//4. sortParticles(cData.deviceCell())
//5. reorderDataAndFindCellStart(cData.deviceCell(),p,v)
//6. collide(cData.deviceCell(),constants)
//7. forceRedistribution(cData.deviceCell(),a)
//8. goto 2

namespace mpd {
	
	__host__ __device__
	uint nCells(threeVector<int> gS)
	{
		return gS.x*gS.y*gS.z;
	}
	
	//std::vector<int> testCount;
	
	template <typename T>
	struct cell {
		cell(uint nP, T rc, threeVector<T> s, T dL):
			nParticles(nP),gridParticleIndex_h(nP),gridParticleHash_h(nP),sortedPos_h(nP),
			sortedVel_h(nP),sortedAcc_h(nP),cutoff(rc),cellBegin_d(),cellEnd_d(NULL),maxCells(0),
			potentialEnergy_h(nP),potentialEnergy_d(NULL)
		{
			gridSize.x=0;
			gridSize.y=0;
			gridSize.z=0;
			resize(s,dL);
			CUDA_API_Errors( cudaMalloc((void **)&gridParticleIndex_d, nP*sizeof(uint)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors( cudaMalloc((void **)&gridParticleHash_d, nP*sizeof(uint)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors( cudaMalloc((void **)&sortedPos_d, nP*sizeof(position<T>)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors( cudaMalloc((void **)&sortedVel_d, nP*sizeof(threeVector<T>)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors( cudaMalloc((void **)&sortedAcc_d, nP*sizeof(threeVector<T>)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors( cudaMalloc((void **)&potentialEnergy_d, nP*sizeof(T)) );
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		~cell()
		{
			if(gridParticleIndex_d!=NULL) CUDA_API_Warnings(cudaFree(gridParticleIndex_d));
			if(gridParticleHash_d!=NULL) CUDA_API_Warnings(cudaFree(gridParticleHash_d));
			if(sortedPos_d!=NULL) CUDA_API_Warnings(cudaFree(sortedPos_d));
			if(sortedVel_d!=NULL) CUDA_API_Warnings(cudaFree(sortedVel_d));
			if(sortedAcc_d!=NULL) CUDA_API_Warnings(cudaFree(sortedAcc_d));
			if(cellBegin_d!=NULL) CUDA_API_Warnings(cudaFree(cellBegin_d));
			if(cellEnd_d!=NULL) CUDA_API_Warnings(cudaFree(cellEnd_d));
			if(potentialEnergy_d!=NULL) CUDA_API_Warnings(cudaFree(potentialEnergy_d));
		}
		
		void resize(threeVector<T> s, T dL)
		{
			size=s;
			gridSize.x=s.x/(cutoff+dL);
			gridSize.y=s.y/(cutoff+dL);
			gridSize.z=s.z/(cutoff+dL);
			cellSize.x=s.x/gridSize.x;
			cellSize.y=s.y/gridSize.y;
			cellSize.z=s.z/gridSize.z;
			uint currentNCells=nCells(gridSize);
			if(currentNCells>maxCells)
			{
				maxCells=currentNCells*2;
				cellBegin_h.resize(maxCells);
				cellEnd_h.resize(maxCells);
				if(cellBegin_d!=NULL) CUDA_API_Errors(cudaFree(cellBegin_d));
				CUDA_API_Errors(cudaDeviceSynchronize());
				if(cellEnd_d!=NULL) CUDA_API_Errors(cudaFree(cellEnd_d));
				CUDA_API_Errors(cudaDeviceSynchronize());
				CUDA_API_Errors( cudaMalloc((void **)&cellBegin_d, maxCells*sizeof(uint)) );
				CUDA_API_Errors(cudaDeviceSynchronize());
				CUDA_API_Errors( cudaMalloc((void **)&cellEnd_d, maxCells*sizeof(uint)) );
				CUDA_API_Errors(cudaDeviceSynchronize());
			}
		}
		
		//host vectors
		std::vector<uint> gridParticleIndex_h;
		std::vector<uint> gridParticleHash_h;
		std::vector<uint> cellBegin_h;
		std::vector<uint> cellEnd_h;
		std::vector<position<T>> sortedPos_h;
		std::vector<threeVector<T>> sortedVel_h;
		std::vector<threeVector<T>> sortedAcc_h;
		std::vector<T> potentialEnergy_h;
		
		//device pointers
		uint *gridParticleIndex_d;
		uint *gridParticleHash_d;
		uint *cellBegin_d;
		uint *cellEnd_d;
		position<T> *sortedPos_d;
		threeVector<T> *sortedVel_d;
		threeVector<T> *sortedAcc_d;
		T *potentialEnergy_d;
		
		//other parameters
		threeVector<int> gridSize;
		threeVector<T> cellSize;
		T cutoff;
		threeVector<T> size;
		uint nParticles;
		uint maxCells;
		
		struct copyCell;
		
		copyCell deviceCell()
		{
			return copyCell(gridParticleIndex_d,gridParticleHash_d,sortedPos_d,sortedVel_d,
					sortedAcc_d,potentialEnergy_d,cellBegin_d,cellEnd_d,gridSize,
					cellSize,size,cutoff,nParticles);
		}
		
		copyCell hostCell()
		{
			return copyCell(gridParticleIndex_h.data(),gridParticleHash_h.data(),
					sortedPos_h.data(),sortedVel_h.data(),sortedAcc_h.data(),
					potentialEnergy_h.data(),cellBegin_h.data(),cellEnd_h.data(),
					gridSize,cellSize,size,cutoff,nParticles);
		}
		
		struct copyCell {
			copyCell(uint *gpi, uint *gph, position<T> *sP, threeVector<T> *sV, 
				 threeVector<T> *sA, T* pE, uint *cB, uint *cE,threeVector<int> gS,
				 threeVector<T> cS, threeVector<T> s, T rc, uint nP):
				gridParticleIndex(gpi),gridParticleHash(gph),p(sP),v(sV),a(sA),
				potentialEnergy(pE),cellBegin(cB),cellEnd(cE),gridSize(gS),
				cellSize(cS),cutoff(rc),size(s),nParticles(nP)
			{}
			
			uint *gridParticleIndex;
			uint *gridParticleHash;
			uint *cellBegin;
			uint *cellEnd;
			position<T> *p;
			threeVector<T> *v;
			threeVector<T> *a;
			T *potentialEnergy;
			threeVector<int> gridSize;
			threeVector<T> cellSize;
			threeVector<T> size;
			T cutoff;
			uint nParticles;
		};
	};
	
	template <typename T>
	__host__ __device__
	void wrapWithOffset(threeVector<int> &neiPos, threeVector<T> &wrapOffset, 
			    threeVector<int> gridSize, threeVector<T> size)
	{
		//This does two things:
		//1. Wraps grid to stay within boundary
		//2. Includes the wraping offset for later calculations
		if(neiPos.x>=gridSize.x) {neiPos.x-=gridSize.x;wrapOffset.x=size.x;}
		if(neiPos.y>=gridSize.y) {neiPos.y-=gridSize.y;wrapOffset.y=size.y;}
		if(neiPos.z>=gridSize.z) {neiPos.z-=gridSize.z;wrapOffset.z=size.z;}
		if(neiPos.x<0) {neiPos.x+=gridSize.x;wrapOffset.x=-size.x;}
		if(neiPos.y<0) {neiPos.y+=gridSize.y;wrapOffset.y=-size.y;}
		if(neiPos.z<0) {neiPos.z+=gridSize.z;wrapOffset.z=-size.z;}
	}
	// calculate position in uniform grid
	template <typename T>
	__host__ __device__ 
	threeVector<int> calcGridPos(position<T> p, threeVector<T> cellSize)
	{
		threeVector<int> gridPos;
		gridPos.x = floor(p.x/cellSize.x);
		gridPos.y = floor(p.y/cellSize.y);
		gridPos.z = floor(p.z/cellSize.z);
		return gridPos;
	}
	
	// calculate address in grid from position (clamping to edges)
	__host__ __device__ 
	uint calcGridHash(threeVector<int> gPos, threeVector<int> gSize)
	{
		//	hashIndices[i].key=x+y*nCells.x+z*nCells.x*nCells.y;
		return gPos.x+gPos.y*gSize.x+gPos.z*gSize.x*gSize.y;
		//return gPos.z*gSize.y*gSize.x+gPos.y*gSize.x+gPos.x;
	}
	
	// calculate grid hash value for each particle
	template <typename T, typename CELL>
	__global__
	void calcHash_kernel(CELL cData, position<T> *p)
	{
		//uint index=blockIdx.x*blockDim.x+threadIdx.x;
		
		//if (index>=cData.nParticles) return;
		uint offset=gridDim.x*blockDim.x;
		for(uint index=blockIdx.x* blockDim.x + threadIdx.x;index<cData.nParticles;index+=offset)
		{
		
		volatile position<T> pos = p[index];
		
		// get address in grid
		threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
		uint hash = calcGridHash(gridPos,cData.gridSize);
		
		// store grid hash and particle index
		cData.gridParticleHash[index] = hash;
		cData.gridParticleIndex[index] = index;
		}
	}
	
	// calculate grid hash value for each particle
	template <typename T, typename CELL>
	void calcHash_host(CELL cData, position<T> *p)
	{
		#pragma omp parallel for
		for(uint index=0;index<cData.nParticles;index++)
		{
			position<T> pos = p[index];
			
			// get address in grid
			threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
			uint hash = calcGridHash(gridPos,cData.gridSize);
			// store grid hash and particle index
			cData.gridParticleHash[index] = hash;
			cData.gridParticleIndex[index] = index;
		}
	}
	
	template <typename CELL>
	void sortParticles_device(CELL cData)
	{
		thrust::sort_by_key(thrust::device_ptr<uint>(cData.gridParticleHash),
				    thrust::device_ptr<uint>(cData.gridParticleHash+cData.nParticles),
				    thrust::device_ptr<uint>(cData.gridParticleIndex));
	}
	
	template <typename CELL>
	void sortParticles_host(CELL cData)
	{
		using pT=std::pair<uint,uint>;
		std::vector<pT> makeSort(cData.nParticles);
		for(int i=0;i<cData.nParticles;i++)
			makeSort[i]=std::make_pair(cData.gridParticleHash[i],cData.gridParticleIndex[i]);
		
		std::sort(makeSort.begin(),makeSort.end(),
				    [](pT a, pT b) {return a.first<b.first;}
		);
		
		for(int i=0;i<cData.nParticles;i++)
		{
			cData.gridParticleHash[i]=makeSort[i].first;
			cData.gridParticleIndex[i]=makeSort[i].second;
		}
	}
	
	template <typename CELL>
	__global__
	void emptyCells_kernel(CELL cData)
	{
		//uint index = blockIdx.x*blockDim.x+threadIdx.x;
		int totalCells=nCells(cData.gridSize);
		//if(index>=totalCells) return;
		uint offset=gridDim.x*blockDim.x;
		for(uint index=blockIdx.x* blockDim.x + threadIdx.x;index<totalCells;index+=offset)
		{
			cData.cellBegin[index]=0xFFFFFFFF;
			cData.cellEnd[index]=0xFFFFFFFF;
		}
	}
	
	template <typename CELL>
	void emptyCells_device(CELL cData)
	{
		int totalCells=nCells(cData.gridSize);
		thrust::fill(thrust::device_ptr<uint>(cData.cellBegin),
			     thrust::device_ptr<uint>(cData.cellBegin+totalCells),
			     0xFFFFFFFF);
		thrust::fill(thrust::device_ptr<uint>(cData.cellEnd),
			     thrust::device_ptr<uint>(cData.cellEnd+totalCells),
			     0xFFFFFFFF);
	}
	
	
	template <typename CELL>
	void emptyCells_host(CELL cData)
	{
		int totalCells=nCells(cData.gridSize);
		for(uint i=totalCells-1;i<totalCells;--i)
			cData.cellBegin[i]=0xFFFFFFFF;
	}
	
	// rearrange particle data into sorted order, and find the start of each cell
	// in the sorted hash array
	template <typename T, typename CELL>
	__global__
	void reorderDataAndFindCellStart_kernel(CELL cData, position<T> *p, threeVector<T> *v)
	{
		//threeVector<uint> blockIdx,threadIdx;
		// Handle to thread block group
		cg::thread_block cta = cg::this_thread_block();
		extern __shared__ uint sharedHash[];	// blockSize + 1 elements
		uint index = blockIdx.x*blockDim.x+threadIdx.x;
		
		uint hash;
		
		// handle case when no. of particles not multiple of block size
		if (index < cData.nParticles)
		{
			hash = cData.gridParticleHash[index];
			
			// Load hash data into shared memory so that we can look
			// at neighboring particle's hash value without loading
			// two hash values per thread
			sharedHash[threadIdx.x+1] = hash;
			
			if (index > 0 && threadIdx.x == 0)
			{
				// first thread in block must load neighbor particle hash
				sharedHash[0] = cData.gridParticleHash[index-1];
			}
		}
		//__syncthreads();
		cg::sync(cta);
		
		if (index < cData.nParticles)
		{
			if (index == 0 || hash != sharedHash[threadIdx.x])
			{
				cData.cellBegin[hash] = index;
				if (index > 0)
					cData.cellEnd[sharedHash[threadIdx.x]] = index;
			}
			
			if (index==cData.nParticles-1)
				cData.cellEnd[hash]=index+1;
			
			// Now use the sorted index to reorder the pos and vel data
			uint sortedIndex = cData.gridParticleIndex[index];
			position<T> pos = p[sortedIndex];
			//threeVector<T> vel = v[sortedIndex];
			
			cData.p[index] = pos;
			//sortedVel[index] = vel;
		}
	}
	
	// rearrange particle data into sorted order, and find the start of each cell
	// in the sorted hash array
	template <typename T, typename CELL>
	void reorderDataAndFindCellStart_host(CELL cData, position<T> *p, threeVector<T> *v)
	{
		uint lastHash=0xFFFFFFFF;
		for(uint index=0;index<cData.nParticles;index++)
		{
			uint hash = cData.gridParticleHash[index];
			if(index==0 || hash !=lastHash)
			{
				cData.cellBegin[hash]=index;
				if(index>0)
					cData.cellEnd[lastHash]=index;
			}
			if(index==cData.nParticles-1)
				cData.cellEnd[hash]=index+1;
			cData.p[index] = p[cData.gridParticleIndex[index]];
			//cData.sortedVel[index] = v[cData.gridParticleIndex[index]];
			lastHash=hash;
		}
	}
	
	//redistribute/accumulate values from sorted grid index
	template <typename GPI, typename V>
	__global__
	void vRedistributionByParticle_kernel(GPI *gridParticleIndex, V *vIn, V *vOut, int nParticles)
	{
		uint index=blockIdx.x*blockDim.x+threadIdx.x;
		if(index>=nParticles) return;
		//uint offset=gridDim.x*blockDim.x;
		//for(uint index=blockIdx.x* blockDim.x + threadIdx.x;index<cData.nParticles;index+=offset)
		//{
		uint sortedIndex=gridParticleIndex[index];
		vOut[sortedIndex]+=vIn[index];
		//}
	}
	
	// calculate grid hash value for each particle
	template <typename GPI, typename V>
	void vRedistributionByParticle_host(GPI *gridParticleIndex, V *vIn, V *vOut, int nParticles)
	{
		#pragma omp parallel for
		for(uint index=0;index<nParticles;index++)
		{	
			uint sortedIndex=gridParticleIndex[index];
			vOut[sortedIndex]+=vIn[index];
		}
	}
	
	// collide a particle against all other particles in a given cell
	template <typename T>
	__host__ __device__
	threeVector<T> collideCellForce(threeVector<int> gridPos,
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
		
		threeVector<T> force(0.0f,0.0f,0.0f);// = make_float3(0.0f);
		if (startIndex != 0xffffffff) // cell is not empty
		{
			// iterate over particles in this cell
			uint endIndex = cellEnd[gridHash];
			
			for (uint j=startIndex; j<endIndex; j++)
			{
				if (j != index) // check not colliding with self
				{
					position<T> pos2 = oldPos[j];
					threeVector<T> d;
					d.x=pos.x-pos2.x-wrapOffset.x;
					d.y=pos.y-pos2.y-wrapOffset.y;
					d.z=pos.z-pos2.z-wrapOffset.z;
					// collide two spheres
					force+= nonBondedF(d, cutoffSqr, constants,
							   pos.type, pos2.type, nTypes);
				}
			}
		}
		
		return force;
	}
	
	// collide a particle against all other particles in a given cell
	template <typename T>
	__host__ __device__
	T collideCellPotential(threeVector<int> gridPos,
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
				if (j != index) // check not colliding with self
				{
					position<T> pos2 = oldPos[j];
					threeVector<T> d;
					d.x=pos.x-pos2.x-wrapOffset.x;
					d.y=pos.y-pos2.y-wrapOffset.y;
					d.z=pos.z-pos2.z-wrapOffset.z;
					// collide two spheres
					potential+= nonBondedP(d, cutoffSqr, constants,
							   pos.type, pos2.type, nTypes)/2.0;
				}
			}
		}
		
		return potential;
	}
	
	// collide a particle against all other particles in a given cell
	template <typename T>
	__host__ __device__
	T collideCellDPotential(threeVector<int> gridPos,
				   threeVector<int> gridSize,
				   threeVector<T> wrapOffset,
				   uint index,
				   position<T> pos,
				   position<T> *oldPos,
				   uint *cellBegin,
				   uint *cellEnd,
				   T *constants,
				   T cutoffSqr,
				   uint nTypes,
				   threeVector<T> scale)
	{
		uint gridHash = calcGridHash(gridPos,gridSize);
		
		// get start of bucket for this cell
		uint startIndex = cellBegin[gridHash];
		
		T dPotential=0;
		if (startIndex != 0xffffffff) // cell is not empty
		{
			// iterate over particles in this cell
			uint endIndex = cellEnd[gridHash];
			
			for (uint j=startIndex; j<endIndex; j++)
			{
				if (j != index) // check not colliding with self
				{
					position<T> pos2 = oldPos[j];
					threeVector<T> d;
					d.x=pos.x-pos2.x-wrapOffset.x;
					d.y=pos.y-pos2.y-wrapOffset.y;
					d.z=pos.z-pos2.z-wrapOffset.z;
					// collide two spheres
					T potentialA= nonBondedP(d, cutoffSqr, constants,
							   pos.type, pos2.type, nTypes)/2.0;
					d.x*=scale.x;
					d.y*=scale.y;
					d.z*=scale.z;
					T potentialB= nonBondedP(d, cutoffSqr, constants,
							   pos.type, pos2.type, nTypes)/2.0;
					dPotential+=(potentialA-potentialB);
				}
			}
		}
		
		return dPotential;
	}
	
	template <typename T, typename CELL>
	__global__
	void collideForce_kernel(CELL cData, T *constants, uint nTypes)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		uint index =blockIdx.x*blockDim.x + threadIdx.x;
		
		if (index >= cData.nParticles) return;
		
		//uint offset=gridDim.x*blockDim.x;
		//for(uint index=blockIdx.x* blockDim.x + threadIdx.x;index<cData.nParticles;index+=offset)
		//{
		// read particle data from sorted arrays
		position<T> pos = cData.p[index];
		
		// get address in grid
		threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
		
		// examine neighbouring cells
		threeVector<T> force(0.0f,0.0f,0.0f);// = make_float3(0.0f);
		for (int z=-1; z<=1; z++)
		{
		for (int y=-1; y<=1; y++)
		{
		for (int x=-1; x<=1; x++)
		{
			threeVector<T> wrapOffset(0.0f,0.0f,0.0f);
			threeVector<int> neiPos = gridPos;
			neiPos.x+=x;
			neiPos.y+=y;
			neiPos.z+=z;
			wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
			force+=collideCellForce(neiPos,cData.gridSize,wrapOffset,index,
						     pos,cData.p,cData.cellBegin,cData.cellEnd,
						     constants,cutoffSqr,nTypes);
		}
		}
		}
		cData.a[index]=force;
		//}
	}
	
	template <typename T, typename CELL>
	void collideForce_host(CELL cData, T *constants, uint nTypes)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		#pragma omp parallel for
		for(uint index=0;index<cData.nParticles;index++)
		{
			// read particle data from sorted arrays
			position<T> pos = cData.p[index];
			
			// get address in grid
			threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
			
			// examine neighbouring cells
			threeVector<T> force(0.0f,0.0f,0.0f);// = make_float3(0.0f);
			
			for (int z=-1; z<=1; z++)
			{
			for (int y=-1; y<=1; y++)
			{
			for (int x=-1; x<=1; x++)
			{
				threeVector<T> wrapOffset(0.0,0.0,0.0);
				threeVector<int> neiPos = gridPos;
				neiPos.x+=x;
				neiPos.y+=y;
				neiPos.z+=z;
				wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
				force+=collideCellForce(neiPos,cData.gridSize,wrapOffset,index,
							     pos,cData.p,cData.cellBegin,cData.cellEnd,
							     constants,cutoffSqr,nTypes);
			}
			}
			}
			cData.a[index]=force;
		}
	}
	
	template <typename T, typename CELL>
	__global__
	void collidePotential_kernel(CELL cData, T *constants, uint nTypes)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		uint index = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (index >= cData.nParticles) return;
		
		//uint offset=gridDim.x*blockDim.x;
		//for(uint index=blockIdx.x* blockDim.x + threadIdx.x;index<cData.nParticles;index+=offset)
		//{
		// read particle data from sorted arrays
		position<T> pos = cData.p[index];
		
		// get address in grid
		threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
		
		// examine neighbouring cells
		T potential=0;
		
		for (int z=-1; z<=1; z++)
		{
		for (int y=-1; y<=1; y++)
		{
		for (int x=-1; x<=1; x++)
		{
			threeVector<T> wrapOffset(0.0,0.0,0.0);
			threeVector<int> neiPos = gridPos;
			neiPos.x+=x;
			neiPos.y+=y;
			neiPos.z+=z;
			wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
			potential+=collideCellPotential(neiPos,cData.gridSize,wrapOffset,index,
						     pos,cData.p,cData.cellBegin,cData.cellEnd,
						     constants,cutoffSqr,nTypes);
		}
		}
		}
		cData.potentialEnergy[index]=potential;
		//}
	}
	
	
	template <typename T, typename CELL>
	void collidePotential_host(CELL cData, T *constants, uint nTypes)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		#pragma omp parallel for
		for(uint index=0;index<cData.nParticles;index++)
		{
			// read particle data from sorted arrays
			position<T> pos = cData.p[index];
			
			// get address in grid
			threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
			
			// examine neighbouring cells
			T potential=0;
			
			for (int z=-1; z<=1; z++)
			{
			for (int y=-1; y<=1; y++)
			{
			for (int x=-1; x<=1; x++)
			{
				threeVector<T> wrapOffset(0.0,0.0,0.0);
				threeVector<int> neiPos = gridPos;
				neiPos.x+=x;
				neiPos.y+=y;
				neiPos.z+=z;
				wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
				potential+=collideCellPotential(neiPos,cData.gridSize,wrapOffset,index,
							     pos,cData.p,cData.cellBegin,cData.cellEnd,
							     constants,cutoffSqr,nTypes);
			}
			}
			}
			cData.potentialEnergy[index]=potential;
		}
	}
	
	template <typename T, typename CELL>
	void collideDPotential_host(CELL cData, T *constants, uint nTypes, threeVector<T> scale)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		#pragma omp parallel for
		for(uint index=0;index<cData.nParticles;index++)
		{
			// read particle data from sorted arrays
			position<T> pos = cData.p[index];
			
			// get address in grid
			threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
			
			// examine neighbouring cells
			T potential=0;
			
			for (int z=-1; z<=1; z++)
			{
			for (int y=-1; y<=1; y++)
			{
			for (int x=-1; x<=1; x++)
			{
				threeVector<T> wrapOffset(0.0,0.0,0.0);
				threeVector<int> neiPos = gridPos;
				neiPos.x+=x;
				neiPos.y+=y;
				neiPos.z+=z;
				wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
				potential+=collideCellDPotential(neiPos,cData.gridSize,wrapOffset,index,
							     pos,cData.p,cData.cellBegin,cData.cellEnd,
							     constants,cutoffSqr,nTypes,scale);
			}
			}
			}
			cData.potentialEnergy[index]=potential;
		}
	}
	
	template <typename T, typename CELL>
	__global__
	void collideDPotential_kernel(CELL cData, T *constants, uint nTypes, threeVector<T> scale)
	{
		T cutoffSqr=cData.cutoff*cData.cutoff;
		uint index = blockIdx.x*blockDim.x + threadIdx.x;
		
		if (index >= cData.nParticles) return;
		
		//uint offset=gridDim.x*blockDim.x;
		//for(uint index=blockIdx.x* blockDim.x + threadIdx.x;index<cData.nParticles;index+=offset)
		//{
		// read particle data from sorted arrays
		position<T> pos = cData.p[index];
		
		// get address in grid
		threeVector<int> gridPos = calcGridPos(pos,cData.cellSize);
		
		// examine neighbouring cells
		T potential=0;
		
		for (int z=-1; z<=1; z++)
		{
		for (int y=-1; y<=1; y++)
		{
		for (int x=-1; x<=1; x++)
		{
			threeVector<T> wrapOffset(0.0,0.0,0.0);
			threeVector<int> neiPos = gridPos;
			neiPos.x+=x;
			neiPos.y+=y;
			neiPos.z+=z;
			wrapWithOffset(neiPos, wrapOffset, cData.gridSize, cData.size);
			potential+=collideCellDPotential(neiPos,cData.gridSize,wrapOffset,index,
						     pos,cData.p,cData.cellBegin,cData.cellEnd,
						     constants,cutoffSqr,nTypes,scale);
		}
		}
		}
		cData.potentialEnergy[index]=potential;
		//}
	}
	
	
	template <typename CELL, typename STATE>
	void cellComputeForce_host(CELL cData, STATE state)
	{
		emptyCells_host(cData);
		calcHash_host(cData,state.p);
		sortParticles_host(cData);
		reorderDataAndFindCellStart_host(cData,state.p,state.v);
		collideForce_host(cData,state.NBFconstants,state.nTypes);
		vRedistributionByParticle_host(cData.gridParticleIndex,cData.a,state.a,state.nParticles);
	}
	
	template <typename CELL, typename STATE, typename DATACOLLECTION>
	void cellComputePotential_host(CELL cData, STATE state, DATACOLLECTION data)
	{
		emptyCells_host(cData);
		calcHash_host(cData,state.p);
		sortParticles_host(cData);
		reorderDataAndFindCellStart_host(cData,state.p,state.v);
		collidePotential_host(cData,state.NBPconstants,state.nTypes);
		vRedistributionByParticle_host(cData.gridParticleIndex,cData.potentialEnergy,
					       data.potentialEnergy,state.nParticles);
	}
	
	template <typename CELL, typename STATE>
	void cellComputeForce_device(CELL cData, STATE state)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(state.nParticles, blockSize, numBlocks, numThreads);
		
		//emptyCells_kernel<<<4096,numThreads>>>(cData);
		emptyCells_device(cData);
		CUDA_API_Errors(cudaDeviceSynchronize());
		calcHash_kernel<<<numBlocks,numThreads>>>(cData,state.p);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		sortParticles_device(cData);
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		uint smemSize=sizeof(uint)*(numThreads+1);
		reorderDataAndFindCellStart_kernel<<<numBlocks,numThreads,smemSize>>>(cData,state.p,state.v);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		collideForce_kernel<<<numBlocks,numThreads>>>(cData,state.NBFconstants,state.nTypes);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		vRedistributionByParticle_kernel<<<numBlocks,numThreads>>>(cData.gridParticleIndex, 
									   cData.a, 
									   state.a, 
									   state.nParticles);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename CELL, typename STATE, typename DATACOLLECTION>
	void cellComputePotential_device(CELL cData, STATE state, DATACOLLECTION data)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(state.nParticles, blockSize, numBlocks, numThreads);
		
		//emptyCells_kernel<<<4096,numThreads>>>(cData);
		emptyCells_device(cData);
		CUDA_API_Errors(cudaDeviceSynchronize());
		calcHash_kernel<<<numBlocks,numThreads>>>(cData,state.p);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		sortParticles_device(cData);
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		uint smemSize=sizeof(uint)*(numThreads+1);
		reorderDataAndFindCellStart_kernel<<<numBlocks,numThreads,smemSize>>>(cData,state.p,state.v);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		collidePotential_kernel<<<numBlocks,numThreads>>>(cData,state.NBPconstants,state.nTypes);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		vRedistributionByParticle_kernel<<<numBlocks,numThreads>>>(cData.gridParticleIndex, 
									   cData.potentialEnergy, 
									   data.potentialEnergy, 
									   state.nParticles);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
	
	template <typename CELL, typename STATE, typename BAROSTAT, typename SCALE>
	void cellComputeDPotential_device(CELL cData, STATE state, BAROSTAT bState, SCALE scale)
	{
		uint numBlocks=0;
		uint numThreads=0;
		uint blockSize=128;
		computeGridSize(state.nParticles, blockSize, numBlocks, numThreads);
		
		//emptyCells_kernel<<<4096,numThreads>>>(cData);
		emptyCells_device(cData);
		CUDA_API_Errors(cudaDeviceSynchronize());
		calcHash_kernel<<<numBlocks,numThreads>>>(cData,state.p);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		sortParticles_device(cData);
		CUDA_API_Errors(cudaDeviceSynchronize());
		
		uint smemSize=sizeof(uint)*(numThreads+1);
		reorderDataAndFindCellStart_kernel<<<numBlocks,numThreads,smemSize>>>(cData,state.p,state.v);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		collideDPotential_kernel<<<numBlocks,numThreads>>>(cData,state.NBPconstants,
								   state.nTypes,scale);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
		//I'm going to reuse "cData.potentialEnergy" since it is overwritten everytime anyways
		vRedistributionByParticle_kernel<<<numBlocks,numThreads>>>(cData.gridParticleIndex, 
									   cData.potentialEnergy, 
									   bState.dPotential, 
									   state.nParticles);
		CUDA_Kernel_Errors();
		CUDA_API_Errors(cudaDeviceSynchronize());
	}
}
#endif
