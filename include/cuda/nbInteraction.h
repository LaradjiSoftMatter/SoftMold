#include <vector>
#include <cstring>
#include "dataTypes.h"
#include "errors.h"
#ifndef MPD_NBINTERACTION
#define MPD_NBINTERACTION

namespace mpd {
	
	//
	template <typename T, int nC>
	struct nbInteraction {
		nbInteraction(int nP, int nT, T *c):
			constants_d(NULL),constants_h(nT*nT*nC,0),elements_d(NULL),
			aTemp_h(nP),aTemp_d(NULL),_nTypes(nT),nbMask_h(nP,0),nbMask_d(NULL)
		{	
			for(int i=0;i<nConstants();i++)
				constants_h[i]=c[i];
		}
		
		~nbInteraction()
		{
			if(elements_d!=NULL) CUDA_API_Warnings(cudaFree(elements_d));
			if(constants_d!=NULL) CUDA_API_Warnings(cudaFree(constants_d));
			if(aTemp_d!=NULL) CUDA_API_Warnings(cudaFree(aTemp_d));
			if(nbMask_d!=NULL) CUDA_API_Warnings(cudaFree(nbMask_d));
		}
		
		void toDevice()
		{
			if(elements_d==NULL)
			{
				CUDA_API_Errors(cudaMalloc((void **)&elements_d, elements_h.size()*sizeof(int)));
				CUDA_API_Errors(cudaDeviceSynchronize());
			}
			if(constants_d==NULL)
			{
				CUDA_API_Warnings(cudaMalloc((void **)&constants_d, nConstants()*sizeof(T)));
				CUDA_API_Errors(cudaDeviceSynchronize());
			}
			if(aTemp_d==NULL)
			{
				CUDA_API_Errors(cudaMalloc((void **)&aTemp_d, aTemp_h.size()*sizeof(threeVector<T>)));
				CUDA_API_Errors(cudaDeviceSynchronize());
			}
			if(nbMask_d==NULL)
			{
				CUDA_API_Errors(cudaMalloc((void **)&nbMask_d, nbMask_h.size()*sizeof(int)));
				CUDA_API_Errors(cudaDeviceSynchronize());
			}
			
			CUDA_API_Errors(cudaMemcpy((char*)constants_d, (char*)constants_h.data(), 
				constants_h.size()*sizeof(T), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			CUDA_API_Errors(cudaMemcpy((char *) elements_d, (char*)elements_h.data(), 
				elements_h.size()*sizeof(int), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			CUDA_API_Errors(cudaMemcpy((char *) aTemp_d, (char*)aTemp_h.data(), 
				aTemp_h.size()*sizeof(threeVector<T>), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			CUDA_API_Errors(cudaMemcpy((char *) nbMask_d, (char*)nbMask_h.data(), 
				nbMask_h.size()*sizeof(int), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		void toHost()
		{
			CUDA_API_Errors(cudaMemcpy((char *) elements_h.data(), (char*)elements_d, 
				elements_h.size()*sizeof(int), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			
			CUDA_API_Errors(cudaMemcpy((char *) aTemp_h.data(), (char*)aTemp_d, 
				aTemp_h.size()*sizeof(threeVector<T>), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			//skipping constants
		}
		
		void addInteraction(int a)
		{
			elements_h.push_back(a);
			nbMask_h[a]=1;
			last++;
		}
		
		struct copyNBInteraction;
		
		copyNBInteraction deviceInteraction()
		{
			return copyNBInteraction(last, _nTypes, elements_d,constants_d,aTemp_d,nbMask_d);
		}
		
		copyNBInteraction hostInteraction()
		{
			return copyNBInteraction(last,_nTypes,elements_h.data(),constants_h.data(),aTemp_h.data(),nbMask_h.data());
		}
		
		constexpr int size() const {return elements_h.size();}
		constexpr int nInteractions() const {return 1;}
		constexpr int nConstants() {return _nTypes*_nTypes*nC;}
		
		//host vector
		std::vector<int> elements_h;
		std::vector<threeVector<T>> aTemp_h;
		std::vector<int> nbMask_h;
		std::vector<T> constants_h;
		
		//device pointers
		int *elements_d;
		threeVector<T> *aTemp_d;
		T *constants_d;
		int *nbMask_d;
		
		//current last element
		int last;
		
		int _nTypes;
		
		struct copyNBInteraction {
			copyNBInteraction(int n, int nT, int *e, T *c, threeVector<T> *aT, int* nbM):
					last(n),_nTypes(nT),elements(e),constants(c),aTemp(aT),nbMask(nbM)
			{}
			
			__host__ __device__
			int endIndex() const {return -1;}
			__host__ __device__
			int size() const {return last;}
			__host__ __device__
			int nInteractions() const {return 1;}
			__host__ __device__
			int nConstants() {return _nTypes*_nTypes*nC;}
			
			int last;
			int _nTypes;
			int *elements;
			T *constants;
			threeVector<T> *aTemp;
			int *nbMask;
		};
	};
	
}

#endif
