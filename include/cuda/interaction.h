#include <vector>
#include <cstring>
#include "dataTypes.h"
#include "errors.h"
#ifndef MPD_INTERACTION
#define MPD_INTERACTION

namespace mpd {
	
	template <typename T,int nI, int nC>
	struct interactionElement {
		__host__ __device__
		interactionElement()
		{
			nextOffset=-1;
			for(int i=0;i<nI;i++)
				index[i]=-1;
			for(int i=0;i<nC;i++)
				constants[i]=0;
		}
		int nextOffset;
		T constants[nC];
		int index[nI];
		int offset;
		__host__ __device__
		int nInteractions() const {return nI;}
		__host__ __device__
		int nConstants() const {return nC;}
		__host__ __device__
		int endIndex() const {return -1;}
	};
		
	//int nI;//2==BOND,  3==BEND, 18==BALL
	//int nC;//per interaction
	template <typename T, int nI, int nC>
	struct interaction {
		interaction(int nP):
			last(nP),elements_h(nP),elements_d(NULL),
			aTemp_h(nP),aTemp_d(NULL)
		{}	
		
		~interaction()
		{
			//if(elements_d!=NULL) delete elements_d;
			if(elements_d!=NULL) CUDA_API_Warnings(cudaFree(elements_d));
			if(aTemp_d!=NULL) CUDA_API_Warnings(cudaFree(aTemp_d));
		}
		
		void toDevice()
		{
			//elements_d=new interactionElement<T,nI,nC>[elements_h.size()];
			//memcpy((char *) elements_d,(char *)elements_h.data(),elements_h.size()*sizeof(interactionElement<T,nI,nC>));
			CUDA_API_Errors(cudaMalloc((void **)&elements_d, elements_h.size()*sizeof(interactionElement<T,nI,nC>)));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) elements_d, (char*)elements_h.data(), 
				elements_h.size()*sizeof(interactionElement<T,nI,nC>), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMalloc((void **)&aTemp_d, aTemp_h.size()*sizeof(threeVector<T>)));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) aTemp_d, (char*)aTemp_h.data(), 
				aTemp_h.size()*sizeof(threeVector<T>), cudaMemcpyHostToDevice));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		void toHost()
		{
			//memcpy((char *) elements_h.data(),(char *)elements_d,elements_h.size()*sizeof(interactionElement<T,nI,nC>));
			CUDA_API_Errors(cudaMemcpy((char *) elements_h.data(), (char*)elements_d, 
				elements_h.size()*sizeof(interactionElement<T,nI,nC>), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
			CUDA_API_Errors(cudaMemcpy((char *) aTemp_h.data(), (char*)aTemp_d, 
				aTemp_h.size()*sizeof(threeVector<T>), cudaMemcpyDeviceToHost));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		void addInteraction(int a, int offset, int *b, T *c)
		{
			int i=a;
			while(i!=endIndex())
			{
				//push an element onto host vector,
				//using nextOffset to indicate an empty slot
				if(elements_h[i].nextOffset==endIndex())
				{
					if(last==elements_h.size())
						elements_h.resize(elements_h.size()*2);
					elements_h[i].nextOffset=last;
					for(int j=0;j<nC;j++)
						elements_h[i].constants[j]=c[j];
					for(int j=0;j<nI;j++)
						elements_h[i].index[j]=b[j];
					elements_h[i].offset=offset;
					i=endIndex();
					center=b[0];
					last++;
				}
				else
					i=elements_h[i].nextOffset;
			}
		}
		
		struct copyInteraction;
		
		copyInteraction deviceInteraction()
		{
			return copyInteraction(last,center,elements_d,aTemp_d);
		}
		
		copyInteraction hostInteraction()
		{
			return copyInteraction(last,center,elements_h.data(),aTemp_h.data());
		}
		
		constexpr int endIndex() const {return -1;}
		constexpr int size() const {return elements_h.size();}
		constexpr int nInteractions() const {return nI;}
		constexpr int nConstants() const {return nC;}
		
		//host vector
		std::vector<interactionElement<T,nI,nC>> elements_h;
		std::vector<threeVector<T>> aTemp_h;
		
		//device pointers
		interactionElement<T,nI,nC> *elements_d;
		threeVector<T> *aTemp_d;
		
		//current last element
		int last,center;
		
		struct copyInteraction {
			copyInteraction(int n, int c, interactionElement<T,nI,nC> *e, threeVector<T> *aT):
					last(n),center(c),elements(e),aTemp(aT)
			{}
			__host__ __device__
			int endIndex() const {return -1;}
			__host__ __device__
			int size() const {return last;}
			__host__ __device__
			int nInteractions() const {return nI;}
			__host__ __device__
			int nConstants() const {return nC;}
			
			int last,center;
			interactionElement<T,nI,nC> *elements;
			threeVector<T> *aTemp;
		};
	};
	
}

#endif
