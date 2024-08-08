
#include <cuda.h>
#include "other.h"
#ifndef MPD_BAROSTAT
#define MPD_BAROSTAT

namespace mpd {
	
	template <typename T>
	struct barostat {
		using value_type=T;
		barostat(threeVector<T> dL, int nP):
			_deltaL(dL),_nParticles(nP),_dPotential_h(nP),_volume(0.0),
			_dist(std::uniform_real_distribution<T>(-dL.x,dL.x),
				std::uniform_real_distribution<T>(-dL.y,dL.y),
				std::uniform_real_distribution<T>(-dL.z,dL.z))
		{
			CUDA_API_Errors(cudaMalloc((void **)&_dPotential_d,nP*sizeof(T)));
			CUDA_API_Errors(cudaDeviceSynchronize());
		}
		
		~barostat()
		{
			if(_dPotential_d!=NULL) CUDA_API_Warnings(cudaFree(_dPotential_d));
		}
		
		//Independent in all 3 directions
		//WARNING: could blow up or collapse volume
		template <typename RANDGEN>
		threeVector<T> fluctuation(RANDGEN &randGen)
		{
			return threeVector<T>(_dist.x(randGen),_dist.y(randGen),_dist.z(randGen));
		}
		
		//independent in 2 directions, 3rd is constrained to maintain volume
		//WARNING: could "linearize" the system, possibly rectangular
		template <typename RANDGEN>
		threeVector<T> independentFluctuationConstantVolume(RANDGEN &randGen, threeVector<T> size)
		{
			if(_volume==0.0)
				_volume=size.x*size.y*size.z;
			threeVector<T> f(0.0,0.0,0.0);
			if(_deltaL.z==0 && _deltaL.x!=0 && _deltaL.y!=0)
			{
				f.x=_dist.x(randGen);
				f.y=_dist.y(randGen);
				//f.z=(size.x*size.y)/((size.x+f.x)*(size.y+f.y));
				f.z=(_volume)/((size.x+f.x)*(size.y+f.y)*size.z);
				f.z=size.z*f.z;
				f.z-=size.z;
			}
			if(_deltaL.x==0 && _deltaL.y!=0 && _deltaL.z!=0)
			{
				f.y=_dist.y(randGen);
				f.z=_dist.z(randGen);
				f.x=(size.z*size.y)/((size.z+f.z)*(size.y+f.y));
				f.x=size.x*f.x;
				f.x-=size.x;
			}
			if(_deltaL.y==0 && _deltaL.x!=0 && _deltaL.z!=0)
			{
				f.x=_dist.x(randGen);
				f.z=_dist.z(randGen);
				f.y=(size.z*size.x)/((size.z+f.z)*(size.x+f.x));
				f.y=size.y*f.y;
				f.y-=size.y;
			}
			return f;
		}
		
		//2 dimensions have the same fluctiuation, 1 dimension is constrained to maintain volume
		//WARNING: could "linearize" the system, only proportional dimensions
		template <typename RANDGEN>
		threeVector<T> dependentFluctuationConstantVolume(RANDGEN &randGen, threeVector<T> size)
		{
			//constrain dimension to maintain volume, but vary projected area
			//Note to anyone who sees this later, this method is very anistropic.
			// If you have any issue, try lowering deltaLXY or exclude any new types
			// from the interaction.
			threeVector<T> f(0.0,0.0,0.0);
			if(_deltaL.z==0 && _deltaL.x!=0 && _deltaL.y!=0)
			{
				f.x=_dist.x(randGen);
				f.y=f.x;
				f.z=(size.x*size.y)/((size.x+f.x)*(size.y+f.y));
				f.z=size.z*f.z;
				f.z-=size.z;
			}
			if(_deltaL.x==0 && _deltaL.y!=0 && _deltaL.z!=0)
			{
				f.y=_dist.y(randGen);
				f.z=f.y;
				f.x=(size.z*size.y)/((size.z+f.z)*(size.y+f.y));
				f.x=size.x*f.x;
				f.x-=size.x;
			}
			if(_deltaL.y==0 && _deltaL.x!=0 && _deltaL.z!=0)
			{
				f.x=_dist.x(randGen);
				f.z=f.x;
				f.y=(size.z*size.x)/((size.z+f.z)*(size.x+f.x));
				f.y=size.y*f.y;
				f.y-=size.y;
			}
			return f;
		}
		
		//dPotential=old-new
		template <typename RANDGEN>
		bool MCtest(T dPotentialTotal, T temperature, T pressure, T tension, threeVector<T> oldSize, 
			    threeVector<T> size, RANDGEN &randGen)
		{
			T D=std::exp(dPotentialTotal/temperature);
			if(tension!=0)
			{
				if(_deltaL.z==0 && _deltaL.x!=0 && _deltaL.y!=0)
					dPotentialTotal+=tension*((size.x*size.y)-(oldSize.x*oldSize.y));
				if(_deltaL.x==0 && _deltaL.y!=0 && _deltaL.z!=0)
					dPotentialTotal+=tension*((size.z*size.y)-(oldSize.z*oldSize.y));
				if(_deltaL.y==0 && _deltaL.x!=0 && _deltaL.z!=0)
					dPotentialTotal+=tension*((size.x*size.z)-(oldSize.x*oldSize.z));
			}
			if(pressure!=0)
				dPotentialTotal+=pressure*((size.x*size.y*size.z)-(oldSize.x*oldSize.y*oldSize.z));
			std::uniform_real_distribution<T> dist(0.0,1.0);
			//D is flippd for whatever reason, mine is old-new which is -(new-old)
			if(D>=dist(randGen) || -dPotentialTotal<=0)
				return true;
			return false;
		}
		
		struct copyState;
		
		copyState deviceState()
		{
			return copyState(_dPotential_d, _deltaL, _nParticles);
		}
		
		copyState host()
		{
			return copyState(_dPotential_h.data(), _deltaL, _nParticles);
		}
		
		struct copyState {
			copyState(T *dP, threeVector<T> dL, int nP):
				  dPotential(dP),deltaL(dL),nParticles(nP)
			{}
			T *dPotential;
			threeVector<T> deltaL;
			int nParticles;
		};
		
		private:
			std::vector<T> _dPotential_h;
			T *_dPotential_d,_tensionL, _pressureL;
			threeVector<T> _deltaL;
			int _nParticles;
			threeVector<std::uniform_real_distribution<T>> _dist;
			T _volume;
	};
	
}

#endif
