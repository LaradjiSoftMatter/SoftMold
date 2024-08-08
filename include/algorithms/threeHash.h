#ifndef THREEHASH
#define THREEHASH

#include "dataTypes.h"

template<typename T, typename I=int>
class threeHash {
	public:
		threeHash(){
			s.x=0;
			s.y=0;
			s.z=0;
			c.x=0;
			c.y=0;
			c.z=0;
			n.x=0;
			n.y=0;
			n.z=0;
		};
		
		/** \brief
		 *  size defines boundaries' size
		 *  cellSize defines the minimum cell size, resulting
		 *  cellSize is the smallest integral splitting cell size.
		 **/
		threeHash(T size, T cellSize)
		{
			s=size;
			n.x=static_cast<I>(s.x/cellSize.x);
			n.y=static_cast<I>(s.y/cellSize.y);
			n.z=static_cast<I>(s.z/cellSize.z);
			c.x=s.x/(n.x);
			c.y=s.y/(n.y);
			c.z=s.z/(n.z);
			maxKey=n.x*n.y*n.z-1;
		};
		
		template <typename U>
		I operator () (U toHash)
		{
			threeVector<I> ci;
			ci.x=static_cast<I>((toHash.x)/c.x)%n.x;
			ci.y=static_cast<I>((toHash.y)/c.y)%n.y;
			ci.z=static_cast<I>((toHash.z)/c.z)%n.z;
			return this->hash(ci);
		}
		
		inline I hash(threeVector<I> toHash)
		{
			return (toHash.x+2*n.x)%n.x+((toHash.y+2*n.y)%n.y)*n.x+((toHash.z+2*n.z)%n.z)*n.x*n.y;
		}
		
		inline T minImg(threeVector<I> toHash)
		{
			T mI;
			mI.x=(toHash.x>=n.x)?s.x:0;
			mI.y=(toHash.y>=n.y)?s.y:0;
			mI.z=(toHash.z>=n.z)?s.z:0;
			mI.x=(toHash.x<0)?-s.x:mI.x;
			mI.y=(toHash.y<0)?-s.y:mI.y;
			mI.z=(toHash.z<0)?-s.z:mI.z;
			return mI;
		}
		
		inline threeVector<I> unHash(I toUnhash)
		{
			threeVector<I> ci;
			ci.x=static_cast<I>(toUnhash)%n.x;
			ci.y=static_cast<I>(toUnhash/n.x)%n.y;
			ci.z=static_cast<I>(toUnhash/(n.x*n.y))%n.z;
			
			return ci;
		}
		
		
		T s,c;
		threeVector<I> n;
		I maxKey;
};

#endif
