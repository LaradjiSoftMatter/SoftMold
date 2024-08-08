#include "../algorithms/dataTypes.h"

#ifndef MPD_DATATYPES
#define MPD_DATATYPES

namespace mpd {
	/*
	template <typename T>
	struct threeVector
	{
		threeVector()=default;
		__device__ __host__
		threeVector(const T &a, const T &b, const T &c):x(a),y(b),z(c){}
		union {
			struct {T x,y,z;};
			T s[3];
		};
		
		__device__ __host__
		threeVector<T>& operator+=(const threeVector<T> &b)
		{
			x+=b.x;
			y+=b.y;
			z+=b.z;
			return *this;
		}
	};
	
	template <typename T>
	struct position
	{
		position()=default;
		position(const T &a, const T &b, const T &c, const int &d):x(a),y(b),z(c),type(d){}
		union {
			struct {T x,y,z;};
			T s[3];
		};
		int type;
	};
	*/
	
	template <typename T>
	struct tvPlus {
	__host__ __device__
	threeVector<T> operator()(const threeVector<T> &a, const threeVector<T> &b)
	{
		threeVector<T> result;
		result.x=a.x+b.x;
		result.y=a.y+b.y;
		result.z=a.z+b.z;
		return result;
	}
	};
	
	template <typename T>
	__host__ __device__
	threeVector<T> difference(const position<T> &pA, const position<T> pB)
	{
		return threeVector<T>(
			pA.x-pB.x,
			pA.y-pB.y,
			pA.z-pB.z);
	}
	
	template <typename T>
	__host__ __device__
	threeVector<T> minImg(threeVector<T> d, const threeVector<T> &size)
	{
		if(d.x>=size.x/2.0)d.x-=size.x;
		if(d.y>=size.y/2.0)d.y-=size.y;
		if(d.z>=size.z/2.0)d.z-=size.z;
		if(d.x<=-size.x/2.0)d.x+=size.x;
		if(d.y<=-size.y/2.0)d.y+=size.y;
		if(d.z<=-size.z/2.0)d.z+=size.z;
		return d;
	}
		
	template <typename T>
	__host__ __device__
	position<T> bounds(position<T> p, const threeVector<T> &size)
	{
		if(p.x>=size.x) p.x-=size.x;
		if(p.x<0.0) p.x+=size.x;
		if(p.y>=size.y) p.y-=size.y;
		if(p.y<0.0) p.y+=size.y;
		if(p.z>=size.z) p.z-=size.z;
		if(p.z<0.0) p.z+=size.z;
		return p;
	}
}

#endif	
