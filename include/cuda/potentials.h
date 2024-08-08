#include <cuda.h>
#include <cmath>
#include "dataTypes.h"
#ifndef MPD_POTENTIALS
#define MPD_POTENTIALS

namespace mpd {
	#define nNonBondedConst 6
	
	template <typename T>
	__host__ __device__
	threeVector<T> nonBondedF(threeVector<T> d, T cutoffSqr, T *constants,
				  uint type1, uint type2, uint nTypes)
	{
		threeVector<T> force(0.0f,0.0f,0.0f);// = make_float3(0.0f);
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		if (dr < cutoffSqr)
		{
			dr=sqrt(dr);
			
			T *c=&constants[nNonBondedConst*((type1*nTypes)+type2)];
			int offset=int(dr/c[0])*nNonBondedConst/2;
			
			T magnitude=c[offset]-dr;
			magnitude=((c[offset+1]-c[offset+2]*magnitude)*magnitude)/dr;
			force.x=d.x*magnitude;
			force.y=d.y*magnitude;
			force.z=d.z*magnitude;
		}
		return force;
	}
	
	template <typename T>
	__host__ __device__
	T nonBondedP(threeVector<T> d, T cutoffSqr, T *constants,
				  uint type1, uint type2, uint nTypes)
	{
		T potential=0;
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		if (dr < cutoffSqr)
		{
			dr=sqrt(dr);
			
			T *c=&constants[nNonBondedConst*((type1*nTypes)+type2)];
			if(dr<=c[0])
			{
				potential=c[0]-dr;
				potential=c[1]*potential*potential+c[2];
			}
			else
			{
				potential=c[3]-dr;
				potential=potential*potential*(c[4]-potential*c[5]);
			}
		}
		return potential;
	}
	
	template <typename T>
	__host__ __device__
	threeVector<T> harmonicF(threeVector<T> d, const int &offset, T *constants)
	{
		T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
		
		T magnitude=dr-constants[0];//1 flop
		//if(dr==0)
		//	std::cerr << d.x << ' ' << d.y << ' ' << d.z;
		magnitude=-magnitude*constants[1]/dr;//2 flops
		if(offset==1) magnitude*=-1.0;
		
		d.x*=magnitude;
		d.y*=magnitude;
		d.z*=magnitude;
		
		return d;
	}
	
	template <typename T>
	__host__ __device__
	T harmonicP(const threeVector<T> &d, const int& offset, T *constants)
	{
		T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
		
		T potential=dr-constants[0];
		potential=0.5*constants[1]*potential*potential;
		
		return potential;
	}
	
	
	template <typename T>
	__host__ __device__
	threeVector<T> bendF(threeVector<T> da, threeVector<T> db, const int &offset, T *constants)
	{
		threeVector<T> f;
		T dra=sqrt(da.x*da.x+da.y*da.y+da.z*da.z);
		T drb=sqrt(db.x*db.x+db.y*db.y+db.z*db.z);
		
		da.x/=dra;
		da.y/=dra;
		da.z/=dra;
		
		db.x/=drb;
		db.y/=drb;
		db.z/=drb;
		
		T costheta=(da.x*db.x)+(da.y*db.y)+(da.z*db.z);
		
		T magnitude=constants[0]-costheta;
		magnitude*=constants[1];
		
		T bufa=magnitude*(db.x-(da.x*costheta))/dra;
		T bufb=magnitude*(da.x-(db.x*costheta))/drb;
		switch(offset){
			case 0: f.x=bufa; break;
			case 1: f.x=bufb-bufa; break;
			case 2: f.x=-bufb; break;
		}
		bufa=magnitude*(db.y-(da.y*costheta))/dra;
		bufb=magnitude*(da.y-(db.y*costheta))/drb;
		switch(offset){
			case 0: f.y=bufa; break;
			case 1: f.y=bufb-bufa; break;
			case 2: f.y=-bufb; break;
		}
		bufa=magnitude*(db.z-(da.z*costheta))/dra;
		bufb=magnitude*(da.z-(db.z*costheta))/drb;
		switch(offset){
			case 0: f.z=bufa; break;
			case 1: f.z=bufb-bufa; break;
			case 2: f.z=-bufb; break;
		}
		
		return f;
	}
	
	template <typename T>
	__host__ __device__
	T bendP(threeVector<T> da, threeVector<T> db, const int& offset, T *constants)
	{
		T dra=sqrt(da.x*da.x+da.y*da.y+da.z*da.z);
		T drb=sqrt(db.x*db.x+db.y*db.y+db.z*db.z);
		
		da.x/=dra;
		da.y/=dra;
		da.z/=dra;
		
		db.x/=drb;
		db.y/=drb;
		db.z/=drb;
		
		T costheta=(da.x*db.x)+(da.y*db.y)+(da.z*db.z);
		
		T potential=constants[0]-costheta;
		potential=constants[1]*potential*potential*0.5;
		
		return potential;
	}
	
	template <typename T>
	threeVector<T> floatingBaseForce(threeVector<T> d, T *constants, threeVector<T> s)
	{
		threeVector<T> a=0;
		if(d.z<=constants[0])
		{
			a.z+=constants[3]*d.z*d.z*d.z-2*constants[3]*constants[0]*d.z*d.z+
				(2*constants[2]+constants[3]*constants[0]*constants[0])*d.z;
		}
		else if(d.z<constants[1])
		{
			T rcz=d.z-constants[1];
			a.z+=constants[5]*(rcz*rcz)*(2*d.z*d.z-(4*constants[1]-7*constants[0])*d.z+
			2*constants[1]*constants[1]+constants[0]*(-7*constants[1]+6*constants[0]));
		}
		return a;
	}
	
	
	template <typename T>
	T floatingBasePotential(threeVector<T> d, T *constants, threeVector<T> s)
	{
		T potential=0;
		if(d.z<=constants[0])
		{
			T rmz=constants[0]-d.z;
			potential=constants[2]*(constants[0]*constants[0]-d.z*d.z)+
				constants[3]*rmz*rmz*rmz*((constants[0]/3.0)-(1.0/4.0)*rmz)+constants[4];
		}
		else if(d.z<constants[1])
		{
			T rcz=constants[1]-d.z;
			potential=constants[5]*rcz*rcz*rcz*
				(rcz*((2.0/5.0)*rcz-(7.0/4.0)*constants[0])+(2.0)*constants[0]*constants[0]);
		}
		return potential;
	}
	
	//Constants (same as harmonic):
	//constants[0]=r_0
	//constants[1]=k
	//Derived from:
	//r_0=R_n+1 where R_n is radius of NP
	//k=2*k_R/(R_n-r_0)^2 where k_R is the maximum potential at surface of NP
	//And if r>r_0, potential=0
	template <typename T>
	__host__ __device__
	threeVector<T> harmonicHalfF(threeVector<T> d, const int &offset, T *constants)
	{
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		
		if(dr<constants[0]*constants[0])
		{
			dr=sqrt(dr);
			T magnitude=dr-constants[0];//1 flop
			//magnitude=magnitude*constants[1]/dr;//2 flops
			magnitude=-magnitude*constants[1]/dr;//2 flops
			
			//if(offset==1) magnitude*=-1.0;
			d.x*=magnitude;
			d.y*=magnitude;
			d.z*=magnitude;
		}
		else
		{
			d.x=0;
			d.y=0;
			d.z=0;
		}
		return d;
	}
	
	template <typename T>
	__host__ __device__
	T harmonicHalfP(threeVector<T> d, const int &offset, T *constants)
	{
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		T potential=0;
		
		if(dr<constants[0]*constants[0])
		{
			dr=sqrt(dr);
			potential=dr-constants[0];
			potential=0.5*constants[1]*potential*potential;
		}
		return potential;
	}
}

#endif
