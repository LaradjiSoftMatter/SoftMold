//algorithm specific includes

#define ERRORS_ENABLED

//dataTypes.h and functions.h are included in all of these anyway 
#include "algorithms/cellOpt.h"
#include "algorithms/verlet.h"

//this is for some lipid specific data collection
//#include "lipids/lipids.h"
#include "fileFormats/xyzFormat.h"

//Common type definitions
#define WALL 1
#define SOLVENT 2
#define MEMBRANE 3

//these are the number of constants used by the Force and Potential functions
//These are for MD

//These are pretty straight forward:
// threeVector d, da, and db are difference vectors (a vector normal to the conservative force between point 1 and 2);
// a1 is the acceleration/force along the difference vector(s);
// a2 is the acceleration/force against the difference vector(s);
// when a3 is present (three-body), a1+a3 are opposite of a2, and any other combination;
// there might eventually be a a2 less one, but I don't care right now.
//They allow you to do your own calculation without the need for a specific particle.
//pair potentials and forces as described by revalee, et. al., forces derived by spangler
//floating point operations (flops) only count required operations

template <typename T>
inline T kernal(T dr, T h)
{
	//T dx=p1.x-p2.x;//1 flop
	//T dy=p1.y-p2.y;//1 flop
	//T dz=p1.z-p2.z;//1 flop
	
	//T dr=dx*dx+dy*dy+dz*dz;//5 flops
	
	T hSqr=h*h;
	dr/=hSqr;	

	T weight=0;
	
	if(dr<4.0 && dr>=1.0)
	{
		dr=sqrt(dr);
		weight=(2.0-dr);
		weight=weight*weight*weight*0.25;
	}
	else if(dr<1.0)
	{
		T drSqrt=sqrt(dr);
		weight=1.0-1.5*dr+0.75*dr*drSqrt;
	}
	return weight/(M_PI*h*hSqr);
}

template <typename T>
inline threeVector<T> gradKernal(threeVector<T> d, T dr, T h)
{
	//T dx=p1.x-p2.x;//1 flop
	//T dy=p1.y-p2.y;//1 flop
	//T dz=p1.z-p2.z;//1 flop
	
	//T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
	
	T hSqr=h*h;
	
	dr/=hSqr;
	
	T coef=0;
	
	if(dr<4.0 && dr>=1.0)
	{
		T drSqrt=sqrt(dr);
		coef=(2.0-dr);
		coef=-coef*coef*0.5/(M_PI*h*hSqr);
	}
	else if(dr<1.0)
	{
		T drSqrt=sqrt(dr);
		coef=-3.0*drSqrt+2.25*dr/(M_PI*h*hSqr);
	}
	
	threeVector<T> gradVector;
	gradVector.x=coef*d.x/dr;
	gradVector.y=coef*d.y/dr;
	gradVector.z=coef*d.z/dr;
}

template <typename T>
inline T partialDensity(T mass, T W)
{
	return mass*W;
}

template <typename T>
inline T pressure(T rho, T rhoZero)
{
	T buf=rho/rhoZero;
	buf*=buf*buf;
	return rhoZero*(buf-1.0);
}

template <typename T>
inline threeVector<T> partialPressure(T mass2, T rho1, T rho2, T P1, T P2, threeVector<T> gradVector)
{
	rho1*=rho1;
	rho2*=rho2;
	T magnitude=mass2*(P1/rho1+P2/rho2);
	threeVector<T> result;
	result.x=magnitude*gradVector.x;
	result.y=magnitude*gradVector.y;
	result.z=magnitude*gradVector.z;
	return result;
}

template <typename T>
inline threeVector<T> partialViscosity(T mass2, T mu1, T mu2, threeVector<T> d, threeVector<T> dv, T dr, threeVector<T> gradVector

template <typename T>
inline threeVector<T> force(T g, T P, T rho, )
{
	threeVector<T> d,dv;
	d.x=p1.x-p2.x;
	d.y=p1.y-p2.y;
	d.z=p1.z-p2.z;
	
	dv.x=v1.x-v2.x;
	dv.y=v1.y-v2.y;
	dv.z=v1.z-v2.z;
	
	T dr=dx*dx+dy*dy+dz*dz;//5 flops
	
	P
}
 
