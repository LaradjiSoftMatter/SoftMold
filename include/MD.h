//These really need their own class/struct per potential/force combo
//Also, these #define constants (as const) need to be part of the class/struct
//algorithm specific includes


//dataTypes.h and functions.h are included in all of these anyway 
#include "algorithms/cell.h"
#include "algorithms/cellOpt.h"
#include "algorithms/dataCollection.h"
#include "algorithms/verlet.h"
#include "algorithms/langevin.h"
#include "algorithms/molecules.h"
#include "algorithms/volume.h"

//this is for some lipid specific data collection
//#include "lipids/lipids.h"
#include "fileFormats/xyzFormat.h"

//other types, note that SOLVENT and SOLVENT_FLAG are not the same...
#define MONOMER 1
#define HEAD 2
#define TAIL 3
#define CYTO 4
#define ANCHOR 5
#define HEAD_ANCHOR 6
#define TAIL_ANCHOR 7
#define SOLVENT 0
//#define SOLVENT_FLAG 0
//#define SOLVENT_TRACER 1

#define TAILA 1
#define HEADA 2
#define TAILB 3
#define HEADB 4

//these are the number of constants used by the Force and Potential functions
//These are for MD
#define nTWOBODYFCONST 6
#define nTWOBODYMFCONST 2
#define nTWOBODYUCONST 6
#define nTWOBODYMUCONST 2
#define nTHREEBODYMFCONST 1
#define nTHREEBODYMUCONST 1
#define nBOUNDARYCONST 4
#define nOFFSET_BOUNDARYCONST 4
#define nFLOATING_BASECONST 6
#define nZTORQUECONST 4
#define nZPOWERPOTENTIALCONST 2

//this is for BEAD type molecules, but requires nTypes like twoBodyFConst and twoBodyUConst
#define BEADBEADOFFSET 11
#define nBEADCONST 22
#define BEADRADIUS 4
//same constants
//#define nBEADUCONST 7

//These are for specialized molecules
#define nCHAINCONST 4
#define nBONDCONST 2
#define nBENDCONST 2
#define nRIGIDBENDCONST 5
#define nPULLBEADCONST 4

//These are offsets for chained types constants
#define CHAINBOND 0 
#define CHAINBEND 2

//These are just names
#define ABOND 0
#define KBOND 1

#define ABEND 0
#define KBEND 1

#define RIGIDBENDX 2
#define RIGIDBENDY 3
#define RIGIDBENDZ 4

#define EPSILON 6.0

#define sign(v) ((v<0)?-1.0:1.0)

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
void nonBondedF(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants, T cutoffSquared)
{
	//3 flops from the difference vector d
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		int offset=int(dr/(constants[0]))*3;
		T magnitude=constants[offset]-dr;//1 flop
		magnitude=((constants[offset+1]-constants[offset+2]*magnitude)*magnitude)/dr;//4 flops
		
		//Expanded:
		//T magnitude=0;
		//if(dr<constants[0])
		//{
		//	//r_m-r
		//	magnitude=constants[0]-dr;//1 flop
		//	//2*(U_max-U_min)/r_m^2*(r_m-r)
		//	magnitude=(constants[1]*magnitude)/dr;//4 flops
		//} 
		//else 
		//{
		//	//r_c-r
		//	magnitude=constants[3]-dr;//1 flop
		//	//(6*U_min/(r_c-r_m)^2-6*U_min/(r_c-r_m)^3*(r_c-r))*(r_c-r)/r
		//	magnitude=((constants[4]-constants[5]*magnitude)*magnitude)/dr;//4 flops
		//}
		
		d.x*=magnitude;//1 flop
		d.y*=magnitude;//1 flop
		d.z*=magnitude;//1 flop
		
		a1.x+=d.x;//1 flop
		a1.y+=d.y;//1 flop
		a1.z+=d.z;//1 flop
		//#pragma omp atomic
		a2.x-=d.x;//1 flop
		//#pragma omp atomic
		a2.y-=d.y;//1 flop
		//#pragma omp atomic
		a2.z-=d.z;//1 flop
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

template <typename T>
T nonBondedP(threeVector<T> d, T *constants, T cutoffSquared)
{
	T potential=0;
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		if(dr<=constants[0])
		{
			potential=constants[0]-dr;
			potential=constants[1]*potential*potential+constants[2];
		}
		else
		{
			potential=constants[3]-dr;
			potential=potential*potential*(constants[4]-potential*constants[5]);
		}
	}
	return potential;
}

template <typename T>
void beadForce(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants, T cutoffSquared)
{
	//3 flops from the difference vector d
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
	T magnitude=0;
	T rminD=constants[4]+constants[5];
	T rminD2=rminD*rminD;
	//is it in range?
	if(rminD2<=dr && dr<cutoffSquared)
	{
		dr=sqrt(dr);
		//constants[0]=rc+r
		//R+r_c-r
		T E=constants[0]-dr;//1 flops
		//E^2
		T E2=E*E;//1 flop
		//E^3
		T E3=E*E2;//1 flop
		//constants[1]=-7/4*rmin
		//constants[2]=2*rmin^2
		//-2*(0.4*E^3-7*E^2*rmin/4+2*rmin^2*E)
		magnitude=2.0*(0.4*E3+constants[1]*E2+constants[2]*E)/(dr);
		magnitude+=4.0*E2+8.0*constants[1]*E+6.0*constants[2];
		//magnitude=4.0*(0.4*E3+constants[1]*E2+constants[2]*E)/dr;//7 flops
		//(2*E^2-4*7*E/4+3*2*rmin^2)
		//magnitude+=4.0*E2+8.0*constants[1]*E+6.0*constants[2];//7 flops
		//constants[3]=M_PI*R*sigma/rmin^3
		//E^2*M_PI*R*sigma/(rmin^3*r^3)
		magnitude*=(E2*constants[3])/(dr*dr);//5 flops
	}
	//is it in range of inner range, r<R+rmin
	else if(dr<rminD2)
	{
		dr=sqrt(dr);
		//constants[4]=R
		//constants[5]=rmin
		//r-R-rmin, There seems to have been a sign error here
		T B=rminD-dr;
		T B2=B*B;
		T B3=B2*B;
		T B4=B2*B2;
		//T B=dr-constants[4]-constants[5];
		
		//constants[6]=-(Umax-Umin)*M_PI*R*sigma/(6*rmin^2)
		//constants[7]=Umin*M_PI*R*sigma
		magnitude=(B4*constants[6]+B3*constants[7]+B2*constants[8]+B*constants[9]+constants[10])/(dr);
		magnitude+=(4.0*B3*constants[6]+3.0*B2*constants[7]+2.0*B*constants[8]+constants[9]);
		magnitude/=dr*dr;
	}
	
	d.x*=magnitude;//1 flop
	d.y*=magnitude;//1 flop
	d.z*=magnitude;//1 flop
	
	a1.x+=d.x;//1 flop
	a1.y+=d.y;//1 flop
	a1.z+=d.z;//1 flop
	#pragma omp atomic
	a2.x-=d.x;//1 flop
	#pragma omp atomic
	a2.y-=d.y;//1 flop
	#pragma omp atomic
	a2.z-=d.z;//1 flop
	
	//41 "useful" flops total, compare to the normal n^2*22 solid pairwise flops
}

template <typename T>
T beadPotential(threeVector<T> d, T *constants, T cutoffSquared)
{
	T potential=0;
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	T rminD=constants[4]+constants[5];
	T rminD2=rminD*rminD;
	//is it in range?
	if(rminD2<=dr && dr<cutoffSquared)
	{
		dr=sqrt(dr);
		//constants[0]=rc+r
		//R+r_c-r
		T E=constants[0]-dr;//1 flops
		//E^2
		T E2=E*E;//1 flop
		//E^3
		T E3=E*E2;//1 flop
		
		//constants[1]=-7/4*rmin
		//constants[2]=2*rmin^2
		//constants[3]=M_PI*R*sigma/rmin^3
		potential=2.0*E3*(0.4*E2+constants[1]*E+constants[2])*constants[3]/(dr);
	}
	//is it in range of inner range, r<R+rmin
	else if(dr<rminD2)
	{
		dr=sqrt(dr);
		//constants[4]=R
		//constants[5]=rmin
		//r-R-rmin, There seems to have been a sign error here
		T B=rminD-dr;
		T B2=B*B;
		T B3=B2*B;
		T B4=B2*B2;
		//constants[6]=-(Umax-Umin)*M_PI*R*sigma/(6*rmin^2)
		//constants[7]=Umin*M_PI*R*sigma
		//potential=(B*B*((3.0*(constants[4]-dr)-constants[5])*constants[6]*B+constants[7])+constants[8])/(dr);
		
		potential=(B4*constants[6]+B3*constants[7]+B2*constants[8]+B*constants[9]+constants[10])/(dr);
	}
	return potential;
}

template <typename T>
void beadBeadForce(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants, T cutoffSquared)
{
	//3 flops from the difference vector d
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
	T magnitude=0;
	T rminD=constants[0+BEADBEADOFFSET];
	T rminD2=rminD*rminD;
	//is it in range?
	if(rminD2<=dr && dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		//2*R+r_c-r
		T E=constants[1+BEADBEADOFFSET]-dr;//1 flops
		//E^2
		T E2=E*E;//1 flop
		//E^3
		T E3=E*E2;//1 flop
		
		magnitude=E3*(6.0*constants[2+BEADBEADOFFSET]*E2+5.0*constants[3+BEADBEADOFFSET]*E+4.0*constants[4+BEADBEADOFFSET]+
		(constants[2+BEADBEADOFFSET]*E3+constants[3+BEADBEADOFFSET]*E2+constants[4+BEADBEADOFFSET]*E)/dr)/(dr*dr);
	}
	//is it in range of inner range, r<R+rmin
	else if(dr<rminD2)
	{
		dr=sqrt(dr);
		//constants[4]=R
		//constants[5]=rmin
		//r-R-rmin, There seems to have been a sign error here
		T B=rminD-dr;
		T B2=B*B;
		T B3=B2*B;
		T B4=B2*B2;
		T B5=B2*B3;
		//constants[6]=-(Umax-Umin)*M_PI*R*sigma/(6*rmin^2)
		//constants[7]=Umin*M_PI*R*sigma
		//potential=(B*B*((3.0*(constants[4]-dr)-constants[5])*constants[6]*B+constants[7])+constants[8])/(dr);
		
		magnitude=5.0*B4*constants[5+BEADBEADOFFSET]+4.0*B3*constants[6+BEADBEADOFFSET]+3.0*B2*constants[7+BEADBEADOFFSET]+
			   2.0*B*constants[8+BEADBEADOFFSET]+constants[9+BEADBEADOFFSET];
		magnitude+=(B5*constants[5+BEADBEADOFFSET]+B4*constants[6+BEADBEADOFFSET]+B3*constants[7+BEADBEADOFFSET]+
			   B2*constants[8+BEADBEADOFFSET]+B*constants[9+BEADBEADOFFSET]+constants[10+BEADBEADOFFSET])/(dr);
		magnitude/=dr*dr;
	}
	//std::cout << magnitude << '\t' << dr << std::endl;
	d.x*=magnitude;//1 flop
	d.y*=magnitude;//1 flop
	d.z*=magnitude;//1 flop
	
	a1.x+=d.x;//1 flop
	a1.y+=d.y;//1 flop
	a1.z+=d.z;//1 flop
	//#pragma omp atomic
	a2.x-=d.x;//1 flop
	//#pragma omp atomic
	a2.y-=d.y;//1 flop
	//#pragma omp atomic
	a2.z-=d.z;//1 flop
	
	//41 "useful" flops total, compare to the normal n^2*22 solid pairwise flops
}

template <typename T>
T beadBeadPotential(threeVector<T> d, T *constants, T cutoffSquared)
{
	T potential=0;
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	T rminD=constants[0+BEADBEADOFFSET];
	T rminD2=rminD*rminD;
	//is it in range?
	if(rminD2<=dr && dr<cutoffSquared)
	{
		dr=sqrt(dr);
		//2*R+r_c-r
		T E=constants[1+BEADBEADOFFSET]-dr;//1 flops
		//E^2
		T E2=E*E;//1 flop
		//E^3
		T E3=E*E2;//1 flop
		//E^4
		T E4=E2*E2;
		potential=E4*(constants[2+BEADBEADOFFSET]*E2+constants[3+BEADBEADOFFSET]*E+constants[4+BEADBEADOFFSET])/dr;
	}
	//is it in range of inner range, r<R+rmin
	
	else if(dr<rminD2)
	{
		dr=sqrt(dr);
		//constants[4]=R
		//constants[5]=rmin
		//r-R-rmin, There seems to have been a sign error here
		T B=rminD-dr;
		T B2=B*B;
		T B3=B2*B;
		T B4=B2*B2;
		T B5=B2*B3;
		//constants[6]=-(Umax-Umin)*M_PI*R*sigma/(6*rmin^2)
		//constants[7]=Umin*M_PI*R*sigma
		//potential=(B*B*((3.0*(constants[4]-dr)-constants[5])*constants[6]*B+constants[7])+constants[8])/(dr);
		
		potential=(B5*constants[5+BEADBEADOFFSET]+B4*constants[6+BEADBEADOFFSET]+B3*constants[7+BEADBEADOFFSET]+
			   B2*constants[8+BEADBEADOFFSET]+B*constants[9+BEADBEADOFFSET]+constants[10+BEADBEADOFFSET])/(dr);
	}
	return potential;
}

template <typename T>
void harmonicF(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants)
{
	T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	
	T magnitude=dr-constants[0];//1 flop
	magnitude=-magnitude*constants[1]/dr;//2 flops
	
	d.x*=magnitude;
	d.y*=magnitude;
	d.z*=magnitude;
	
	a1.x+=d.x;
	a1.y+=d.y;
	a1.z+=d.z;
	
	a2.x-=d.x;
	a2.y-=d.y;
	a2.z-=d.z;
}

template <typename T>
void harmonicF(threeVector<T> &d, threeVector<T> &a, T *constants)
{
	T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	
	T magnitude=dr-constants[0];//1 flop
	magnitude=-magnitude*constants[1]/dr;//2 flops
	
	d.x*=magnitude;
	d.y*=magnitude;
	d.z*=magnitude;
	
	a.x+=d.x;
	a.y+=d.y;
	a.z+=d.z;
}

template <typename T>
T harmonicP(threeVector<T> d, T *constants)
{
	T dr=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	
	T potential=dr-constants[0];
	potential=0.5*constants[1]*potential*potential;
	
	return potential;
}

template <typename T>
void harmonicFZ(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants)
{
	d.x*=constants[0];
	d.y*=constants[0];
	d.z*=constants[0];
	
	a1.x+=d.x;
	a1.y+=d.y;
	a1.z+=d.z;
	
	a2.x-=d.x;
	a2.y-=d.y;
	a2.z-=d.z;
}

template <typename T>
T harmonicPZ(threeVector<T> d, T *constants)
{
	
	T potential=0.5*constants[0]*(d.x*d.x+d.y*d.y+d.z*d.z);
	
	return potential;
}

template <typename T>
void floatingBaseForce(threeVector<T> d, threeVector<T> &a, T *constants, threeVector<T> s)
{
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



/** \brief Boundary force
 * 
 * 
 **/
template <typename T>
void offsetBoundaryF(position<T> p, threeVector<T> &a, T *constants, threeVector<T> s)
{
	//For U_boundary=(k/2)*(|z|-z_0)^2
	//c[0]==dim {0,1,2} cooresponding to {x,y,z}
	//c[1]==center point (scaled with size, s vector)
	//c[2]==cutoff
	//c[3]==k/2
	/*
	T magnitude=0;
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	
	T dir=sign(d);
	d=abs(d)-offset;
	if(d>constants[2] && d<2.0*constants[2])
	{
		d=((2.0*d-constants[2])-(2.0*constants[2]-d))*(2.0*constants[2]-d);
		magnitude=(2.0*constants[3]/pow(constants[2],3.0))*d;
	}
	*/
	//a.s[dim]+=magnitude*dir;
	T magnitude=0;
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	T dir=sign(d);//this is +/-1.0
	d=fabs(d)-constants[2];
	T d2=d*d;
	if(d2<=2.0)
	{
		magnitude=constants[3]*(4.0/(d2*d2*d)-2.0/(d2*d));
	}
	//std::cout << d << ' ' << magnitude << '\n';
	a.s[dim]+=magnitude*dir;
}

/** \brief Boundary force
 * 
 * 
 **/
template <typename T>
void boundaryF(position<T> p, threeVector<T> &a, T *constants, threeVector<T> s)
{
	//For U_boundary=(k/2)*(|z|-z_0)^2
	//c[0]==dim {0,1,2} cooresponding to {x,y,z}
	//c[1]==center point (scaled with size, s vector)
	//c[2]==cutoff
	//c[3]==k/2
	/*
	T magnitude=0;
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	T dir=sign(d);
	d=abs(d);
	if(d>constants[2] && d<2.0*constants[2])
	{
		d=((2.0*d-constants[2])-(2.0*constants[2]-d))*(2.0*constants[2]-d);
		magnitude=(2.0*constants[3]/pow(constants[2],3.0))*d;
	}
	*/
	//For U_boundary=k*(1/z^4-1/z^2) if z<=sqrt(2)
	//    U_boundary=0               if z>sqrt(2)
	T magnitude=0;
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	T dir=sign(d);//this is +/-1.0
	d=fabs(d);//constants[2] is not used currently
	T d2=d*d;
	if(d2<=2.0)
	{
		magnitude=constants[3]*(4.0/(d2*d2*d)-2.0/(d2*d));
	}
	
	a.s[dim]+=magnitude*dir;
}

/** \brief Boundary potential
 * 
 * 
 **/
template <typename T>
T boundaryP(position<T> p, T *constants, threeVector<T> s)
{
	//For U_boundary=(k/2)*(|z|-z_0)^2
	//c[0]==dim {0,1,2} cooresponding to {x,y,z}
	//c[1]==center point
	//c[2]==cutoff
	//c[3]==k/2
	/*
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	d=abs(d);
	T potential=0;
	if(d>0 && d<constants[2])
	{
		potential=constants[3];
	}
	if(d>constants[2] && d<2.0*constants[2])
	{
		d=(2.0*d-constants[2])*pow((2.0*constants[2]-d),2.0);
		potential=(constants[3]/pow(constants[2],3.0))*d;
	}
	*/
	//For U_boundary=k*(1/z^4-1/z^2) if z<=sqrt(2)
	//    U_boundary=0               if z>sqrt(2)
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	d=fabs(d);
	T potential=0;
	d=d*d;
	if(d<=2.0)
	{
		potential=constants[3]*(1.0/(d*d)-1.0/d);
	}
	return potential;
}

/** \brief Boundary potential
 * 
 * 
 **/
template <typename T>
T offsetBoundaryP(position<T> p, T *constants, threeVector<T> s)
{
	//For U_boundary=(k/2)*(|z|-z_0)^2
	//c[0]==dim {0,1,2} cooresponding to {x,y,z}
	//c[1]==center point
	//c[2]==cutoff
	//c[3]==k/2
	/*
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	
	d=abs(d)-offset;
	T potential=0;
	if(d>0 && d<constants[2])
	{
		potential=constants[3];
	}
	if(d>constants[2] && d<2.0*constants[2])
	{
		d=(2.0*d-constants[2])*pow((2.0*constants[2]-d),2.0);
		potential=(constants[3]/pow(constants[2],3.0))*d;
	}
	*/
	//For U_boundary=k*(1/z^4-1/z^2) if z<=sqrt(2)
	//    U_boundary=0               if z>sqrt(2)
	//T d=abs(z);
	//T potential=0;
	//d=d*d;
	//if(d<=2.0)
	//{
	//	potential=constants[3]*(1.0/(d*d)-1.0/d);
	//}
	int dim=static_cast<int>(constants[0]);
	T d=p.s[dim]-constants[1];
	d-=(d>s.s[dim]/2.0)?s.s[dim]:0;
	d+=(d<-s.s[dim]/2.0)?s.s[dim]:0;
	d=fabs(d)-constants[2];
	T potential=0;
	d=d*d;
	if(d<=2.0)
	{
		potential=constants[3]*(1.0/(d*d)-1.0/d);
	}
	return potential;
}

template <typename T>
void bendF(threeVector<T> da, threeVector<T> db, threeVector<T> &a1, 
		   threeVector<T> &a2,threeVector<T> &a3,T *constants)
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
	
	T magnitude=constants[0]-costheta;
	magnitude*=constants[1];
	
	T bufa=magnitude*(db.x-(da.x*costheta))/dra;
	T bufb=magnitude*(da.x-(db.x*costheta))/drb;
	
	a1.x+=bufa;
	a2.x+=(bufb-bufa);
	a3.x-=bufb;
	
	bufa=magnitude*(db.y-(da.y*costheta))/dra;
	bufb=magnitude*(da.y-(db.y*costheta))/drb;
	
	a1.y+=bufa;
	a2.y+=(bufb-bufa);
	a3.y-=bufb;
	
	bufa=magnitude*(db.z-(da.z*costheta))/dra;
	bufb=magnitude*(da.z-(db.z*costheta))/drb;
	
	a1.z+=bufa;
	a2.z+=(bufb-bufa);
	a3.z-=bufb;
}

template <typename T>
void bendF(threeVector<T> &da, threeVector<T> &db, threeVector<T> &a, T *constants, int &offset)
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
	
	T magnitude=constants[0]-costheta;
	magnitude*=constants[1];
	
	T tmpa=magnitude/dra;
	T tmpb=magnitude/drb;
	
	//we will see how fast this little lookup table is, there should be no warp divergence
	twoVector<T> signs[3];
	signs[0].s[0]=1;
	signs[0].s[1]=0;
	signs[1].s[0]=-1;
	signs[1].s[1]=1;
	signs[2].s[0]=0;
	signs[2].s[1]=-1;
	
	T bufa=tmpa*(db.s[0]-(da.s[0]*costheta));
	T bufb=tmpb*(da.s[0]-(db.s[0]*costheta));
	a.s[0]+=signs[offset].s[0]*bufa+signs[offset].s[1]*bufb;
	
	bufa=tmpa*(db.s[1]-(da.s[1]*costheta));
	bufb=tmpb*(da.s[1]-(db.s[1]*costheta));
	a.s[1]+=signs[offset].s[0]*bufa+signs[offset].s[1]*bufb;
	
	bufa=tmpa*(db.s[2]-(da.s[2]*costheta));
	bufb=tmpb*(da.s[2]-(db.s[2]*costheta));
	a.s[2]+=signs[offset].s[0]*bufa+signs[offset].s[1]*bufb;
}

template <typename T>
T bendP(threeVector<T> da, threeVector<T> db, T *constants)
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

//These are ultra optimized for use in CellOpt
//Non-bonded pair force function as described by Revalee, et. al., derived by Spangler.
//It is optimized for utilization in cellOpt.
//Note that there is a slow down using this function. CellOpt is capable of being 4 times faster without this function.
template <typename T>
inline void Force(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)
{
	T dx=p1.x-p2.x;//1 flop
	T dy=p1.y-p2.y;//1 flop
	T dz=p1.z-p2.z;//1 flop
	
	T dr=dx*dx+dy*dy+dz*dz;//5 flops
	
	//if(dr+=dx*dx<cutoffSquared)
	//	if(dr+=dy*dy<cutoffSquared)
	//		if(dr+=dz*dz<cutoffSquared)
	//if(dx+dy+dz<=fC[4])//triangle inequality
	//is it in range?
	#ifndef SOLVENT_FLAG
		if(dr<cutoffSquared)
	#else
		int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		//also, don't do solvent solvent interactions
		if(dr<cutoffSquared && cindex!=nTWOBODYFCONST*((SOLVENT_FLAG*nT)+SOLVENT_FLAG))
	#endif
	{
		dr=sqrt(dr);
		#ifndef SOLVENT_FLAG
			int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		#endif
		//floor() is worse than int()
		cindex+=(int(dr/(fC[cindex]))*3);
		T magnitude=fC[cindex]-dr;//1 flop
		magnitude=((fC[cindex+1]-fC[cindex+2]*magnitude)*magnitude)/dr;//4 flops
		
		dx*=magnitude;//1 flop
		dy*=magnitude;//1 flop
		dz*=magnitude;//1 flop
		
		//newton's third law
		#ifndef FORCE_COUNT
			a1.x+=dx;//1 flop
			a1.y+=dy;//1 flop
			a1.z+=dz;//1 flop
			a2.x-=dx;//1 flop
			a2.y-=dy;//1 flop
			a2.z-=dz;//1 flop
		#else
			a1.x+=1;//1 flop
			a1.y+=1;//1 flop
			a1.z+=1;//1 flop
			a2.x+=1;//1 flop
			a2.y+=1;//1 flop
			a2.z+=1;//1 flop
		#endif
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

template <typename T>
inline void ForceSelf(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, T *fC)
{
	T dx=p1.x-p2.x;//1 flop
	T dy=p1.y-p2.y;//1 flop
	T dz=p1.z-p2.z;//1 flop
	
	T dr=dx*dx+dy*dy+dz*dz;//5 flops
	
	//if(dr+=dx*dx<cutoffSquared)
	//	if(dr+=dy*dy<cutoffSquared)
	//		if(dr+=dz*dz<cutoffSquared)
	//if(dx+dy+dz<=fC[4])//triangle inequality
	//is it in range?
	#ifndef SOLVENT_FLAG
		if(dr<cutoffSquared)
	#else
		int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		//also, don't do solvent solvent interactions
		if(dr<cutoffSquared && cindex!=nTWOBODYFCONST*((SOLVENT_FLAG*nT)+SOLVENT_FLAG))
	#endif
	{
		dr=sqrt(dr);
		#ifndef SOLVENT_FLAG
			int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		#endif
		//floor() is worse than int()
		cindex+=(int(dr/(fC[cindex]))*3);
		T magnitude=fC[cindex]-dr;//1 flop
		magnitude=((fC[cindex+1]-fC[cindex+2]*magnitude)*magnitude)/dr;//4 flops
		
		dx*=magnitude;//1 flop
		dy*=magnitude;//1 flop
		dz*=magnitude;//1 flop
		
		//only apply force to self
		a1.x+=dx;//1 flop
		a1.y+=dy;//1 flop
		a1.z+=dz;//1 flop
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

///Non-bonded pair potential function as described by Revalee, et. al..
///It is optimized for use in cellOpt.
template <typename T>
inline T Potential(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC)
{
	T dx,dy,dz,dr,potential=0;
	dx=p1.x-p2.x;
	dy=p1.y-p2.y;
	dz=p1.z-p2.z;
	
	dr=dx*dx+dy*dy+dz*dz;
	//is it in range?
	#ifndef SOLVENT_FLAG
		if(dr<cutoffSquared)
	#else
		int cindex=nTWOBODYUCONST*((p1.type*nT)+p2.type);
		//also, don't do solvent solvent interactions
		if(dr<cutoffSquared && cindex!=nTWOBODYUCONST*((SOLVENT_FLAG*nT)+SOLVENT_FLAG))
	#endif
	{
		dr=sqrt(dr);
		#ifndef SOLVENT_FLAG
			int cindex=nTWOBODYUCONST*((p1.type*nT)+p2.type);
		#endif
		
		if(dr<=uC[cindex])
		{
			potential=uC[cindex]-dr;
			potential=uC[cindex+1]*potential*potential+uC[cindex+2];
		}
		else
		{
			potential=uC[cindex+3]-dr;
			potential=potential*potential*(uC[cindex+4]-potential*uC[cindex+5]);
		}
	}
	return potential;
}

///Counter non-bonded pair force. To subtract bonds if it isn't easy to do otherwise. Optimized for cellOpt.
template <typename T>
inline void CounterForce(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, threeVector<T> &a1, threeVector<T> &a2, T *fC)
{
	T dx,dy,dz,dr;

	dx=p1.x-p2.x;//1 flop
	dy=p1.y-p2.y;//1 flop
	dz=p1.z-p2.z;//1 flop
	
	dr=dx*dx+dy*dy+dz*dz;//5 flops
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		int cindex=nTWOBODYFCONST*((p1.type*nT)+p2.type);
		cindex+=(int(dr/(fC[cindex]))*3);
		T magnitude=fC[cindex]-dr;//1 flop
		magnitude=((fC[cindex+1]-fC[cindex+2]*magnitude)*magnitude)/dr;//4 flops
		
		dx*=magnitude;//1 flop
		dy*=magnitude;//1 flop
		dz*=magnitude;//1 flop
		
		a1.x-=dx;//1 flop
		a1.y-=dy;//1 flop
		a1.z-=dz;//1 flop
		a2.x+=dx;//1 flop
		a2.y+=dy;//1 flop
		a2.z+=dz;//1 flop
		
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

///Counter non-bonded pair force. To subtract bonds if it isn't easy to do otherwise. Optimized for cellOpt.
template <typename T>
inline void CounterForce(T cutoffSquared, int nT, threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *fC, int type1, int type2)
{
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		
		int cindex=nTWOBODYFCONST*((type1*nT)+type2);
		cindex+=(int(dr/(fC[cindex]))*3);
		T magnitude=fC[cindex]-dr;//1 flop
		magnitude=((fC[cindex+1]-fC[cindex+2]*magnitude)*magnitude)/dr;//4 flops
		
		d.x*=magnitude;//1 flop
		d.y*=magnitude;//1 flop
		d.z*=magnitude;//1 flop
		
		a1.x-=d.x;//1 flop
		a1.y-=d.y;//1 flop
		a1.z-=d.z;//1 flop
		a2.x+=d.x;//1 flop
		a2.y+=d.y;//1 flop
		a2.z+=d.z;//1 flop
		
	}
	//22 "useful" flops total, compare to a lennard jones 23
}

///Counter non-bonded pair potential. To subtract bonds if it isn't easy to do otherwise. Optimized for cellOpt.
template <typename T>
inline T CounterPotential(T cutoffSquared, int nT, position<T> &p1, position<T> &p2, T *uC)
{
	T dx,dy,dz,dr,potential=0;
	dx=p1.x-p2.x;
	dy=p1.y-p2.y;
	dz=p1.z-p2.z;
	
	dr=dx*dx+dy*dy+dz*dz;
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		int cindex=nTWOBODYUCONST*((p1.type*nT)+p2.type);
		
		if(dr<=uC[cindex])
		{
			potential=uC[cindex]-dr;
			potential=uC[cindex+1]*potential*potential+uC[cindex+2];
		}
		else
		{
			potential=uC[cindex+3]-dr;
			potential=potential*potential*(uC[cindex+4]-potential*uC[cindex+5]);
		}
	}
	return 0-potential;
}

///Counter non-bonded pair potential. To subtract bonds if it isn't easy to do otherwise. Optimized for cellOpt.
template <typename T>
inline T CounterPotential(T cutoffSquared, int nT, threeVector<T> d, T *uC, int type1, int type2)
{
	T potential=0;
	
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	//is it in range?
	if(dr<cutoffSquared)
	{
		dr=sqrt(dr);
		int cindex=nTWOBODYUCONST*((type1*nT)+type2);
		
		if(dr<=uC[cindex])
		{
			potential=uC[cindex]-dr;
			potential=uC[cindex+1]*potential*potential+uC[cindex+2];
		}
		else
		{
			potential=uC[cindex+3]-dr;
			potential=potential*potential*(uC[cindex+4]-potential*uC[cindex+5]);
		}
	}
	return 0-potential;
}


/** \b Maximum torque with a 3rd order polynomial dropoff
 *  This force is non-conservative. It will raise the temperature of the system.
 *  constants[0], constants[1], and constants[2] defines the x, y, and z components, respectively,
 *  of the orientation vector. It should be normalized. constants[3] and constants[4] define the
 *  torque constant and maximum torque onset angle. The orientation vector sits at the center of 
 *  mass. This force does not maintain the distance between the two particles, but you can add 
 *  a harmonic spring with the two body bonding force or really any two body force that attracts.
 * 
 *        p--->F
 *       / \
 *        |\ theta
 *        |--->z
 *        |
 *        |
 *  -F<---p
 */
template <typename T>
threeVector<T> kmaxTorqueF(threeVector<T> r, T *constants)//threeVector<T> z, T k, T thetaD)
{
	//normalized displacement vector
	T dr=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
	r.x/=dr;
	r.y/=dr;
	r.z/=dr;
	
	//sine of the angle formed between the displaced and prefered directions
	//this is negative because we are using p--->p--->p convention rather than p<---p--->p convention
	T theta=asin(-((r.x*constants[0])+(r.y*constants[1])+(r.z*constants[2])));
	//T thetaC=acos(-((r.x*z.x)+(r.y*z.y)+(r.z*z.z)));
	T thetaSqr=theta*theta;
	//T thetaAbs=abs(theta);
	
	//magnitude of the force vector
	T magnitude=0;
	if(theta<constants[4])
		magnitude=constants[3];
	else
		magnitude=constants[3]*(-2.0*theta*thetaSqr+3.0*(M_PI/2.0+constants[4])*thetaSqr-6.0*(M_PI/2.0)*constants[4]*theta+(M_PI*M_PI/4.0)*(3.0*constants[4]-M_PI/2.0))/pow(constants[4]-(M_PI/2.0),3.0);
	
	//vector in the direction of the torqued axis
	threeVector<T> crossProd;
	crossProd.x=(r.y*constants[2]-r.z*constants[1]);
	crossProd.y=-(r.x*constants[2]-r.z*constants[0]);
	crossProd.z=(r.x*constants[1]-r.y*constants[0]);
	
	//obtaining the correct length to find the vector tangent to rotation
	dr=sqrt(crossProd.x*crossProd.x+crossProd.y*crossProd.y+crossProd.z*crossProd.z);
	crossProd.x/=dr;
	crossProd.y/=dr;
	crossProd.z/=dr;
	
	//the vector tangent to the direction of rotation
	threeVector<T> f;
	f.x=(crossProd.y*r.z-crossProd.z*r.y)*magnitude;
	f.y=-(crossProd.x*r.z-crossProd.z*r.x)*magnitude;
	f.z=(crossProd.x*r.y-crossProd.y*r.x)*magnitude;
	
	return f;
}

template <typename T>
T zTorquePotential(threeVector<T> d, T z, T *constants)
{
	T cosTheta=std::abs(d.z)/magnitude(d);
	//T epsilon=a-b*tanh(c*(z-d));
	T epsilon=constants[0]-constants[1]*tanh(constants[2]*(z-constants[3]));
	T potential=epsilon*(1.0-cosTheta*cosTheta)/2.0;
	return potential;
}

template <typename T>
void zTorqueForce(threeVector<T> d, T z, threeVector<T> &a1, threeVector<T> &a2, T *constants)
{
	T mag=magnitude(d);
	T cosTheta=std::abs(d.z)/mag;
	T epsilon=constants[0]-constants[1]*tanh(constants[2]*(z-constants[3]));
	if(z>constants[3]) epsilon=0;
	T m=-epsilon*(d.z/(mag*mag*mag*mag));
	a1.x+=d.z*m*d.x;
	a1.y+=d.z*m*d.y;
	a1.z+=-m*(d.x*d.x+d.y*d.y);
	a2.x+=-d.z*m*d.x;
	a2.y+=-d.z*m*d.y;
	a2.z+=m*(d.x*d.x+d.y*d.y);
}

template <typename T>
T zPowerPotential(threeVector<T> d, T *constants)
{
	return -constants[0]*std::pow(d.z,constants[1]+1)/(constants[1]+1);
}

template <typename T>
void zPowerForce(threeVector<T> d, threeVector<T> &a, T *constants)
{
	//constants[0]=-3.0*0.000085/4.0;
	//constants[1]=2.1;
	a.z+=constants[0]*std::pow(d.z,constants[1]);
}

//Constants (same as harmonic):
//constants[0]=r_0
//constants[1]=k
//Derived from:
//r_0=R_n+1 where R_n is radius of NP
//k=2*k_R/(R_n-r_0)^2 where k_R is the maximum potential at surface of NP
//And if r>r_0, potential=0
template <typename T>
void harmonicHalfF(threeVector<T> d, threeVector<T> &a1, threeVector<T> &a2, T *constants)
{
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	
	if(dr<constants[0]*constants[0])
	{
		dr=sqrt(dr);
		T magnitude=dr-constants[0];//1 flop
		magnitude=-magnitude*constants[1]/dr;//2 flops
		
		d.x*=magnitude;
		d.y*=magnitude;
		d.z*=magnitude;
		
		a1.x+=d.x;
		a1.y+=d.y;
		a1.z+=d.z;
		
		a2.x-=d.x;
		a2.y-=d.y;
		a2.z-=d.z;
	}
}

template <typename T>
void harmonicHalfF(threeVector<T> &d, threeVector<T> &a, T *constants)
{
	T dr=d.x*d.x+d.y*d.y+d.z*d.z;
	
	if(dr<constants[0]*constants[0])
	{
		dr=sqrt(dr);
		T magnitude=dr-constants[0];//1 flop
		magnitude=-magnitude*constants[1]/dr;//2 flops
		
		d.x*=magnitude;
		d.y*=magnitude;
		d.z*=magnitude;
		
		a.x+=d.x;
		a.y+=d.y;
		a.z+=d.z;
	}
}

template <typename T>
T harmonicHalfP(threeVector<T> d, T *constants)
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
