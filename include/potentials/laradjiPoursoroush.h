#ifndef LARADJIPOURSOROUSH_MPD
#define LARADJIPOURSOROUSH_MPD

#include <cmath>
#include <vector>

namespace mpd {
	
	const int nLaradjiPoursoroushFC=6;
	const int nLaradjiPoursoroushPC=6;
	
	/** \b Conservative force between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Asma Poursoroush in 
	 *  the journal of Chemical Physics doi:10.1063/1.4981008
	 **/
	template <typename T>
	threeVector<T> laradjiPoursoroushF(const threeVector<T> &d, const T &cutoffSquared, 
					const std::vector<T> &constants, const int &type)
	{
		threeVector<T> a=0;
		T *c=&constants[type*nLaradjiPoursoroushFC];
		if(d.z<=c[0])
		{
			a.z=c[3]*d.z*d.z*d.z-2*c[3]*c[0]*d.z*d.z+
				(2*c[2]+c[3]*c[0]*c[0])*d.z;
		}
		else if(d.z<c[1])
		{
			T rcz=d.z-c[1];
			a.z=c[5]*(rcz*rcz)*(2*d.z*d.z-(4*c[1]-7*c[0])
				*d.z+2*c[1]*c[1]+c[0]*(-7*c[1]+6*c[0]));
		}
		return a;
	}
	
	/** \b Conservative potential between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Asma Poursoroush in 
	 *  the journal of Chemical Physics doi:10.1063/1.4981008
	 **/
	template <typename T>
	T laradjiPoursoroushP(const threeVector<T> &d, const T &cutoffSquared, 
			const std::vector<T> &constants, const int &type)
	{
		T potential=0;
		T *c=&constants[type*nLaradjiPoursoroushPC];
		if(d.z<=c[0])
		{
			T rmz=c[0]-d.z;
			potential=c[2]*(c[0]*c[0]-d.z*d.z)+
				c[3]*rmz*rmz*rmz*((c[0]/3.0)-(1.0/4.0)*rmz)+c[4];
		}
		else if(d.z<c[1])
		{
			T rcz=c[1]-d.z;
			potential=c[5]*rcz*rcz*rcz*
				(rcz*((2.0/5.0)*rcz-(7.0/4.0)*c[0])+(2.0)*c[0]*c[0]);
		}
		return potential;
	}
	
	/** \b Conservative force constants between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Asma Poursoroush in 
	 *  the journal of Chemical Physics doi:10.1063/1.4981008
	 *  Generates a flat matrix of values of the form:
	 *  a(0),b(0),c(0),d(0),e(0),f(0),
	 *  a(1),b(1),c(1),d(1),e(1),f(1),...
	 *  a(n),b(n),c(n),d(n),e(n),f(n)
	 *  Where Umax and Umin are also flat matrices with size n, and
	 *  n is the number of types for the interaction.
	 **/
	template <typename T>
	std::vector<T> laradjiPoursoroushFC(const std::vector<T> &Umax, const std::vector<T> &Umin, const T &rc, const T &rm, const T &density)
	{
		std::vector<T> output;
		#ifdef ERROR_MPD
		if(Umax.size()!=Umin.size())
			error_mpd("Umax and Umin sizes are not the same! ");
		#endif
		auto nTypes=Umax.size();
		double rcm=rc-rm;
		for(int i=0;i<nTypes;i++)
		{	
			output.push_back(rm);//0, Z direction
			output.push_back(rc);//1
			output.push_back(M_PI*density*Umin[i]);//2
			output.push_back(2*M_PI*density*(Umax[i]-Umin[i])/(rm*rm));//3
			output.push_back(2*M_PI*density*Umin[i]*((2.0/5.0)*rcm*rcm*rcm*rcm*rcm-
				rm*(7.0/4.0)*rcm*rcm*rcm*rcm+2.0*rm*rm*rcm*rcm*rcm)/(rm*rm*rm));//4
			output.push_back(2*M_PI*density*Umin[i]/(rm*rm*rm));//5
		}
		return output;
	}
}

#endif
