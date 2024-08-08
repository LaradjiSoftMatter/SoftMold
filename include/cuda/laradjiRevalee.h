#include <cmath>
#include <vector>
#include "dataTypes.h"

#ifndef LARADJIREVALEE_MPD
#define LARADJIREVALEE_MPD


namespace mpd {
	
	const int nLaradjiRevaleeFC=6;
	const int nLaradjiRevaleePC=6;
	
	//int constantOffset(const int &type1, const int &type2, const int &nTypes)
	//{
	//	return type1*nTypes+type2;
	//}
	
	/** \b Conservative force between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Joel Revalee in 
	 *  the journal of Chemical Physics doi:10.1063/1.2825300
	 **/
	template <typename T>
	threeVector<T> laradjiRevaleeF(const threeVector<T> &d, const T &cutoffSquared, 
					T* constants, const int &type1, 
					const int &type2, const int &nTypes)
	{
		threeVector<T> a=0;
		//3 flops from the difference vector d
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
		//is it in range?
		if(dr<cutoffSquared)
		{
			dr=sqrt(dr);
			T *c=&constants[nLaradjiRevaleeFC*((type1*nTypes)+type2)];
			
			int offset=int(dr/(c[0]))*3;
			T magnitude=c[offset]-dr;//1 flop
			magnitude=((c[offset+1]-c[offset+2]*magnitude)*magnitude)/dr;//4 flops
			
			a=d;
			a.x*=magnitude;//1 flop
			a.y*=magnitude;//1 flop
			a.z*=magnitude;//1 flop
		}
		//16 "useful" flops total, compare to a lennard jones 17
		return a;
	}
	
	/** \b Conservative potential between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Joel Revalee in 
	 *  the journal of Chemical Physics doi:10.1063/1.2825300
	 **/
	template <typename T>
	T laradjiRevaleeP(const threeVector<T> &d, const T &cutoffSquared, 
			T* constants, const int &type1, 
			const int &type2, const int &nTypes)
	{
		T potential=0;
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		//is it in range?
		if(dr<cutoffSquared)
		{
			dr=std::sqrt(dr);
			T *c=&c[nLaradjiRevaleePC*((type1*nTypes)+type2)];
			
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
	
	/** \b Conservative force constants between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Joel Revalee in 
	 *  the journal of Chemical Physics doi:10.1063/1.2825300
	 *  Generates a flat matrix of values of the form:
	 *  a(0,0),b(0,0),c(0,0),d(0,0),e(0,0),f(0,0),
	 *  a(1,0),b(1,0),c(1,0),d(1,0),e(1,0),f(1,0),...
	 *  a(n,n),b(n,n),c(n,n),d(n,n),e(n,n),f(n,n)
	 *  Where Umax and Umin are also flat matrices with size n*n, and
	 *  n is the number of types for the interaction.
	 **/
	template <typename T>
	std::vector<T> laradjiRevaleeFC(const std::vector<T> &Umax, const std::vector<T> &Umin, const T &rc, const T &rm)
	{
		std::vector<T> output;
		#ifdef ERROR_MPD
		if(Umax.size()!=Umin.size())
			error_mpd("Umax and Umin sizes are not the same! ");
		#endif
		auto nTypes=std::sqrt(Umax.size());
		#ifdef ERROR_MPD
		if(Umax.size()%nTypes!=0)
			error_mpd("Umax and Umin matrix must be square! ");
		#endif
		for(int i=0;i<nTypes;i++)
		{
			for(int j=0;j<nTypes;j++)
			{
				int g=i+j*nTypes;
				
				output.push_back(rm);//C0
				output.push_back((2.0*(Umax[g]-Umin[g]))/(rm*rm));//C1
				output.push_back(0);//C2, part of index trick
				output.push_back(rc);//C3
				output.push_back((6.0*Umin[g])/((rc-rm)*(rc-rm)));//C4
				output.push_back((6.0*Umin[g])/((rc-rm)*(rc-rm)*(rc-rm)));//C5
			}
		}
		return output;
	}
	
	/** \b Conservative potential constants between non-bonded particles.
	 *  Model as described by Mohamed Laradji and Joel Revalee in 
	 *  the journal of Chemical Physics doi:10.1063/1.2825300
	 *  Generates a flat matrix of values of the form:
	 *  a(0,0),b(0,0),c(0,0),d(0,0),e(0,0),f(0,0),
	 *  a(1,0),b(1,0),c(1,0),d(1,0),e(1,0),f(1,0),...
	 *  a(n,n),b(n,n),c(n,n),d(n,n),e(n,n),f(n,n)
	 *  Where Umax and Umin are also flat matrices with size n*n, and
	 *  n is the number of types for the interaction.
	 **/
	template <typename T>
	std::vector<T> laradjiRevaleePC(const std::vector<T> &Umax, const std::vector<T> &Umin, const T &rc, const T &rm)
	{
		std::vector<T> output;
		#ifdef ERROR_MPD
		if(Umax.size()!=Umin.size())
			error_mpd("Umax and Umin sizes are not the same! ");
		#endif
		auto nTypes=std::sqrt(Umax.size());
		#ifdef ERROR_MPD
		if(Umax.size()%nTypes!=0)
			error_mpd("Umax and Umin matrix must be square! ");
		#endif
		for(int i=0;i<nTypes;i++)
		{
			for(int j=0;j<nTypes;j++)
			{
				int g=i+j*nTypes;
				
				output.push_back(rm);//C0
				output.push_back((Umax[g]-Umin[g])/(rm*rm));//C1
				output.push_back(Umin[g]);//C2
				output.push_back(rc);//C3
				output.push_back((3.0*Umin[g])/((rc-rm)*(rc-rm)));//C4
				output.push_back((2.0*Umin[g])/((rc-rm)*(rc-rm)*(rc-rm)));//C5
			}
		}
		return output;
	}
}

#endif
