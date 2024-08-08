#ifndef LARADJISPANGLER_MPD
#define LARADJISPANGLER_MPD

#include <cmath>
#include <vector>

namespace mpd {
	
	const int nLaradjiSpanglerFC=22;
	const int nLaradjiSpanglerPC=22;
	
	/** \b Conservative force between a non-bonded nanoparticle and normal particle.
	 *  Model as described by Mohamed Laradji and Eric Spangler in 
	 *  the journal of Chemical Physics doi:10.1063/1.5138897
	 **/
	template <typename T>
	threeVector<T> laradjiSpanglerF(const threeVector<T> &d, const T &cutoffSquared, 
					const std::vector<T> &constants, const int &type1, 
					const int &type2, const int &nTypes)
	{
		threeVector<T> a=0;
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
		T *c=&constants[nLaradjiSpanglerFC*((type1*nTypes)+type2)];
		T magnitude=0;
		T rminD=c[4]+c[5];
		T rminD2=rminD*rminD;
		//is it in range?
		if(rminD2<=dr && dr<cutoffSquared)
		{
			dr=std::sqrt(dr);
			T E=c[0]-dr;//1 flops
			T E2=E*E;//1 flop
			T E3=E*E2;//1 flop
			magnitude=2.0*(0.4*E3+c[1]*E2+c[2]*E)/(dr);
			magnitude+=4.0*E2+8.0*c[1]*E+6.0*c[2];
			magnitude*=(E2*c[3])/(dr*dr);//5 flops
		}
		else if(dr<rminD2)
		{
			dr=std::sqrt(dr);
			T B=rminD-dr;
			T B2=B*B;
			T B3=B2*B;
			T B4=B2*B2;
			magnitude=(B4*c[6]+B3*c[7]+B2*c[8]+B*c[9]+c[10])/(dr);
			magnitude+=(4.0*B3*c[6]+3.0*B2*c[7]+2.0*B*c[8]+c[9]);
			magnitude/=dr*dr;
		}
		
		a=d;
		a.x*=magnitude;//1 flop
		a.y*=magnitude;//1 flop
		a.z*=magnitude;//1 flop
		return a;
	}
	
	/** \b Conservative potential between a non-bonded nanoparticle and normal particle.
	 *  Model as described by Mohamed Laradji and Eric Spangler in 
	 *  the journal of Chemical Physics doi:10.1063/1.5138897
	 **/
	template <typename T>
	T laradjiSpanglerP(const threeVector<T> &d, const T &cutoffSquared, 
			const std::vector<T> &constants, const int &type1, 
			const int &type2, const int &nTypes)
	{
		T potential=0;
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		const T *c=&constants[nLaradjiSpanglerPC*((type1*nTypes)+type2)];
		T rminD=c[4]+c[5];
		T rminD2=rminD*rminD;
		if(rminD2<=dr && dr<cutoffSquared)
		{
			dr=std::sqrt(dr);
			T E=c[0]-dr;//1 flops
			T E2=E*E;//1 flop
			T E3=E*E2;//1 flop
			potential=2.0*E3*(0.4*E2+c[1]*E+c[2])*c[3]/(dr);
		}
		else if(dr<rminD2)
		{
			dr=std::sqrt(dr);
			T B=rminD-dr;
			T B2=B*B;
			T B3=B2*B;
			T B4=B2*B2;
			potential=(B4*c[6]+B3*c[7]+B2*c[8]+B*c[9]+c[10])/(dr);
		}
		return potential;
	}
	
	/** \b Conservative force between non-bonded nanoparticles.
	 *  Model as described by Mohamed Laradji and Eric Spangler in 
	 *  the journal of Chemical Physics doi:10.1063/1.5138897
	 **/
	template <typename T>
	threeVector<T> laradjiSpanglerBBF(const threeVector<T> &d, const T &cutoffSquared, 
					const std::vector<T> &constants, const int &type1, 
					const int &type2, const int &nTypes)
	{
		threeVector<T> a=0;
		const T *c=&constants[nLaradjiSpanglerFC*((type1*nTypes)+type2)+11];//11=beadbeadOffset
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;//5 flops
		T magnitude=0;
		T rminD=c[0];
		T rminD2=rminD*rminD;
		if(rminD2<=dr && dr<cutoffSquared)
		{
			dr=std::sqrt(dr);
			T E=c[1]-dr;//1 flops
			T E2=E*E;//1 flop
			T E3=E*E2;//1 flop
			magnitude=E3*(6.0*c[2]*E2+5.0*c[3]*E+4.0*c[4]+
			(c[2]*E3+c[3]*E2+c[4]*E)/dr)/(dr*dr);
		}
		else if(dr<rminD2)
		{
			dr=std::sqrt(dr);
			T B=rminD-dr;
			T B2=B*B;
			T B3=B2*B;
			T B4=B2*B2;
			T B5=B2*B3;
			magnitude=5.0*B4*c[5]+4.0*B3*c[6]+3.0*B2*c[7]+
				   2.0*B*c[8]+c[9];
			magnitude+=(B5*c[5]+B4*c[6]+B3*c[7]+
				   B2*c[8]+B*c[9]+c[10])/(dr);
			magnitude/=dr*dr;
		}
		
		a=d;
		a.x*=magnitude;//1 flop
		a.y*=magnitude;//1 flop
		a.z*=magnitude;//1 flop
		return a;
	}
	
	/** \b Conservative force between non-bonded nanoparticles.
	 *  Model as described by Mohamed Laradji and Eric Spangler in 
	 *  the journal of Chemical Physics doi:10.1063/1.5138897
	 **/
	template <typename T>
	T laradjiSpanglerBBP(const threeVector<T> &d, const T &cutoffSquared, 
			const std::vector<T> &constants, const int &type1, 
			const int &type2, const int &nTypes)
	{
		T potential=0;
		const T *c=&constants[nLaradjiSpanglerFC*((type1*nTypes)+type2)+11];//11=beadbeadOffset
		T dr=d.x*d.x+d.y*d.y+d.z*d.z;
		T rminD=c[0];
		T rminD2=rminD*rminD;
		if(rminD2<=dr && dr<cutoffSquared)
		{
			dr=std::sqrt(dr);
			T E=c[1]-dr;//1 flops
			T E2=E*E;//1 flop
			T E3=E*E2;//1 flop
			T E4=E2*E2;
			potential=E4*(c[2]*E2+c[3]*E+c[4])/dr;
		}
		else if(dr<rminD2)
		{
			dr=std::sqrt(dr);
			T B=rminD-dr;
			T B2=B*B;
			T B3=B2*B;
			T B4=B2*B2;
			T B5=B2*B3;
			potential=(B5*c[5]+B4*c[6]+B3*c[7]+
				   B2*c[8]+B*c[9]+c[10])/(dr);
		}
		return potential;
	}
	
	/** \b Conservative force constants between non-bonded nanoparticles and normal particles.
	 *  Model as described by Mohamed Laradji and Eric Spangler in 
	 *  the journal of Chemical Physics doi:10.1063/1.5138897
	 *  Generates a flat matrix of values of the form:
	 *  a(0,0),b(0,0),c(0,0),d(0,0),e(0,0),f(0,0),
	 *  a(1,0),b(1,0),c(1,0),d(1,0),e(1,0),f(1,0),...
	 *  a(n,n),b(n,n),c(n,n),d(n,n),e(n,n),f(n,n)
	 *  Where Umax and Umin are also flat matrices with size n*n, and
	 *  n is the number of types for the interaction.
	 *  This one is slightly more flexible by allow for different radii and
	 *  explicit rows.
	 **/
	template <typename T>
	std::vector<T> laradjiSpanglerFCrow(T rc, T rm, T rad, T density, T Umax, T Umin)
	{
		std::vector<T> output;
		//bead-particle constants
		output.push_back(rc+rad);//0
		output.push_back(-7.0/4.0*rm);//1
		output.push_back(2.0*rm*rm);//2
		output.push_back(Umin*M_PI*rad*density/(rm*rm*rm));//3
		output.push_back(rad);//4
		output.push_back(rm);//5
		T D=density*M_PI*rad;
		T A=Umax-Umin;
		
		output.push_back(-D*A/(2.0*rm*rm));//6,0@B^4
		output.push_back(2.0*D*A/(3.0*rm));//7,1,@B^3
		output.push_back(-D*Umin);//8,2,@B^2
		output.push_back(2.0*D*Umin*rm);//9,3,@B^1
		output.push_back(D*1.3*Umin*rm*rm);//10,4
		
		//bead-bead constants
		D=pow(2.0*M_PI*density*rad,2.0)/(rm*rm*rm);
		
		output.push_back(2.0*rad+rm);//0
		output.push_back(2.0*rad+rc);//1
		output.push_back(Umin*D*2.0/30.0);//2, x^6
		output.push_back(-Umin*D*7.0*rm/20.0);//3, x^5
		output.push_back(Umin*D*rm*rm/2.0);//4, x^4
		
		output.push_back(-D*A*rm/20.0);//5,0@B^5
		output.push_back(D*A*rm*rm/12.0);//6,1,@B^4
		output.push_back(-D*Umin*pow(rm,3.0)/6.0);//7,2,@B^3
		output.push_back(D*Umin*pow(rm,4.0)/2.0);//8,3,@B^2
		output.push_back(D*1.3*Umin*pow(rm,5.0)/2.0);//9,@B
		output.push_back(D*13.0*Umin*pow(rm,6.0)/60.0);//10
		return output;
	}
	
	/** \brief Sets up laradji-spangler force constants for many different radii
	 *  I don't recommend usig this one directly as we currently don't have a 
	 *  working bead-bead interaction for different sizes
	 **/
	template <typename T>
	std::vector<T> laradjiSpanglerFC(const std::vector<T> &Umax, const std::vector<T> &Umin, const T &rc, 
					const T &rm, const std::vector<T> &nanoRadius, const T &density)
	{
		std::vector<T> output;
		#ifdef ERROR_MPD
		if(Umax.size()!=Umin.size())
			error_mpd("Umax and Umin sizes are not the same! ");
		if(Umax.size()!=nanoRadius.size())
			error_mpd("Umax and Umin sizes are not the same! ");
		if(rc/2.0<rm-0.0001 || rc/2.0>rm+0.0001)
			error_mpd("Cutoff (rc) must be at least twice rm! ");
		#endif
		auto nTypes=std::sqrt(Umax.size());
		#ifdef ERROR_MPD
		if(Umax.size()!=nTypes*nTypes)
			error_mpd("Umax and Umin matrix must be square! ");
		#endif
		for(int i=0;i<nTypes;i++)
		{
			for(int j=0;j<nTypes;j++)
			{
				int g=i+j*nTypes;
				auto row=laradjiSpanglerFCrow(rc, rm, nanoRadius[g], density, Umax[g], Umin[g]);
				output.insert(output.end(),row.begin(),row.end());
			}
		}
		return output;
	}
	
	/** \brief Sets up laradji-spangler force constants for 1 nanoparticle radius.
	 *  I do recommend usig this one directly as we currently do have a 
	 *  working bead-bead interaction for same nanoparticle sizes.
	 **/
	template <typename T>
	std::vector<T> laradjiSpanglerFC(const std::vector<T> &Umax, const std::vector<T> &Umin, const T &rc, 
					const T &rm, const T &nanoRadius, const T &density)
	{
		std::vector<T> nRadius;
		auto nTypes=std::sqrt(Umax.size());
		for(int i=0;i<nTypes*nTypes;i++) nRadius.emplace_back(nanoRadius);
		return laradjiSpanglerFC(Umax, Umin, rc, rm, nRadius, density);
	}
}

#endif
