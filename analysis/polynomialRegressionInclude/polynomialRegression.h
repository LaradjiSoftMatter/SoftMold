#ifndef POLYNOMIAL_REGRESSION_MPD
#define POLYNOMIAL_REGRESSION_MPD

#include <vector>
#include <cmath>

#ifdef VERBOSE_MPD
	#include <iostream>
#endif

#include "pivotize.h"
#include "gaussianElimination.h"

namespace mpd {
	
	#ifdef VERBOSE_MPD
	template <typename T>
	std::ostream &operator << (std::ostream &stream, std::vector<std::vector<T>> m)
	{
		for(auto n:m)
		{
			for(auto i:n)
				stream << i << '\t';
			stream << std::endl;
		}
		return stream;
	}
	#endif
	
	/** \b Polynomial regression of vector_n(y)=vector_n(m)*matrix_n_degree(x)+vector_n(b).
	 *
	 *
	**/
	template <typename T>
	std::vector<T> polynomialRegression(std::vector<T> x, std::vector<T> y, unsigned degree)
	{
		unsigned dp1=degree+1;
		unsigned dp2=degree+2;
		unsigned tp=2*degree+1;
		std::vector<T> X(tp,0);
		for(unsigned i=0;i<tp;++i)
			for(unsigned j=0;j<x.size();++j)
				X[i]+=std::pow(x[j],i);
		
		std::vector<T> out(dp1,0);
		
		std::vector<std::vector<T>> m(dp1,std::vector<T>(dp2,0));
		for(unsigned i=0;i<dp1;++i)
			for(unsigned j=0;j<dp1;++j)
				m[i][j]=X[i+j];
		
		#ifdef VERBOSE_MPD
			std::cout << "Initial Matrix:" << std::endl;
			std::cout << m << std::endl;
		#endif
		
		std::vector<T> Y(dp1,0);
		for(unsigned i=0;i<dp1;++i)
			for(unsigned j=0;j<y.size();++j)
				Y[i]+=std::pow(x[j],i)*y[j];
		
		for(unsigned i=0;i<dp1;++i)
			m[i][dp1]=Y[i];
		
		#ifdef VERBOSE_MPD
			std::cout << "Full Matrix:" << std::endl;
			std::cout << m << std::endl;
		#endif
		
		m=pivotize(m);
		
		#ifdef VERBOSE_MPD
			std::cout << "Pivotized Matrix:" << std::endl;
			std::cout << m << std::endl;
		#endif
		
		m=gaussianElimination(m);
		
		#ifdef VERBOSE_MPD
			std::cout << "Gaussian Eliminated Matrix:" << std::endl;
			std::cout << m << std::endl;
		#endif
		
		for(unsigned i=m.size()-1;i<m.size();--i)
		{
			out[i]=m[i].back();
			for(unsigned j=0;j<dp1;++j)
				if(i!=j)
					out[i]-=m[i][j]*out[j];
			out[i]/=m[i][i];
		}
		
		#ifdef VERBOSE_MPD
			std::cout << "returning..." << std::endl;
		#endif
		
		return out;
	}

}

#endif
