#ifndef GAUSSIAN_ELIMINATION_MPD
#define GAUSSIAN_ELIMINATION_MPD

#include <vector>

namespace mpd {
	
	template <typename T>
	std::vector<std::vector<T>> gaussianElimination(std::vector<std::vector<T>> in)
	{
		unsigned n=in.size();
		for(unsigned i=0;i<n-1;++i)
		{
			for(unsigned k=i+1;k<n;++k)
			{
				T t=in[k][i]/in[i][i];
				for(unsigned j=0;j<=n;j++)
					in[k][j]=in[k][j]-t*in[i][j];
			}
		}
		return in;
	}
	
}

#endif
