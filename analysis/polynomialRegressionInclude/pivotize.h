#ifndef PIVOTIZE_MPD
#define PIVOTIZE_MPD

#include <vector>
#include <cmath>

namespace mpd {
	
	template <typename T>
	std::vector<std::vector<T>> pivotize(std::vector<std::vector<T>> in)
	{
		unsigned n=in.size();
		//check that it is an nxn matrix...
		//if...
		
		std::vector<std::vector<T>> out=in;
		for(unsigned i=0;i<n;++i)
		{
			for(unsigned k=i+1;k<n;++k)
			{
				if(std::abs(in[i][i])<std::abs(in[k][i]))
				{
					for(int j=0;j<n+1;++j)
					{
						T temp=out[i][j];
						out[i][j]=out[k][j];
						out[k][j]=temp;
					}
				}
			}
		}
		return out;
	}
	
}

#endif
