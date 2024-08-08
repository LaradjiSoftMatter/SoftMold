#include <iostream>
#include <vector>


int main()
{
	std::vector<double> f,x;
	double ff,xx;
	while(std::cin >> xx >> ff)
	{
		f.push_back(ff);
		x.push_back(xx);
	}
	
	double mean=0;
	for(auto &y:f)
		mean+=y;
	mean/=f.size();
	
	std::vector<double> fbar(f.size(),0);
	double n=0;
	for(int i=0;i<f.size();i++)
	{
		double yi=f[i]-mean;
		for(int j=i;j<f.size();j++)
		{
			double yj=f[j]-mean;
			fbar[j-i]+=yi*yj;
		}
		n+=yi*yi;
	}
	for(int i=0;i<fbar.size();i++)
	{
		std::cout << x[i] << ' ' << fbar[i]/n << std::endl;
	}
	
	return 0;
}
