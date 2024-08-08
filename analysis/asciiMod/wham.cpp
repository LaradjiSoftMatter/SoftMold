//program for averaging column data
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <iomanip>
#include <random>
#include <chrono>

#define PARTIAL_SUM_LIMIT 5
double pairSum(double *list, int start, int end)
{
	if(end-start<=PARTIAL_SUM_LIMIT)
	{
		double sum=0;
		for(int i=start;i<end;i++)
			sum+=list[i];
		return sum;
	}
	else
	{
		return pairSum(list, start, start+int((end-start)/2))+pairSum(list, start+int((end-start)/2), end);
	}
}

std::vector<double> normalize(std::vector<double> &histogram, double binWidth)
{
	double sum=0;
	std::vector<double> P;
	for(int i=0;i<histogram.size();i++)
		sum+=histogram[i];
	for(int i=0;i<histogram.size() && sum>0 && binWidth>0;i++)
		P.push_back(histogram[i]/(sum*binWidth));
	return P;
}

double potential(double k, double d, double x)
{
	return (k/2.0)*(d-x)*(d-x);
}

int main(int argc, char **argv)
{
	int nCol=2, col=1, nBootStrap=0,throwAway=0;
	double xStart=0, xEnd=0, xWidth=0, beta=1, tol=0.1;
	std::vector<int> fileArgN;
	std::vector<double> dBias, kBias;
		
	//std::stringstream cmdArg;
	//for(int i=1;i<argc;i++)
	//	cmdArg << argv[i] << ' ';
	
	//options
	for(int i=1;i<argc;i++)
	{
		if(argv[i][0]=='-' && argv[i][1]=='-')
		{
			bool found=false;
			if((strcmp(argv[i],"--h")==0) || (strcmp(argv[i],"--help")==0))
			{
				found=true;
				std::cerr << "Usage: " << argv[0] << " File1 dBias1 kBias1 ... Options\n";
				std::cerr << "Options:\n";
				std::cerr << "--nCol [integer]\n\tNumber of columns present in data.\n";
				std::cerr << "--col [integer]\n\tColumn to average.\n";
				std::cerr << "--tol [float]\n\tConvergence tolerance.\n";
				std::cerr << "--xStart [float]\n\tHistogram start. Required!\n";
				std::cerr << "--xEnd [float]\n\tHistogram end. Required!\n";
				std::cerr << "--xWidth [float]\n\tHistogram bin width. Required!\n";
				std::cerr << "--beta [float]\n\t1/(k_b*T*epsilon). default=1!\n";
				std::cerr << "--throwAway [int]\n\tNumber of initial points to throw away. default=0!\n";
				std::cerr << "--nBootStrap [integer] \n\tNumber of times to do bootstrap analysis. default=1!\n";
				
				return 0;
			}
			if(strcmp(argv[i],"--nCol")==0)
			{
				found=true;
				i++;
				nCol=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--throwAway")==0)
			{
				found=true;
				i++;
				throwAway=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--col")==0)
			{
				found=true;
				i++;
				col=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--tol")==0)
			{
				found=true;
				i++;
				tol=atof(argv[i]);
			}
			if(strcmp(argv[i],"--xStart")==0)
			{
				found=true;
				i++;
				xStart=atof(argv[i]);
			}
			if(strcmp(argv[i],"--xEnd")==0)
			{
				found=true;
				i++;
				xEnd=atof(argv[i]);
			}
			if(strcmp(argv[i],"--xWidth")==0)
			{
				found=true;
				i++;
				xWidth=atof(argv[i]);
			}
			if(strcmp(argv[i],"--beta")==0)
			{
				found=true;
				i++;
				beta=atof(argv[i]);
			}
			if(strcmp(argv[i],"--nBootStrap")==0)
			{
				found=true;
				i++;
				nBootStrap=atoi(argv[i]);
			}
			if(!found)
			{
				std::cerr << "Unrecognized option \"" << argv[i] << "\".\n\n";
				
				std::cerr << "Usage: " << argv[0] << " File1 dBias1 kBias1 ... Options\n";
				std::cerr << "Options:\n";
				std::cerr << "--nCol [integer]\n\tNumber of columns present in data.\n";
				std::cerr << "--col [integer]\n\tColumn to average.\n";
				std::cerr << "--tol [float]\n\tConvergence tolerance.\n";
				std::cerr << "--xStart [float]\n\tHistogram start. Required!\n";
				std::cerr << "--xEnd [float]\n\tHistogram end. Required!\n";
				std::cerr << "--xWidth [float]\n\tHistogram bin width. Required!\n";
				std::cerr << "--beta [float]\n\t1/(k_b*T*epsilon). default=1!\n";
				std::cerr << "--nBootStrap [integer] \n\tNumber of times to do bootstrap analysis. default=1!\n";
				
				return 0;
			}
		}
		else
		{
			if(fileArgN.size()==dBias.size() && fileArgN.size()==kBias.size())
			{
				fileArgN.push_back(i);
				//std::cerr << "A " << fileArgN.size() << std::endl;
			}
			else if(dBias.size()==kBias.size() && fileArgN.size()>dBias.size())
			{
				dBias.push_back(atof(argv[i]));
				//std::cerr << "B " << dBias.size() << std::endl;
			}
			else if(kBias.size()<dBias.size())
			{
				kBias.push_back(atof(argv[i]));
				//std::cerr << "C " << kBias.size() << std::endl;
			}
			else
			{
				std::cerr << "File or dBias or kBias is out of order!\n";
				//std::cerr << fileArgN.size() << ' ' << dBias.size() << ' ' << kBias.size() << std::endl;
				return 0;
			}
		}
	}
	
	if(fileArgN.size()==0)
	{
		std::cerr << "No input files!\n-help for list of options and format.\n";
		return 0;
	}
	
	if(fileArgN.size()!=dBias.size() || fileArgN.size()!=kBias.size())
	{
		std::cerr << "dBias or kBias not specified for all files!\n";
		return 0;
	}

	if(xStart>=xEnd || xWidth==0 || xWidth>xEnd-xStart)
	{
		std::cerr << "Please check xStart, xEnd, and xWidth for consistency!\n";
		std::cerr << "xStart:" << xStart << " xEnd:" << xEnd << " xWidth:" << xWidth << std::endl;
		return 0;
	}
	
	//preprocessing histograms
	int xBins=static_cast<int>((xEnd-xStart)/xWidth)+1;//xBins==1 shouldn't happen if last check worked
	std::vector< std::vector<double> > histograms;
	
	for(int file=0;file<fileArgN.size();file++)
	{
		std::vector<double> histogram(xBins, 0);
		std::fstream input;
		
		input.open(argv[fileArgN[file]],std::ios::in);
		
		if(!input.is_open())
		{
			std::cerr << argv[fileArgN[file]] << " not found!\n";
			return 0;
		}
		int nRows=0;
		do 
		{
			double value, buf;
			for(int i=0;i<=col;i++)
			{
				if(input.eof())
				{
					std::cerr << "Unexpected end of file encountered before column!\n";
					std::cerr << "column:" << i << " last value:" << buf << std::endl;
					input.close();
					return 0;
				}
				if(i==col)//this is what we want
					input >> value;
				else//throw this out
					input >> buf;
				if(input.eof())
					break;
			}
			for(int i=col+1;i<nCol && !input.eof();i++)
			{
				if(input.eof())
				{
					std::cerr << "Unexpected end of file encountered after column!\n";
					std::cerr << "column:" << i << " last value:" << buf << std::endl;
					input.close();
					return 0;
				}
				input >> buf;
			}
			
			if(!input.eof() && nRows>throwAway)
			{
				int bin=(value-xStart)/xWidth;
				if(bin<0 || bin>=xBins)
				{
					//std::cerr << "Value in " << argv[fileArgN[file]] << " is out of range!\n";
					//std::cerr << "Value:" << value << " bin:" << bin << std::endl;
					//input.close();
					//return 0;
				}
				else
				{
					histogram[bin]++;
				}
			}
			nRows++;
		} while(!input.eof());
		
		//std::cout << file << std::endl;
		//std::cin.get();
		if(input.is_open())
			input.close();
		//std::cout << "Hello!" << std::endl;
		//std::cin.get();
		
		histograms.push_back(histogram);
		
		//dump it to look at it
		//for(int i=0;i<histogram.size();i++)
		//	std::cout << static_cast<double>(i)*xWidth+xStart << '\t' << histogram[i] << std::endl;
		//std::cout << std::endl;
	}
	
	//normalize(histograms[0],xWidth);
	//dump it to look at it
	//for(int i=0;i<histograms[0].size();i++)
	//	std::cout << static_cast<double>(i)*xWidth+xStart << '\t' << histograms[0][i] << std::endl;
	
	//normalize the histograms
	//for(int i=0;i<histograms.size();i++)
	//	normalize(histograms[i], xWidth);
	//leave them unbiased to keep the weight factor	

	std::vector<double> xValues, nValues;
	
	for(int i=0;i<xBins;i++)
		xValues.push_back(static_cast<double>(i)*xWidth+xStart+xWidth/2.0);
	
	//total number of values for each histogram bin
	for(int i=0;i<histograms.size();i++)
	{
		double sum=0;
		for(int j=0;j<histograms[i].size();j++)
			sum+=histograms[i][j];
		nValues.push_back(sum);
	}
	
	//our free energy
	//std::vector<double> F(xBins,0.0), Flast(xBins,0.0);
	std::vector<double> F(histograms.size(),1.0), Flast(histograms.size(),1.0);
	
	std::vector<double> P(xBins,0.0);
	bool consistent=false;
	//iterate through wham a few thousand times
	for(int iteration=0;iteration<10000 && !consistent;iteration++)
	{
		for(int i=0;i<histograms.size();i++)
		{
			Flast[i]=F[i];
			F[i]=0;
		}
		//go through each bin
		for(int i=0;i<xBins;i++)
		{
			double upperSum=0, lowerSum=0;
			for(int j=0;j<histograms.size();j++)
			{
				double upperTrial=histograms[j][i];
				double lowerTrial=exp(beta*(Flast[j]-potential(kBias[j],dBias[j],xValues[i])));
				upperSum+=upperTrial;
				lowerSum+=static_cast<double>(nValues[j])*lowerTrial;
			}
			P[i]=upperSum/lowerSum;
		}
		P=normalize(P, xWidth);
		
		for(int i=0;i<xBins;i++)
			for(int j=0;j<histograms.size();j++)
				F[j]+=P[i]*exp(-beta*potential(kBias[j],dBias[j],xValues[i]));
		double avgError=0;
		for(int i=0;i<histograms.size();i++)
			F[i]=-log(F[i])/beta;
		//free energy perturbation method
		for(int i=histograms.size()-1;i>=0;i--)
			avgError+=fabs(F[i]-Flast[i]);
		avgError/=static_cast<double>(histograms.size());
		std::cerr << "iteration: " << iteration << ' ' << std::setprecision(15) << "average error: " << avgError << std::endl;
		if(avgError<tol)
			consistent=true;
	}
	
	P=normalize(P, xWidth);
	
	for(int i=0;i<P.size();i++)
		P[i]=-log(P[i])/beta;
	
	std::mt19937 randI(std::chrono::system_clock::now().time_since_epoch().count());
	
	//bootstrap
	std::vector<std::vector<double>> Pp;
	for(int o=0;o<nBootStrap;o++)
	{
		//shuffle the bootstrap
		std::vector<int> bootstrap;
		for(int i=0;i<histograms.size();i++)
			bootstrap.push_back(randI()%histograms.size());
		
		std::vector<double> nValuesb;
		
		//total number of values for each histogram bin
		for(int i=0;i<bootstrap.size();i++)
		{
			double sum=0;
			double l=bootstrap[i];
			for(int j=0;j<histograms[l].size();j++)
				sum+=histograms[l][j];
			nValuesb.push_back(sum);
		}
		
		//our free energy
		std::vector<double> Fb(histograms.size(),1.0), Flastb(histograms.size(),1.0);
		
		std::vector<double> Pb(xBins,0.0);
		consistent=false;
		//iterate through wham a few thousand times
		for(int iteration=0;iteration<10000 && !consistent;iteration++)
		{
			for(int i=0;i<histograms.size();i++)
			{
				Flastb[i]=F[i];
				Fb[i]=0;
			}
			//go through each bin
			for(int i=0;i<xBins;i++)
			{
				double upperSum=0, lowerSum=0;
				for(int j=0;j<histograms.size();j++)
				{
					double upperTrial=histograms[bootstrap[j]][i];
					double lowerTrial=exp(beta*(Flastb[j]-potential(kBias[j],dBias[j],xValues[i])));
					upperSum+=upperTrial;
					lowerSum+=static_cast<double>(nValuesb[j])*lowerTrial;
				}
				Pb[i]=upperSum/lowerSum;
				
			}
			Pb=normalize(Pb, xWidth);
			
			for(int i=0;i<xBins;i++)
				for(int j=0;j<histograms.size();j++)
					Fb[j]+=Pb[i]*exp(-beta*potential(kBias[j],dBias[j],xValues[i]));
			double avgError=0;
			for(int i=0;i<histograms.size();i++)
				Fb[i]=-log(Fb[i])/beta;
			//free energy perturbation method
			for(int i=histograms.size()-1;i>=0;i--)
				avgError+=fabs(Fb[i]-Flastb[i]);
			avgError/=static_cast<double>(histograms.size());
			if(avgError<tol)
				consistent=true;
		}
		Pb=normalize(Pb, xWidth);
		
		for(int i=0;i<P.size();i++)
			Pb[i]=-log(Pb[i])/beta;
		Pp.push_back(Pb);
		std::cerr << "bootstrap iteration: " << o << std::endl;
	}
	std::vector<double> Pavg(P.size(),0);
	for(int i=0;i<Pp.size();i++)
		for(int j=0;j<Pp[i].size();j++)
			Pavg[j]+=Pp[i][j];
	if(Pp.size()!=0)
		for(int i=0;i<Pavg.size();i++)
			Pavg[i]/=Pp.size();
	std::vector<double> Perror(P.size(),0);
	for(int i=0;i<Pp.size();i++)
		for(int j=0;j<Pp[i].size();j++)
			Perror[j]+=pow(Pp[i][j]-Pavg[j],2);
	if(Pp.size()!=0)
		for(int i=0;i<Perror.size();i++)
			Perror[i]/=Pp.size();
	
	for(int i=0;i<P.size();i++)
	{
		//print it out
		if(P[i]!=0)
			std::cout << xValues[i] << '\t' << P[i] << '\t' << Perror[i] << std::endl;
		else
			std::cout << xValues[i] << '\t' << 0 << '\t' << 0 << std::endl;
	}
	//std::cout << std::endl;
	return 0;
}
