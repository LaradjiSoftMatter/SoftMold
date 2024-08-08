#include <iostream>
#include <sstream>
#include <random>
#include <chrono>
#include <vector>
#include "../include/algorithms/cellAlgorithm.h"

int main(int argc, char **argv)
{
	if(argc!=10)
	{
		std::cerr << "Usage: " << argv[0] << " nParticles xSize ySize zSize cutoff range periodic nTrials seed" << std::endl;
		std::cerr << "\tNote 1=true and 0=false for periodic." << std::endl;
		return 0;
	}
	
	threeVector<double> size;
	int nParticles, nTrials, seed;
	double cutoff, range;
	bool periodic;
	
	std::stringstream cmdArg;
	for(int arg=1;arg<argc;arg++)
		cmdArg << argv[arg] << ' ';
	cmdArg >> nParticles >> size.x >> size.y >> size.z >> cutoff >> range >> periodic >> nTrials >> seed;
	
	//std::chrono::system_clock 
	std::mt19937 randNum(seed);//std::chrono::system_clock::now());
	
	typedef std::chrono::high_resolution_clock clock;
	
	clock::time_point start = clock::now();
	
	std::vector< position<double> > p;
	std::vector<bool> flag;
	//std::cout << nParticles << "\ntest\n";
	for(int i=0;i<nParticles;i++)
	{
		position<double> toPlace;
		toPlace.x=size.x*static_cast<double>(randNum())/static_cast<double>(randNum.max());
		toPlace.y=size.y*static_cast<double>(randNum())/static_cast<double>(randNum.max());
		toPlace.z=size.z*static_cast<double>(randNum())/static_cast<double>(randNum.max());
		//std::cout << 1 << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
		p.push_back(toPlace);
		flag.push_back(false);
	}
	clock::time_point end=clock::now();
	std::chrono::milliseconds timeSpan=std::chrono::duration_cast<std::chrono::milliseconds > (end-start);;
	
	std::cerr << timeSpan.count() << " milliseconds to generate " << p.size() << " targets." << std::endl;
	
	Cell<double> cellList(&p[0], p.size(), cutoff, size);
	
	start=clock::now();
	
	cellList.build();
	
	int count=0;
	double cutoffSqr=cutoff*cutoff;
	if(periodic)
	{
		for(int i=0;i<p.size();i++)
		{
			threeVector<double> minImg;
			for(int j=cellList.query(i,minImg);j!=-1;j=cellList.query(i,minImg))
			{
				if(i!=j)
				{
					threeVector<double> d;
					d.x=p[i].x-p[j].x+minImg.x;
					d.y=p[i].y-p[j].y+minImg.y;
					d.z=p[i].z-p[j].z+minImg.z;
					if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					{
						flag[j]=true;
						count++;
					}
				}
			}
		}
	}
	else
	{
		for(int i=0;i<p.size();i++)
		{
			for(int j=cellList.query(i);j!=-1;j=cellList.query(i))
			{
				if(i!=j)
				{
					threeVector<double> d;
					d.x=p[i].x-p[j].x;
					d.y=p[i].y-p[j].y;
					d.z=p[i].z-p[j].z;
					if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					{
						flag[j]=true;
						count++;
					}
				}
			}
		}
	}
	
	end=clock::now();
	
	timeSpan=std::chrono::duration_cast<std::chrono::milliseconds > (end-start);;
	
	std::cerr << timeSpan.count() << " milliseconds to count " << count << " targets." << std::endl;
	
	for(int i=0;i<flag.size();i++)
		if(!flag[i])
			std::cerr << i << std::endl;
	start=clock::now();
	
	cellList.build();
	
	count=0;
	
	if(periodic)
	{
		for(int i=0;i<p.size();i++)
		{
			threeVector<double> minImg;
			position<double> pi=p[i];
			for(int j=cellList.query(p[i],minImg);j!=-1;j=cellList.query(p[i],minImg))
			{
				if(i!=j)
				{
					threeVector<double> d;
					d.x=pi.x-p[j].x+minImg.x;
					d.y=pi.y-p[j].y+minImg.y;
					d.z=pi.z-p[j].z+minImg.z;
					if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					{
						flag[j]=true;
						count++;
					}
				}
			}
		}
	}
	else
	{
		for(int i=0;i<p.size();i++)
		{
			position<double> pi=p[i];
			for(int j=cellList.query(p[i]);j!=-1;j=cellList.query(p[i]))
			{
				if(i!=j)
				{
					threeVector<double> d;
					d.x=pi.x-p[j].x;
					d.y=pi.y-p[j].y;
					d.z=pi.z-p[j].z;
					if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					{
						flag[j]=true;
						count++;
					}
				}
			}
		}
	}
	
	end=clock::now();
	
	timeSpan=std::chrono::duration_cast<std::chrono::milliseconds > (end-start);;
	
	std::cerr << timeSpan.count() << " milliseconds to count " << count << " targets." << std::endl;
	
	return 0;
}
