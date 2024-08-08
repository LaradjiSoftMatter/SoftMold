#include <iostream>
#include <sstream>
#include <random>
#include <chrono>
#include <vector>
#include "../include/algorithms/cellAlgorithm.h"

threeVector<float> halfHarmonic(threeVector<float> d, float dist, float springForce)
{
        //this is a basic harmonic force, taken from -grad(k*(r-r0)^2)
        // where grad is the gradient operator, k is the spring constant,
        // r is the distance between particles, and r0 is the equilibrium distance
        //it is assumed that the spring constant takes the 2 term implicitly

        //threeVector<float> d=a-b;
        //float dist=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	
        float mag=-springConstant*(dist-equilibrium)/dist;

        threeVector<float> force;
        force.x=mag*d.x;
        force.y=mag*d.y;
        force.z=mag*d.z;
        return force;
}

int main(int argc, char **argv)
{
	if(argc!=11)
	{
		std::cerr << "Usage: " << argv[0] << " nParticles xSize ySize zSize radius periodic nSteps seed mass springForce" << std::endl;
		std::cerr << "\tNote 1=true and 0=false for periodic." << std::endl;
		return 0;
	}
	
	threeVector<float> size;
	int nParticles, nSteps, seed;
	float radius, mass, springForce;
	bool periodic;
	
	std::stringstream cmdArg;
	for(int arg=1;arg<argc;arg++)
		cmdArg << argv[arg] << ' ';
	cmdArg >> nParticles >> size.x >> size.y >> size.z >> radius >> range >> periodic >> nSteps >> seed >> mass >> springForce;
	
	//std::chrono::system_clock 
	std::mt19937 randNum(seed);//std::chrono::system_clock::now());
	
	typedef std::chrono::high_resolution_clock clock;
	
	clock::time_point start = clock::now();
	
	std::vector< position<float> > p;
	std::vector< threeVector<float> > v,f;
	//std::cout << nParticles << "\ntest\n";
	threeVector<int> n;
	n.x=size.x/(2.0*radius);
	n.y=size.y/(2.0*radius);
	n.z=size.z/(2.0*radius);
	int count=0;
	for(float x=radius;x<size.x;x+=2.0*radius)
	{
		for(float y=radius;y<size.y;y+=2.0*radius)
		{
			for(float z=radius;z<size.z;z+=2.0*radius)
			{
				if(count<nParticles)
				{
					position<float> toPlace;
					threeVector<float> toVel;
					threeVector<float> toFor=0;
					toPlace.x=x;
					toPlace.y=y;
					toPlace.z=z;
					toVel.x=static_cast<float>(randNum())/static_cast<float>(randNum.max());
					toVel.y=static_cast<float>(randNum())/static_cast<float>(randNum.max());
					toVel.z=static_cast<float>(randNum())/static_cast<float>(randNum.max());
					//std::cout << 1 << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
					p.push_back(toPlace);
					v.push_back(toVel);
					f.push_back(toFor);
					count++;
				}
			}
		}
	}
	clock::time_point end=clock::now();
	std::chrono::milliseconds timeSpan=std::chrono::duration_cast<std::chrono::milliseconds > (end-start);;
	
	std::cerr << timeSpan.count() << " milliseconds to generate " << p.size() << " targets." << std::endl;
	
	Cell<float> cellList(&p[0], p.size(), radius, size);
	
	start=clock::now();
	for(int step=0;step<nSteps;step++)
	{
		cellList.build();
		
		for(int i=0;i<a.size();i++)
			a[i]=0;
		
		float radiusSqr=radius*radius;
		if(periodic)
		{
			for(int i=0;i<p.size();i++)
			{
				threeVector<float> minImg;
				for(int j=cellList.queryHalf(i,minImg);j!=-1;j=cellList.queryHalf(i,minImg))
				{
					if(i!=j)
					{
						threeVector<float> d;
						d.x=p[i].x-p[j].x+minImg.x;
						d.y=p[i].y-p[j].y+minImg.y;
						d.z=p[i].z-p[j].z+minImg.z;
						float dr=d.x*d.x+d.y*d.y+d.z*d.z;
						if(dr<radiusSqr)
						{
							dr=sqrt(dr);
							threeVector<float> buf=springForce(d, dr, springConstant);
							f[i]+=buf;
							f[j]-=buf;
						}
					}
				}
			}
		}
		else
		{
			for(int i=0;i<p.size();i++)
			{
				for(int j=cellList.queryHalf(i);j!=-1;j=cellList.queryHalf(i))
				{
					if(i!=j)
					{
						threeVector<float> d;
						d.x=p[i].x-p[j].x;
						d.y=p[i].y-p[j].y;
						d.z=p[i].z-p[j].z;
						if(d.x*d.x+d.y*d.y+d.z*d.z<radiusSqr)
						{
							dr=sqrt(dr);
							threeVector<float> buf=springForce(d, dr, springConstant);
							f[i]+=buf;
							f[j]-=buf;
						}
					}
				}
			}
		}
	}
	end=clock::now();
	
	std::chrono::seconds executeTime=std::chrono::duration_cast<std::chrono::seconds > (end-start);;
	
	std::cerr << executeTime.count() << " seconds to execute " << nSteps << " steps." << std::endl;
	std::cerr << "Approximately " << static_cast<float>(nSteps)/static_cast<float>(executeTime.count()) << 
		     " steps per second." << std::endl;
	
	return 0;
}
