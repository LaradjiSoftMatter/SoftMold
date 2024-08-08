
#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#include "cell.h"

int main()
{
	double cutoffRad=7;
	threeVector<double> size;
	size.x=29.7989898732233;
	size.y=29.7989898732233;
	size.z=280;
	std::mt19937 randG(948373);
	std::uniform_real_distribution<double> randD(-1.0,1.0);
	std::vector<position<double>> p;
	for(int i=0;i<10000;i++)
	{
		position<double> in;
		in.x=(randD(randG)+1.0)*size.x/2.0;
		in.y=(randD(randG)+1.0)*size.y/2.0;
		in.z=(randD(randG)+1.0)*size.z/2.0;
		p.push_back(in);
	}
	
	for(int step=0;step<100;step++)
	{
		threeVector<int> nCells;
		nCells.x=size.x/cutoffRad;
		nCells.y=size.y/cutoffRad;
		nCells.z=size.z/cutoffRad;
		threeVector<double> cSize;
		cSize.x=size.x/nCells.x;
		cSize.y=size.y/nCells.y;
		cSize.z=size.z/nCells.z;
		//std::cout << nCells.x << '\t' << nCells.y << '\t' << nCells.z << std::endl;
		auto cMap=createCMap(&(p.front()),&(p.back()),nCells,cSize);
		//int cCount=0;
		//for(auto &iM:iMap)
		//	for(auto &k:iM.second)
		//		cCount++;
		//std::cout << iMap.size() << '\t' << cCount << std::endl;
		
		
		double cutoffRadSqr=cutoffRad*cutoffRad;
		
		double avgCount=0;
		
		#pragma omp parallel for reduction(+:avgCount)
		for(int j=0;j<p.size();j++)
		{
			int first=j;
			auto iCell=getCell(p[first],cSize);
			int cHash=hash(iCell,nCells);
			auto neigh=neighIndices(cHash,nCells);
			double count=0;
			for(auto nHash:neigh)
			for(auto &second:cMap[nHash])
			{
				//int second=k;
				//if(first!=second)
				{
					threeVector<double> d;
					d.x=p[first].x-second.x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[first].y-second.y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[first].z-second.z;
					if(d.x>size.z/2.0) d.z-=size.z;
					if(d.x<-size.z/2.0) d.z+=size.z;
					
					double distSqr=d.x*d.x+d.y*d.y+d.z*d.z;
					if(distSqr<cutoffRadSqr && distSqr!=0)
					{
						count++;
					}
				}
			}
			//std::cout << std::endl;
			avgCount+=count;
		}
		
		//std::cout << p.size() << "\nTest\n";
		for(auto &pE:p)
		{
			//std::cout << "1 " << pE.x << ' ' << pE.y << ' ' << pE.z << '\n';
			pE.x+=randD(randG);
			pE.y+=randD(randG);
			pE.z+=randD(randG);
			while(pE.x>=size.x)pE.x-=size.x;
			while(pE.y>=size.y)pE.y-=size.y;
			while(pE.z>=size.z)pE.z-=size.z;
			while(pE.x<0.0)pE.x+=size.x;
			while(pE.y<0.0)pE.y+=size.y;
			while(pE.z<0.0)pE.z+=size.z;
		}
		//std::cout << avgCount << std::endl;
	}
	
	
	



	return 0;
}
