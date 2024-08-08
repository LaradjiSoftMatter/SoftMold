#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>

int main(int argc, char **argv)
{
	double threshold=0.00001;
	int offset=0;
	if(argc==2)
	{
		std::stringstream cmdArg;
		cmdArg << argv[1];
		cmdArg >> offset;
	}
	if(argc==3)
	{
		std::stringstream cmdArg;
		cmdArg << argv[1] << argv[2];
		cmdArg >> offset >> threshold;
	}
	//data is split into [frame], [lines], and [datum]
	std::vector<std::vector<std::vector<double>>> data;
	std::vector<std::vector<double>> frame;
	std::string line;
	int frameIndex=0;
	while(std::getline(std::cin,line,'\n'))
	{
		if(line.size()<2)
		{
			if(frameIndex>=offset)
				data.push_back(frame);
			frame.clear();
			frameIndex++;
		}
		else
		{
			std::stringstream splitLine;
			splitLine << line;
			std::vector<double> dataLine;
			double datum;
			while(splitLine >> datum)
				dataLine.push_back(datum);
			frame.push_back(dataLine);
		}
	}
	
	std::map<int,std::vector<std::vector<double>>> groups;
	for(auto thisFrame:data)
	{
		for(auto dataLine:thisFrame)
		{
			int index=static_cast<int>(dataLine.front()/threshold);
			auto current=groups.find(index);
			if(current!=groups.end())
				current->second.push_back(dataLine);
			else
				groups[index].push_back(dataLine);
		}	
	}
	std::map<double,std::vector<double>> averageData;
	for(auto g=groups.begin();g!=groups.end();g++)
	{
		double avgX=0,nX=0;
		std::vector<double> avgY,nY;
		for(auto dataLine:g->second)
		{
			avgX+=dataLine[0]; nX++;
			for(int i=1;i<dataLine.size();i++)
			{
				while(avgY.size()<i)
					avgY.push_back(0);
				while(nY.size()<i)
					nY.push_back(0);
				avgY[i-1]+=dataLine[i]; nY[i-1]++;
			}
		}
		double nearX=avgX/nX;
		if(g!=groups.end())
		{
			auto nextG=g; nextG++;
			if(nextG!=groups.end())
			{
				for(auto dataLine:nextG->second)
				{
					if(abs(avgX-nearX)<threshold)
					{
						avgX+=dataLine[0]; nX++;
						for(int i=1;i<dataLine.size();i++)
						{
							while(avgY.size()<i)
								avgY.push_back(0);
							while(nY.size()<i)
								nY.push_back(0);
							avgY[i-1]+=dataLine[i]; nY[i-1]++;
						}
					}
				}
			}
		}
		if(g!=groups.begin())
		{
			auto previousG=g; previousG--;
			if(previousG!=groups.end())
			{
				for(auto dataLine:previousG->second)
				{
					if(abs(avgX-nearX)<threshold)
					{
						avgX+=dataLine[0]; nX++;
						for(int i=1;i<dataLine.size();i++)
						{
							while(avgY.size()<i)
								avgY.push_back(0);
							while(nY.size()<i)
								nY.push_back(0);
							avgY[i-1]+=dataLine[i]; nY[i-1]++;
						}
					}
				}
			}
		}
		avgX/=nX;
		for(int i=0;i<avgY.size();i++)
			avgY[i]/=nY[i];
		averageData[avgX]=avgY;
	}
	
	for(auto avgD:averageData)
	{
		std::cout << avgD.first << '\t';
		for(auto &avgDatum:avgD.second)
		{
			std::cout << avgDatum;
			if(&avgDatum!=&avgD.second.back())
				std::cout << '\t';
			else
				std::cout << std::endl;
		}
	}
	return 0;
}
