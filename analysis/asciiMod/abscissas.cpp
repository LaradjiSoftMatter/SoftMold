#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <map>

std::vector<std::string> splitLine(std::string line)
{
	std::vector<std::string> split;
	int comment=0;
	for(;comment<line.size();comment++)
		if(line[comment]=='#')
			break;
	while(line.size()!=comment)
		line.pop_back();
	if(line.size()>0)
	{
		std::stringstream processLine;
		processLine << line;
		std::string input;
		while(processLine >> input)
			split.push_back(input);
	}
	return split;
}

int main(int argc, char **argv)
{
	if(argc<3)
	{
		std::cerr << "Usage: " << argv[0] << " fileA fileB [options]" << std::endl;
		std::cerr << "Result: Obtains B's lines with column values matching A's column values." << std::endl;
		std::cerr << "[options]:" << std::endl;
		std::cerr << "\t--offsetA [float]" << std::endl;
		std::cerr << "\t\tOffset of the A data set." << std::endl;
		std::cerr << "\t--offsetB [float]" << std::endl;
		std::cerr << "\t\tOffset of the B data set." << std::endl;
		std::cerr << "\t--columnA [integer]" << std::endl;
		std::cerr << "\t\tColumn of the A data set." << std::endl;
		std::cerr << "\t--columnB [integer]" << std::endl;
		std::cerr << "\t\tColumn of the B data set." << std::endl;
		std::cerr << "\t--threshold [float]" << std::endl;
		std::cerr << "\t\tThreshold for binning." << std::endl;
		return 0;
	}
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	double offsetA=0,offsetB=0,threshold=0;
	int columnA=0,columnB=0;
	std::string fileNameA="",fileNameB="";
	int fCount=0;
	std::string option;
	while(cmdArg >> option)
	{
		if(option=="--offsetA")
		{
			cmdArg >> offsetA;
			if(cmdArg.bad() || cmdArg.fail())
			{
				std::cerr << "Bad --offsetA value!" << std::endl;
				return 1;
			}
		}
		else if(option=="--offsetB")
		{
			cmdArg >> offsetB;
			if(cmdArg.bad() || cmdArg.fail())
			{
				std::cerr << "Bad --offsetB value!" << std::endl;
				return 2;
			}
		}
		else if(option=="--columnA")
		{
			cmdArg >> columnA;
			if(cmdArg.bad() || cmdArg.fail())
			{
				std::cerr << "Bad --columnA value!" << std::endl;
				return 3;
			}
		}
		else if(option=="--columnB")
		{
			cmdArg >> columnB;
			if(cmdArg.bad() || cmdArg.fail())
			{
				std::cerr << "Bad --columnB value!" << std::endl;
				return 4;
			}
		}
		else if(option=="--threshold")
		{
			cmdArg >> threshold;
			if(cmdArg.bad() || cmdArg.fail())
			{
				std::cerr << "Bad --threshold value!" << std::endl;
				return 5;
			}
		}
		else
		{
			if(fCount==0)
			{
				fileNameA=option;
				if(cmdArg.bad() || cmdArg.fail())
				{
					std::cerr << "No file defined!" << std::endl;
					return 6;
				}
				fCount++;
			}
			else 
			{
				fileNameB=option;
				if(cmdArg.bad() || cmdArg.fail())
				{
					std::cerr << "No file defined!" << std::endl;
					return 7;
				}
			}
		}
	}
	std::map<int, std::vector<std::string> > A,B;
	std::map<std::string, std::vector<std::string> > Aexact,Bexact; 
	std::fstream fileA,fileB;
	fileA.open(fileNameA,std::ios::in);
	fileB.open(fileNameB,std::ios::in);
	if(fileA.is_open() && fileB.is_open())
	{
		std::string line;
		while(!fileA.eof())
		{
			std::getline(fileA,line,'\n');
			auto sLine=splitLine(line);
			if(sLine.size()>columnA && !fileA.eof() && !fileA.bad() && !fileA.fail())
			{
				if(threshold==0)
				{
					std::stringstream convert;
					convert << sLine[columnA];
					double value;
					convert >> value;
					value=value+offsetA;
					convert << value;
					convert >> sLine[columnA];
					Aexact[sLine[columnA] ]=sLine;
				}
				else
				{
					std::stringstream convert;
					convert << sLine[columnA];
					double value;
					convert >> value;
					value=value+offsetA;
					convert << value;
					convert >> sLine[columnA];
					A[int(value/threshold)]=sLine;
				}
			}
		}
		while(!fileB.eof())
		{
			std::getline(fileB,line,'\n');
			auto sLine=splitLine(line);
			if(sLine.size()>columnB && !fileB.eof() && !fileB.bad() && !fileB.fail())
			{
				if(threshold==0)
				{
					std::stringstream convert;
					convert << sLine[columnB];
					double value;
					convert >> value;
					value=value+offsetB;
					convert << value;
					convert >> sLine[columnB];
					Bexact[sLine[columnB] ]=sLine;
				}
				else
				{
					std::stringstream convert;
					convert << sLine[columnB];
					double value;
					convert >> value;
					value=value+offsetB;
					convert << value;
					convert >> sLine[columnB];
					B[int(value/threshold)]=sLine;
				}
			}
		}
		if(threshold==0)
		{
			for(auto a:Aexact)
			{
				auto b=Bexact.find(a.first);
				if(b!=Bexact.end())
				{
					for(auto &bC:b->second)
						std::cout << bC << ' ';
					std::cout << std::endl;
				}
			}
		}
		else
		{
			for(auto a:A)
			{
				std::stringstream convert;
				convert << a.second[columnA];
				double aValue;
				convert >> aValue;
				for(auto &bC:a.second)
					std::cout << bC << '_';
				std::cout << std::endl;
				for(int i=a.first-1;i<a.first+2;i++)
				{
					auto b=B.find(i);
					if(b!=B.end())
					{
						convert << b->second[columnB];
						double bValue;
						convert >> bValue;
						for(auto &bC:b->second)
							std::cout << bC << '|';
						std::cout << std::endl;
						if(abs(aValue-bValue)<threshold)
						{
							for(auto &bC:b->second)
								std::cout << bC << ' ';
							std::cout << std::endl;
						}
					}
				}
			}
		}
	}
	else
	{
		if(!fileA.is_open())
			std::cerr << "File " << fileNameA << " not found!" << std::endl;
		if(!fileB.is_open())
			std::cerr << "File " << fileNameB << " not found!" << std::endl;
	}
	
	return 0;
}
