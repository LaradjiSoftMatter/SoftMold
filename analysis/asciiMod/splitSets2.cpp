#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

void printOptions(char *command)
{
	std::cerr << "Usage: " << command << " [options]" << std::endl;
	std::cerr << "[options]:" << std::endl;
	std::cerr << "\t--name\n\t\tName of input file. Output files are suffixed with threshold bin. Required!" << std::endl;
	std::cerr << "\t--threshold\n\t\tBin size. Required!" << std::endl;
	std::cerr << "\t--offset\n\t\tOffset to start using data sets. Optional." << std::endl;
	std::cerr << "\t--column\n\t\tColumn to use in data sets. Default=1" << std::endl;
}

int main(int argc, char **argv)
{
	double threshold=-1;
	int offset=0, column=1;
	std::string name;
	
	if(argc<3)
	{
		printOptions(argv[0]);
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	std::string option;
	while(cmdArg >> option)
	{
		if(option=="--name")
		{
			cmdArg >> name;
		}
		if(option=="--threshold")
		{
			cmdArg >> threshold;
		}
		if(option=="--offset")
		{
			cmdArg >> offset;
		}
		if(option=="--column")
		{
			cmdArg >> column;
		}
	}
	
	if(name.size()==0)
	{
		std::cerr << "--name is required!" << std::endl;
		return 1;
	}
	
	if(threshold<0)
	{
		std::cerr << "--threshold is required and positive!" << std::endl;
		return 2;
	}
	
	//data is split into [frame], [lines], and [datum]
	std::map<int,std::vector<std::vector<double>>> data;
	std::string line;
	int frameIndex=0;
	std::fstream inputFile(name,std::ios::in);
	if(!inputFile.is_open())
	{
		std::cerr << "File " << name << " not found!" << std::endl;
		return 3;
	}
	
	while(std::getline(inputFile,line,'\n'))
	{
		std::stringstream splitLine;
		splitLine << line;
		std::vector<double> dataLine;
		double datum;
		while(splitLine >> datum)
			dataLine.push_back(datum);
		if(dataLine.size()>column)
		{
			int bin=(int(dataLine[column])+offset)/threshold;
			data[bin].push_back(dataLine);
		}
	}
	inputFile.close();
	
	for(auto thisFrame:data)
	{
		std::stringstream outputFileName;
		outputFileName << name << '.' << thisFrame.first;
		std::fstream outputFile(outputFileName.str(),std::ios::out|std::ios::app);
		for(auto dataLine:thisFrame.second)
		{
			for(auto datum:dataLine)
			//	std::cout << datum << '\t';
			//std::cout << std::endl;
				outputFile << datum << '\t';
			outputFile << std::endl;
		}
		//std::cout << std::endl;
		outputFile << std::endl;
		outputFile.close();
	}
	return 0;
}
