//extracts single frames into several files with format ####.xyz
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		cout << "usage: " << argv[0] << " filename.xyz" << endl;
		return 0;
	}
	fstream inputFile, outputFile;
	inputFile.open(argv[1], ios::in);
	char line[1024];
	int n, lineNumber=0;
	if(!inputFile.is_open())
	{
		cout << argv[1] << " not found!" << endl;
		return 0;
	}
	
	for(int i=0;!inputFile.eof();i++)
	{
		sprintf(line,"%d.xyz",i);
		outputFile.open(line,ios::out);
		inputFile.getline(line,1024);
		if(!atoi(line))
		{
			inputFile.close();
			return 0;
		}
		n=atoi(line);
		outputFile << line << endl;
		inputFile.getline(line,1024);
		outputFile << line << endl;
		for(int i=0;i<n;i++)
		{
			if(inputFile.eof())
			{
				cout << "File format not recognized!" << endl;
				inputFile.close();
				return 0;
			}
			inputFile.getline(line,1024);
			outputFile << line << endl;
		}
		outputFile.close();
	}
	inputFile.close();

	return 0;
}
