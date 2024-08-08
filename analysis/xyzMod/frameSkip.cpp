//skips to every nSkip frame
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		cout << "usage: " << argv[0] << " inFile.xyz nSkip > outFile.xyz" << endl;
		return 0;
	}
	fstream inputFile;
	inputFile.open(argv[1], ios::in);
	int skip=atoi(argv[2]);
	char line[1024];
	int n, lineNumber=0;
	if(!inputFile.is_open())
	{
		cout << argv[1] << " not found!" << endl;
		return 0;
	}
	
	while(!inputFile.eof())
	{
		inputFile.getline(line,1024);
		if(!atoi(line))
		{
			inputFile.close();
			return 0;
		}
		n=atoi(line);
		if(lineNumber%skip==0)
			cout << line << '\n';
		inputFile.getline(line,1024);
		if(lineNumber%skip==0)
			cout << line << '\n';
		for(int i=0;i<n;i++)
		{
			if(inputFile.eof())
			{
				cout << "File format not recognized!" << endl;
				inputFile.close();
				return 0;
			}
			inputFile.getline(line,1024);
			if(lineNumber%skip==0)
				cout << line << '\n';
		}
		lineNumber++;
	}
	inputFile.close();

	return 0;
}
