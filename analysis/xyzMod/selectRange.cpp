//Takes an xyz file and selects a range of frames to output.
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
//	for(int i=0;i<argc;i++)
//		cout << argv[i] << endl;
	if(argc!=4)
	{
		cout << "usage: " << argv[0] << " filename.xyz firstFrame lastFrame" << endl;
		return 0;
	}
	fstream inputFile;
	inputFile.open(argv[1], ios::in);
	int firstFrame=atoi(argv[2]);
	int lastFrame=atoi(argv[3]);

	if(!inputFile.is_open())
	{
		cout << argv[1] << " not found!" << endl;
		return 0;
	}
	
	for(int frameNumber=0;frameNumber<lastFrame;frameNumber++)
	{
		char line[1024];
		int n;
		inputFile.getline(line,1024);
		if(!atoi(line))
		{
			cout << "File format not recognized!" << endl;
			inputFile.close();
			return 0;
		}
		n=atoi(line);
		if(frameNumber<=lastFrame && frameNumber>=firstFrame)
			cout << line << endl;
		for(int i=0;i<=n;i++)
		{
			if(inputFile.eof())
			{
				cout << "End of file encountered!" << endl;
				inputFile.close();
				return 0;
			}			
			inputFile.getline(line,1024);
			if(frameNumber<=lastFrame && frameNumber>=firstFrame)
				cout << line << endl;
		}
	}
	inputFile.close();

	return 0;
}
