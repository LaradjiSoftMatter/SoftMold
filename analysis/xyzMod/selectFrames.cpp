//Selects specific frames in an xyz file to output into a continuous output stream.

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
//	for(int i=0;i<argc;i++)
//		cout << argv[i] << endl;
	if(argc==1)
	{
		cout << "usage: " << argv[0] << " filename.xyz frame# frame# frame# ..." << endl;
		return 0;
	}
	fstream inputFile;
	inputFile.open(argv[1], ios::in);
	int nFrames=argc-2;
	int *selected=new int[nFrames];
	
	for(int i=0;i<nFrames;i++)
		selected[i]=atoi(argv[2+i]);

	if(!inputFile.is_open())
	{
		cout << argv[1] << " not found!" << endl;
		return 0;
	}
	
	for(int frameNumber=0,nextFrame=0;nextFrame<nFrames;frameNumber++)
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
		if(frameNumber==selected[nextFrame])
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
			if(frameNumber==selected[nextFrame])
				cout << line << endl;
		}
		if(frameNumber==selected[nextFrame])
			nextFrame++;
	}
	inputFile.close();

	return 0;
}
