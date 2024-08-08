#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		cout << "Usage: " << argv[0] << " prefix suffix start end\n";
		return 0;
	}

	stringstream input;
	input << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
	string prefix, suffix;
	int start, end;
	input >> prefix >> suffix >> start >> end;

	if(start>end)
	{
		int buf=start;
		start=end;
		end=buf;		
	}
	
	for(int i=start;i<end+1;i++)
	{
		fstream inputFile;
		stringstream filenamestream;
		filenamestream << prefix << i << suffix;
		//cout << filenamestream.str() << '\n';
		inputFile.open(filenamestream.str().c_str(), ios::in);

		if(inputFile.is_open())
		{
			double time, value;
			while(!inputFile.eof())
				inputFile >> time >> value;
			cout << time << '\t' << i << '\t' << value << '\n';
			inputFile.close();
		}
		else
		{
			//cout << filenamestream.str() <<  " not found!\n";
		}
	}

	return 0;
}
