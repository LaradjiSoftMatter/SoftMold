#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
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
	vector< vector<double> > averageSet;
	vector<int> nValues;

	if(start>end)
	{
		int buf=start;
		start=end;
		end=buf;		
	}
	
	for(int i=start;i<=end;i++)
	{
		fstream inputFile;
		stringstream filenamestream;
		filenamestream << prefix << i << suffix;
		//cout << filenamestream.str() << '\n';
		inputFile.open(filenamestream.str().c_str(), ios::in);
		
		if(inputFile.is_open())
		{
			double time, value;
			for(int j=0;!inputFile.eof();j++)
			{	
				inputFile >> time >> value;
				if(!inputFile.eof())
				{
					if(averageSet.size()<=j)
					{
						vector<double> coord;
						coord.push_back(time);
						coord.push_back(value);
						averageSet.push_back(coord);
						nValues.push_back(1);
					}
					else
					{
						averageSet[j][0]+=time;
						averageSet[j][1]+=value;
						nValues[j]++;
					}
				}
			}
			//cout << time << '\t' << i << '\t' << value << '\n';
			inputFile.close();
		}
		else
		{
			//cout << filenamestream.str() <<  " not found!\n";
		}
	}

	for(int i=0;i<averageSet.size();i++)
	{
		for(int j=0;j<averageSet[i].size();j++)
		{
			averageSet[i][j]/=static_cast<double>(nValues[i]);
			cout << averageSet[i][j] << '\t';
		}
		cout << endl;
	}
	
	return 0;
}
