//Takes an xyz file and outputs only selected types.

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
	if(argc==1)
	{
		cout << "usage: " << argv[0] << " filename.xyz type# type# type# ..." << endl;
		return 0;
	}
	fstream inputFile;
	inputFile.open(argv[1], ios::in);
	int nTypes=argc-2;
	int *selected=new int[nTypes];//type index
	
	//cout << nTypes << endl;

	for(int i=0;i<nTypes;i++)
		selected[i]=atoi(argv[2+i]);

	if(!inputFile.is_open())
	{
		cout << argv[1] << " not found!" << endl;
		return 0;
	}
	
	
	char line[1024];
	int n, nT;//n is number of particles, nT is the total number of the given types
	int type;
	double x,y,z;
	
	inputFile >> n >> line;//read in header
	nT=0;//I need to be explicit, the compiler will overload nT within the scope of the for statement
	for(int i=0;i<n;i++)
	{
		if(inputFile.eof())
		{
			cout << "End of file encountered!" << endl;
			inputFile.close();
			return 0;
		}			
		inputFile >> type >> x >> y >> z;
		for(int j=0;j<nTypes;j++)
			if(type==selected[j])
				nT++;
	}
	
	inputFile.seekg(ios::beg);//rewind()

	while(!inputFile.eof())
	{
		inputFile >> n >> line;//old header
		//cout << nT << endl << line << endl;//new header
		for(int i=0;i<n;i++)
		{
			if(inputFile.eof())
			{
				cerr << "End of file encountered!" << endl;
				inputFile.close();
				return 0;
			}			
			inputFile >> type >> x >> y >> z;
			
			for(int j=0;j<nTypes;j++)
				if(type==selected[j])
					//cout << type << '\t' << x << '\t' << y << '\t' << z << endl;
					cout << x << '\t' << y << '\t' << z << endl;
		}
	}
	inputFile.close();

	return 0;
}
