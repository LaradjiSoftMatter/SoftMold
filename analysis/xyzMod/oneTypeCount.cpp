//Counts a specific type.

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		cout << "Usage: " << argv[0] << " file type_to_output" << endl;
		return 0;
	}
	int type;
	double x,y,z;
	int typeOutput=atoi(argv[2]);
	fstream input;
	input.open(argv[1], ios::in);
	int count=0;
	while(!input.eof())
	{
		input >> type >> x >> y >> z;
		if(type==typeOutput)
			count++;
	}
	cout << count << endl;
	input.close();
	return 0;
}
