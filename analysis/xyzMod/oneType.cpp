//Takes a single xyz file and outputs position vectors for one type.

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
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
	if(!input.is_open())
	{
		std::cout << "Cannot find " << argv[1] << std::endl;
		return -1;
	}
	int n;
	std::string name;
	input >> n >> name;
	while(!input.eof())
	{
		input >> type >> x >> y >> z;
		if(type==typeOutput)
			cout << x << '\t' << y << '\t' << z << endl;
	}
	input.close();
	return 0;
}
