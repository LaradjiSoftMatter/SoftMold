//converts a target file into an xyz file

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		cout << "Usage: " << argv[0] << " file scale" << endl;
		return 0;
	}
	double scale=atof(argv[2]);
	fstream input;
	input.open(argv[1], ios::in);
	while(!input.eof())
	{
		double index, x,y,z, iX, iY,iZ;
		input >> index >> x >> y >> z >> iX >> iY >> iZ;
		cout << "1\t" << x*scale << '\t' << y*scale << '\t' << z*scale << endl;
	}
	input.close();
	return 0;
}
