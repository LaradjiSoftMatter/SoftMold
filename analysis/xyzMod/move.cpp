//Takes a single xyz file and outputs position vectors for one type.

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=8)
	{
		cout << "Usage: " << argv[0] << " file xOffset yOffset zOffset xSize ySize zSize" << endl;
		return 0;
	}
	double xO,yO,zO;
	xO=atof(argv[2]);
	yO=atof(argv[3]);
	zO=atof(argv[4]);
	double xS,yS,zS;
	xS=atof(argv[5]);
	yS=atof(argv[6]);
	zS=atof(argv[7]);
	int type;
	double x,y,z;
	fstream input;
	input.open(argv[1], ios::in);
	int nParticles;
	input >> nParticles;
	cout << nParticles << '\n';
	char name[512];
	input >> name;
	cout << name << '\n';
	while(!input.eof())
	{
		input >> type >> x >> y >> z;
		x+=xO;
		y+=yO;
		z+=zO;
		x-=(x>=xS)?xS:0;
		y-=(y>=yS)?yS:0;
		z-=(z>=zS)?zS:0;
		x+=(x<0)?xS:0;
		y+=(y<0)?yS:0;
		z+=(z<0)?zS:0;
		
		cout << type << '\t' << x << '\t' << y << '\t' << z << endl;
	}
	input.close();
	return 0;
}
