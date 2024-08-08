//strips types out of xyz file

#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		cout << "Usage: " << argv[0] << " file.xyz\n";
		return 0;
	}
	int type;
	double x,y,z;
	int n;
	char name[256];

	fstream input;
	input.open(argv[1], ios::in);
	if(!input.is_open())
	{
		cout << argv[1] << " not found!\n";
		return 0;
	}
	input >> n >> name;
	for(int i=0;i<n;i++)
	{
		input >> type >> x >> y >> z;
		cout << x << '\t' << y << '\t' << z << endl;
	}
	input.close();
	return 0;
}
