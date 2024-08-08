#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

typedef struct {
	double x,y,z;
	int type;
} position;

bool compare(position a, position b)
{
	return a.type>b.type;
}

int main(int argc, char **argv)
{
	if(argc!=2)
	{
		cout << "Usage: " << argv[0] << " frames.xyz\n";
		return 0;
	}

	position *p;
	fstream input;
	input.open(argv[1],ios::in);

	if(!input.is_open())
	{
		cout << "File not found!\n";
		return 0;
	}

	int n;
	char buf[256];
	input >> n >> buf;
	p=new position[n];

	while(!input.eof())
	{
		for(int i=0;i<n;i++)
			input >> p[i].type >> p[i].x >> p[i].y >> p[i].z;
		sort(&p[0],&p[n],compare);
		cout << n << endl << buf << endl;
		for(int i=0;i<n;i++)
			cout << p[i].type << '\t' << p[i].x << '\t' << p[i].y << '\t' << p[i].z << endl;
		input >> n >> buf;
	}

	input.close();
	return 0;
}
