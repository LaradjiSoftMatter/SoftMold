#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

typedef struct {
	float x,y,z;
} position;

/**
  * \brief Creates binary file from xyz file.
  * Binary block is set up as:
  * {unsigned int nP, unsigned int nF, int t1, ..., int t[nP], 
  *  float x[1,1], float y[1,1], float z[1,1], ..., float x[nP,nF]], float y[nP,nF], float z[nP,nF]}
  * Basically, first 4 bytes is the number of particles, second 4 bytes is the number of frames.
  * The next block of 4*nP bytes are the types. The next block of 4*nP*nF bytes are the coordinates.
 **/
int main(int argc, char **argv)
{
	if(argc!=3)
	{
		cout << "Usage: " << argv[0] << " frames.xyz frames.bin\n";
		return 0;
	}

	fstream input;
	input.open(argv[1],ios::in);
	fstream output;
	output.open(argv[2], ios::out | ios::binary);

	if(!input.is_open())
	{
		cout << "File not found!\n";
		return 0;
	}

	unsigned int nP=0;
	unsigned int nF=0;
	string buf;
	input >> nP >> buf;
	position *p=new position[nP];
	int *type=NULL;
	
	while(!input.eof())
	{
		if(type==NULL)
		{
			type=new int[nP];
			for(int i=0;i<nP;i++)
				input >> type[i] >> p[i].x >> p[i].y >> p[i].z;
			if(output.is_open())
			{
				output.write((char *)(&nP), sizeof(unsigned int));
				output.write((char *)(&nF), sizeof(unsigned int));
				output.write((char *)type, nP*sizeof(int));
			}
		}
		else
		{
			for(int i=0;i<nP;i++)
				input >> type[i] >> p[i].x >> p[i].y >> p[i].z;
		}
		output.write((char *)p,nP*sizeof(position));
		nF++;
		
		if(output.eof())
		{
			cerr << "Cannot write output file! Check storage limitations!" << endl;
			input.close();
			output.close();
			return -1;
		}
		input >> nP >> buf;
	}
	
	output.seekg(sizeof(unsigned int));
	output.write((char *)&nF, sizeof(unsigned int));
	
	if(type!=NULL)
		delete type;
	delete p;
	input.close();
	output.close();
	return 0;
}
