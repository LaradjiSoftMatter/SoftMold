//strips types out of xyz file

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
using namespace std;

template <typename T>
struct threeVector {
	T x,y,z;
};

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		cout << "Usage: " << argv[0] << " file.xyz radius\n";
		return 0;
	}
	double radius=atof(argv[2]);
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
	vector<threeVector<double> > p;
	for(int i=0;i<n;i++)
	{
		int type;
		threeVector<double> pIn;
		input >> type >> pIn.x >> pIn.y >> pIn.z;
		//if(type==4)
			p.push_back(pIn);
	}
	input.close();
	
	for(int i=0;i<p.size();i++)
	{
		double rSqr=1.5*radius*radius;
		
		int nS=100;
		vector<threeVector<double> > sP;
		double s=3.6/sqrt(double(nS));
		double length=0;
		double dz=2.0/double(nS);
		double z=1.0-dz/2.0;
		for(int k=0;k<nS;k++)
		{
			threeVector<double> toPlace;
			double r=sqrt(1.0-z*z);
			toPlace.x=p[i].x+r*cos(length)*radius;
			toPlace.y=p[i].y+r*sin(length)*radius;
			toPlace.z=p[i].z+z*radius;
			sP.push_back(toPlace);
			z=z-dz;
			length=length+s/r;
		}
		
		for(int j=0;j<p.size();j++)
		{
			if(i!=j)
			{
				threeVector<double> d;
				d.x=p[j].x-p[i].x;
				d.y=p[j].y-p[i].y;
				d.z=p[j].z-p[i].z;
				
				if(d.x*d.x+d.y*d.y+d.z*d.z<rSqr)
				{
					for(int k=0;k<sP.size();k++)
					{
						d.x=p[j].x-sP[k].x;
						d.y=p[j].y-sP[k].y;
						d.z=p[j].z-sP[k].z;
						if(d.x*d.x+d.y*d.y+d.z*d.z<rSqr)
						{
							sP.erase(sP.begin()+k);
							k--;
						}
					}
				}
			}
		}
		for(int k=0;k<sP.size();k++)
			cout << sP[k].x << '\t' << sP[k].y << '\t' << sP[k].z << endl;
		if(i%int(p.size()/100)==0)
			cerr << int(100*i/p.size()) << endl;
	}
	
	return 0;
}
