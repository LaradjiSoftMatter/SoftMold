//strips types out of xyz file

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
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
		bool open[6];//{-x,x,-y,y,-z,z}
		open[0]=true;
		open[1]=true;
		open[2]=true;
		open[3]=true;
		open[4]=true;
		open[5]=true;
		double rSqr=radius*radius;
		double rrSqr=radius+radius;
		rrSqr*=rrSqr;

		

		for(int j=0;j<p.size();j++)
		{
			if(i!=j)
			{
				threeVector<double> d;
				d.x=p[j].x-p[i].x;
				d.y=p[j].y-p[i].y;
				d.z=p[j].z-p[i].z;
				double dr=d.x*d.x+d.y*d.y+d.z*d.z;
				
				if(rSqr<dr && dr<rrSqr)
				{
					if(radius<-d.x && -d.x<2.0*radius) open[0]=false;
					if(radius<d.x && d.x<2.0*radius) open[1]=false;
					if(radius<-d.y && -d.y<2.0*radius) open[2]=false;
					if(radius<d.y && d.y<2.0*radius) open[3]=false;
					if(radius<-d.z && -d.z<2.0*radius) open[4]=false;
					if(radius<d.z && d.z<2.0*radius) open[5]=false;
					//if(i==3000)
					//	cout << 1 << '\t' << p[j].x << '\t' << p[j].y << '\t' << p[j].z << endl;
				}
			}
		}
		bool oneOpen=false;
		for(int j=0;j<6;j++)
			if(open[j])
				oneOpen=true;
		if(oneOpen)
			cout << p[i].x << '\t' << p[i].y << '\t' << p[i].z << endl;
		if(i%int(p.size()/100)==0)
			cerr << int(100*i/p.size()) << endl;
	}
	
	return 0;
}
