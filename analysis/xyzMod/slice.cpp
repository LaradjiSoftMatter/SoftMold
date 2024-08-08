#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

typedef struct {
	double x,y,z;
	int t;
} position;

struct compareIt {
	bool operator() (position i, position j)
	{
		double dotI= i.x*normal.x+i.y*normal.y+i.z*normal.z;
		double dotJ= j.x*normal.x+j.y*normal.y+j.z*normal.z;
		return dotI<dotJ;
	}
	position normal;
};

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		cout << "Usage: " << argv[0] << " name.xyz coord NSlices chain\n";
		cout << "coord can be X Y or Z\n";
		return 0;
	}

	fstream input;
	input.open(argv[1],ios::in);

	if(!input.is_open())
	{
		cout << argv[1] << " not found!\n";
		return 0;
	}
	
	int nSlices=atoi(argv[3]);
	int chain=atoi(argv[4]);
	
	position normal;
	switch(argv[2][0])
	{
		case 'X':
			normal.x=1;
			normal.y=0;
			normal.z=0;
			break;
		case 'Y':
			normal.x=0;
			normal.y=1;
			normal.z=0;
			break;
		case 'Z':
			normal.x=0;
			normal.y=0;
			normal.z=1;
			break;
		default:
			cout << argv[2][0] << " not recognized!\n";
			return 0;
			break;
	}

	compareIt surface;
	surface.normal=normal;

	position *list=NULL;
	position *molCom=NULL;
	int n=0;
	char buf[256];

	while(!input.eof())
	{
		int nBuf,start,end;
		double min=1000000000,max=0,dot;

		input >> nBuf >> buf;
		if(input.eof())
			break;
		if(nBuf%chain)
		{
			cout << "Molecules incorrect.\n";
			break;
		}
			
		if(nBuf!=n || list==NULL || molCom==NULL)
		{
			n=nBuf;
			if(list!=NULL)
				delete list;
			if(molCom!=NULL)
				delete molCom;
			list=new position[n];
			molCom=new position[n/chain];
		}
		for(int i=0,mol=0;i<n;i++)
		{
			input >> list[i].t >> list[i].x >> list[i].y >> list[i].z;
			
			if(i%chain)
			{
				molCom[mol].x+=list[i].x;
				molCom[mol].y+=list[i].y;
				molCom[mol].z+=list[i].z;
			}
			else
			{
				if(i!=0)
				{
					molCom[mol].x/=double(chain);
					molCom[mol].y/=double(chain);
					molCom[mol].z/=double(chain);
					mol++;
				}
				molCom[mol].x=list[i].x;
				molCom[mol].y=list[i].y;
				molCom[mol].z=list[i].z;
				molCom[mol].t=mol;
			}

			dot=normal.x*list[i].x+normal.y*list[i].y+normal.z*list[i].z;
			if(dot<min)
				min=dot;
			if(dot>max)
				max=dot;
		}

		sort(&molCom[0],&molCom[(n/chain)-1],surface);
		
		double width=(max-min)/(double)nSlices;
		
		for(int slice=0, current=0;slice<nSlices;slice++,current=nBuf)
		{
			fstream sliceFile;
			sprintf(buf,"slice_%d_%d.xyz",n,slice);
			
			sliceFile.open(buf,ios::out);
			
			if(!sliceFile.is_open())
			{
				cout << "Cannot open " << buf << endl;
				if(list!=NULL)
					delete list;
				if(molCom!=NULL)
					delete molCom;
				input.close();
				return 0;
			}
			
			dot=molCom[current].x*normal.x+molCom[current].y*normal.y+molCom[current].z*normal.z;
			for(nBuf=current;nBuf<n/chain && dot<min+(double)(slice+1)*width;nBuf++)
				dot=molCom[nBuf].x*normal.x+molCom[nBuf].y*normal.y+molCom[nBuf].z*normal.z;
			
			sliceFile << chain*(nBuf-current) << '\n' << buf << '\n';
			for(int i=current;i<nBuf;i++)
			{
				int index=chain*molCom[i].t;
				sliceFile << list[index].t << ' ' << list[index].x << ' ' << list[index].y << ' ' << list[index].z << '\n';
				sliceFile << list[index+1].t << ' ' << list[index+1].x << ' ' << list[index+1].y << ' ' << list[index+1].z << '\n';
				sliceFile << list[index+2].t << ' ' << list[index+2].x << ' ' << list[index+2].y << ' ' << list[index+2].z << '\n';
			}
			sliceFile.close();
		}
	}

	if(list!=NULL)
		delete list;
	if(molCom!=NULL)
		delete molCom;
	
	input.close();
	return 0;
}
