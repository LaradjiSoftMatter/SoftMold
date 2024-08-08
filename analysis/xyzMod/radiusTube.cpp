#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <deque>

struct position {
	double x,y,z;
	int type;
};

#define MAX_CELL_LENGTH 1024

int hash(position input, double r)
{
	return static_cast<int>(floor(input.x/r)*MAX_CELL_LENGTH*MAX_CELL_LENGTH+floor(input.y/r)*MAX_CELL_LENGTH+floor(input.z/r));
}
	

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		std::cerr << "Usage: " << argv[0] << " file.xyz headType axis cutoff" << std::endl;
		return 0;
	}
	std::stringstream cmdArg;
	cmdArg << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
	std::string name,axis;
	int headType;
	double cutoff;
	cmdArg >> name >> headType >> axis >> cutoff;
	std::fstream inputFile(name,std::ios::in);
	if(inputFile.is_open())
	{
		int n;
		std::string comment;
		std::vector<position> p;
		int frameNumber=0;
		while(inputFile >> n)
		{
			inputFile >> comment;
			std::vector<position> test;
			position s;
			for(int i=0;i<n;i++)
			{
				inputFile >> s.type >> s.x >> s.y >> s.z;
				test.push_back(s);
				if(inputFile.eof())
					std::cerr << "eof: " << inputFile.eof() << std::endl;
				if(inputFile.fail())
					std::cerr << "fail: " << inputFile.fail() << std::endl;
				if(inputFile.bad())
					std::cerr << "bad: " << inputFile.bad() << std::endl;
			}
			if(test.size()==n)
			{
				p=test;
				frameNumber++;
			}
		}
		std::cout << "Found " << frameNumber << " frames... using last one..." << std::endl;
		std::unordered_map<int,std::vector<int> >  cells;
		
		for(int i=0;i<p.size();i++)
		{
			if(p[i].type==headType)
			{
				int h=hash(p[i],cutoff);
				cells[h].push_back(i);
			}
		}
		
		std::vector<std::vector<int> > layers;
		std::vector<int> sFlag(p.size(),-1);
		for(auto& cell:cells)
		{
			std::deque<int> stack;
			stack.push_back(cell.second[0]);
			if(sFlag[cell.second[0]]==-1)
			{
				sFlag[cell.second[0]]=layers.size();
				std::vector<int> layer;
				while(stack.size()>0)
				{
					int pI=stack.back();
					stack.pop_back();
					int currentCell=hash(p[pI],cutoff);
					int cx=(currentCell/(MAX_CELL_LENGTH*MAX_CELL_LENGTH))%MAX_CELL_LENGTH;
					int cy=(currentCell/(MAX_CELL_LENGTH))%MAX_CELL_LENGTH;
					int cz=(currentCell)%MAX_CELL_LENGTH;
					layer.push_back(pI);
					for(int i=-1;i<2;i++)
					{
						for(int j=-1;j<2;j++)
						{
							for(int k=-1;k<2;k++)
							{
								int nearCellHash=((cx+i)%MAX_CELL_LENGTH)*MAX_CELL_LENGTH*MAX_CELL_LENGTH+
									         ((cy+j)%MAX_CELL_LENGTH)*MAX_CELL_LENGTH+
									         ((cz+k)%MAX_CELL_LENGTH);
								auto nearCellIter=cells.find(nearCellHash);
								if(nearCellIter!=cells.end())
								{
									for(auto &pN:nearCellIter->second)
									{
										double dx=p[pI].x-p[pN].x;
										double dy=p[pI].y-p[pN].y;
										double dz=p[pI].z-p[pN].z;
										double dr=dx*dx+dy*dy+dz*dz;
										if(sFlag[pN]==-1 && dr<cutoff*cutoff)
										{
											stack.push_back(pN);
											sFlag[pN]=layers.size();
										}
									}
								}
							}
						}
					}
				}
				layers.push_back(layer);
			}
								
		}
		std::cout << "#radius density (1 line per layer)" << std::endl;
		for(auto& layer:layers)
		{
			double comx=0,comy=0,comz=0;
			double minx=MAX_CELL_LENGTH*cutoff,maxx=0;
			double miny=MAX_CELL_LENGTH*cutoff,maxy=0;
			double minz=MAX_CELL_LENGTH*cutoff,maxz=0;
			for(auto& i:layer)
			{
				comx+=p[i].x;
				comy+=p[i].y;
				comz+=p[i].z;
				if(p[i].x>maxx) maxx=p[i].x+cutoff/2.0;
				if(p[i].y>maxy) maxy=p[i].y+cutoff/2.0;
				if(p[i].z>maxz) maxz=p[i].z+cutoff/2.0;
				if(p[i].x<minx) minx=p[i].x-cutoff/2.0;
				if(p[i].y<miny) miny=p[i].y-cutoff/2.0;
				if(p[i].z<minz) minz=p[i].z-cutoff/2.0;
			}
			comx/=double(layer.size());
			comy/=double(layer.size());
			comz/=double(layer.size());
			double r=0,min,max;
			if(axis=="x")
			{
				min=minx;
				max=maxx;
				for(auto& i:layer)
					r+=sqrt(pow(p[i].y-comy,2.0)+pow(p[i].z-comz,2.0));
			}
			if(axis=="y")
			{
				min=miny;
				max=maxy;
				for(auto& i:layer)
					r+=sqrt(pow(p[i].x-comx,2.0)+pow(p[i].z-comz,2.0));
			}
			if(axis=="z")
			{
				min=minz;
				max=maxz;
				for(auto& i:layer)
					r+=sqrt(pow(p[i].x-comx,2.0)+pow(p[i].y-comy,2.0));
			}
			r/=double(layer.size());
			
			//std::cout << "(radius=" << r << " density=" << double(layer.size())/(2.0*M_PI*(max-min)*r) << ") ";
			std::cout << r << double(layer.size())/(2.0*M_PI*(max-min)*r) << std::endl;
		}
		std::cout << std::endl;
	}
	else
	{
		std::cerr << "File " << name << " not found!" << std::endl;
		return 1;
	}
	return 0;
	
}
