//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char **argv)
{
	if(argc!=4)
	{
		std::cout << argv[0] << " name centerType sliceSize\n";
		return 0;
	}

	std::stringstream cmdArg;
	cmdArg << argv[2] << '\t' << argv[3];
	int centerType;
	double sliceSize;
	cmdArg >> centerType >> sliceSize;

	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	std::vector< std::vector<double> > profile;
	for(int i=0;i<System.readNTypes();i++)
		profile.push_back(std::vector<double>());

	std::vector< std::vector< position<double> > > cellList;
	threeVector<int> n;
	n.x=floor(System.readSize().x/System.readCutoff());
	n.y=floor(System.readSize().y/System.readCutoff());
	n.z=floor(System.readSize().z/System.readCutoff());
	
	threeVector<double> r;
	r.x=System.readSize().x/static_cast<double>(n.x);
	r.y=System.readSize().y/static_cast<double>(n.y);
	r.z=System.readSize().z/static_cast<double>(n.z);

	position<double> *p=System.getPositions();
	for(int i=0;i<System.readNParticles();i++)
	{
		threeVector<int> c;
		c.x=floor(p[i].x/r.x);
		c.y=floor(p[i].y/r.y);
		c.z=floor(p[i].z/r.z);
		int hash=c.x+c.y*n.x;
		while(cellList.size()<=hash)
			cellList.push_back(std::vector< position<double> >());
		cellList[hash].push_back(p[i]);
	}
	
	double minZ=0;
	double maxZ=System.readSize().z;
	for(int i=0;i<cellList.size();i++)
	{
		double avgCenterZ=0;
		int nCenter=0;
		for(int j=0;j<cellList[i].size();j++)
		{
			if(cellList[i][j].type==centerType)
			{
				avgCenterZ+=cellList[i][j].z;
				nCenter++;
			}
		}
		if(nCenter>0)
			avgCenterZ/=static_cast<double>(nCenter);
		for(int j=0;j<cellList[i].size();j++)
		{
			cellList[i][j].z-=avgCenterZ;
			if(cellList[i][j].z>maxZ)
				maxZ=cellList[i][j].z;
			if(cellList[i][j].z<minZ)
				minZ=cellList[i][j].z;
		}
	}

	for(int i=0;i<cellList.size();i++)
	{
		for(int j=0;j<cellList[i].size();j++)
		{
			int pIndex=floor((cellList[i][j].z-minZ)/sliceSize);
			while(profile[cellList[i][j].type].size()<=pIndex)
				profile[cellList[i][j].type].push_back(0);
			profile[cellList[i][j].type][pIndex]+=1.0/(n.x*n.y*r.x*r.y*sliceSize);
		}
	}

	for(int i=0;i<profile.size();i++)
	{
		for(int j=0;j<profile[i].size();j++)
			std::cout << (static_cast<double>(j)*sliceSize) << '\t' << profile[i][j] << '\n';
		std::cout << '\n';
	}

	return 0;
}
