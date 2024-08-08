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
	if(argc<2)
	{
		std::cout << argv[0] << " name type1 type2 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	std::vector<int> types;
	
	for(int i=2;i<argc;i++)
	{
		int type;
		cmdArg >> type;
		types.push_back(type);
	}
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	std::vector< std::vector< position<double> > > cellList;
	threeVector<int> n;
	n.x=1024;//floor(System.readSize().x/System.readCutoff());
	n.y=1024;//floor(System.readSize().y/System.readCutoff());
	
	threeVector<double> r;
	r.x=System.readSize().x/static_cast<double>(n.x);
	r.y=System.readSize().y/static_cast<double>(n.y);
	
	std::vector<double> heightMap;
	std::vector<int> countMap;
	for(int i=0;i<n.x*n.y;i++)
	{
		heightMap.push_back(0);
		countMap.push_back(0);
	}
	
	position<double> *p=System.getPositions();
	for(int i=0;i<System.readNParticles();i++)
	{
		bool match=false;
		for(int tIndex=0;!match && tIndex<types.size();tIndex++)
			match=(p[i].type==types[tIndex]);
		
		if(match)
		{
			threeVector<int> c;
			c.x=floor(p[i].x/r.x);
			c.y=floor(p[i].y/r.y);
			int hash=c.x+c.y*n.x;
			heightMap[hash]+=p[i].z;
			countMap[hash]++;
		}
	}
	
	for(int i=0;i<heightMap.size();i++)
	{
		if(countMap[i]!=0)
			heightMap[i]/=static_cast<double>(countMap[i]);
		std::cout << heightMap[i] << std::endl;
		//if(i%n.y==n.y-1)
		//	std::cout << std::endl;
	}

	return 0;
}
