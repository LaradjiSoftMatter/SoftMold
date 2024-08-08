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
	if(argc!=5)
	{
		std::cout << argv[0] << " name fromTime centerType sliceSize\n";
		return 0;
	}

	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	int centerType;
	double sliceSize,fromTime;
	cmdArg >> fromTime >> centerType >> sliceSize;

	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//this is preallocated
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	//This one is now writing to the 'p' pointer.
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	//grab sizes
	threeVector<double> size=System.readSize();
	
	double time=0;
	
	std::string newName("size_");
	newName+=argv[1];
	newName+=".dat";
		
	std::fstream sizeFile;
	sizeFile.open(newName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		if(sTime>-1)
			time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	std::vector< std::vector<double> > profile;
	for(int i=0;i<System.readNTypes();i++)
		profile.push_back(std::vector<double>());
	
	int frames=0;
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		
		if(time>=fromTime)
		{
			std::vector< std::vector< position<double> > > cellList;
			threeVector<int> n;
			n.x=floor(System.readSize().x/System.readCutoff());
			n.y=floor(System.readSize().y/System.readCutoff());
			n.z=floor(System.readSize().z/System.readCutoff());
			
			threeVector<double> r;
			r.x=System.readSize().x/static_cast<double>(n.x);
			r.y=System.readSize().y/static_cast<double>(n.y);
			r.z=System.readSize().z/static_cast<double>(n.z);
		
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
			frames++;
		}
		time+=System.readStoreInterval();
	}
	
	for(int i=0;i<profile.size();i++)
	{
		for(int j=0;j<profile[i].size();j++)
			std::cout << (static_cast<double>(j)*sliceSize) << '\t' << profile[i][j]/frames << '\n';
		std::cout << '\n';
	}

	return 0;
}
