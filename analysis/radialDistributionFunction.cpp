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
	if(argc!=6 && argc!=5)
	{
		std::cerr << "Usage: " << argv[0] << " name fromTime deltaR typeCenter typeDist \n";
		std::cerr << "Usage: " << argv[0] << " name fromTime deltaR typeCenter \n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime,deltaR;
	int typeCenter,typeDist;
	cmdArg >> fromTime >> deltaR >> typeCenter;
	if(argc==6)
		cmdArg >> typeDist;
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//set time
	double time=0;

	//grab sizes
	threeVector<double> size=System.readSize();
	
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
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
		
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	int frames=0;
	std::vector<double> dist;
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
			position<double> s;
			s.x=size.x;
			s.y=size.y;
			s.z=size.z;
			std::vector<int> typeCenterInd,typeDistInd;
			for(int i=0;i<nParticles;i++)
			{
				if(p[i].type==typeCenter)
					typeCenterInd.push_back(i);
				if(p[i].type==typeDist && argc==6)
					typeDistInd.push_back(i);
			}
			
			if(argc==6)
				for(const auto &i:typeCenterInd)
				for(const auto &j:typeDistInd)
				{
					int dInd=distance(p[i],p[j],s)/deltaR;
					while(dist.size()<dInd) dist.push_back(0);
					dist[dInd]++;
				}
			else
				for(auto it=typeCenterInd.begin();it!=typeCenterInd.end();it++)
				for(auto jt=it+1;jt!=typeCenterInd.end();jt++)
				{
					int dInd=distance(p[*it],p[*jt],s)/deltaR;
					while(dist.size()<dInd) dist.push_back(0);
					dist[dInd]++;
				}
				
			frames++;
			std::cerr << time << std::endl;
		}
		time+=System.readStoreInterval();
	}
	double count=0;
	for(const auto& d:dist)
		count+=d;
	for(int i=0;i<dist.size();i++)
	{
		double r=deltaR*(i+0.5);
		double v=4.0*M_PI*(pow(r+deltaR/2.0,3.0)-pow(r-deltaR/2.0,3.0))/3.0;
		std::cout << r << '\t' << dist[i]/(v*frames) << std::endl;
	}
		
	return 0;
}

