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
	if(argc<4)
	{
		std::cout << argv[0] << " name pIndex1 pIndex2 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	std::vector<int> pIndex;
	
	for(int i=2;i<argc;i++)
	{
		int pIndexVal;
		cmdArg >> pIndexVal;
		pIndex.push_back(pIndexVal);
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
	
		//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string newName("size_");
	newName+=argv[1];
	newName+=".dat";
	
	double time=0;
	
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
	
	int frame=0;
	
	while(xyzFile.load())
	{
		std::vector< std::vector< position<double> > > cellList;
		time=System.readStoreInterval()*static_cast<double>(frame);
		
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
			{
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				//std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
			size=sCurrent;
		}
		
		for(int i=0;i<pIndex.size();i++)
		{
			for(int j=i+1;j<pIndex.size();j++)
			{
				threeVector<double> d;
				d.x=p[pIndex[i]].x-p[pIndex[j]].x;
				d.y=p[pIndex[i]].y-p[pIndex[j]].y;
				d.z=p[pIndex[i]].z-p[pIndex[j]].z;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.x<-size.x/2.0) d.x+=size.x;
				if(d.y<-size.y/2.0) d.y+=size.y;
				if(d.z<-size.z/2.0) d.z+=size.z;
				std::cout << time << '\t' << sqrt(d.x*d.x+d.y*d.y+d.z*d.z) << std::endl;
			}
		}
		std::cerr << "Time: " << time << '\t' << "Frame: " << frame << std::endl;
		frame++;
	}

	return 0;
}
