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
	if(argc<5)
	{
		std::cerr << argv[0] << " name fromTime dimension nDivisions type1 type2 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime=0;
	threeVector<int> n;
	std::string dChar;
	int nDivisions=1;
	cmdArg >> fromTime >> dChar >> nDivisions;
	
	int dimension=0;
	if(dChar=="0")
		dimension=0;
	else if(dChar=="1")
		dimension=1;
	else if(dChar=="2")
		dimension=2;
	else if(dChar=="X")
		dimension=0;
	else if(dChar=="Y")
		dimension=1;
	else if(dChar=="Z")
		dimension=2;
	else if(dChar=="x")
		dimension=0;
	else if(dChar=="y")
		dimension=1;
	else if(dChar=="z")
		dimension=2;
	else
	{
		std::cerr << "dimension must be one of 0, 1, 2, X, Y, Z, x, y, or z!" << std::endl;
		return -1;
	}

	std::vector<int> types;
	
	for(int i=5;i<argc;i++)
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
	
	double widthAvg=0,thickAvg=0;
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
			std::vector< std::vector< position<double> > > cellList(nDivisions*nDivisions);
			threeVector<double> cellSize=size/(nDivisions-1);
			double avgHeight=0, hCount=0;
			for(int i=0;i<System.readNParticles();i++)
			{
				bool match=false;
				for(int tIndex=0;tIndex<types.size();tIndex++)
					if(p[i].type==types[tIndex]) match=true;
				
				if(match)
				{
					threeVector<int> cellI;
					cellI.x=int(p[i].x/cellSize.x)%nDivisions;
					cellI.y=int(p[i].y/cellSize.y)%nDivisions;
					cellI.z=int(p[i].z/cellSize.z)%nDivisions;
					cellList[cellI.s[(dimension+1)%3]+cellI.s[(dimension+2)%3]*nDivisions].push_back(p[i]);
					avgHeight+=p[i].s[dimension];
					hCount++;
				}
			}
			if(hCount!=0)
				avgHeight/=hCount;
			double width1=0,width2=0;
			double thick=0;
			for(auto cell:cellList)
			{
				double height=0;
				double mx=0,mn=size.s[dimension];
				for(auto pi:cell)
				{
					width1+=std::pow(avgHeight-pi.s[dimension],2.0);
					if(pi.s[dimension]>mx) mx=pi.s[dimension];
					if(pi.s[dimension]<mn) mn=pi.s[dimension];
				}
				if(cell.size()>0)
				{
					
					//width1+=avgHeight-(height/cell.size());
					//width2+=std::pow(height/cell.size(),2.0);
					thick+=mx-mn;
				}
			}
			width1/=hCount;
			width1=std::sqrt(width1);
			//width2/=nDivisions*nDivisions;
			//width1=width2-width1*width1;
			thick/=nDivisions*nDivisions;
			frames++;
			std::cerr << time << std::endl;
			//std::cout << time << ' ' << width1 << ' ' << thick << std::endl;
			widthAvg+=width1;
			thickAvg+=thick;
		}
		time+=System.readStoreInterval();
	}
	std::cout << argv[1] << ' ' << (widthAvg/frames) << ' ' << (thickAvg/frames) << std::endl;
	
	return 0;
}

