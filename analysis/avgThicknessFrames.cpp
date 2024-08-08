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
		std::cerr << argv[0] << " name fromTime dimension cutoff type1 type2 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime=0,cutoff=0;
	threeVector<int> n;
	std::string dChar;
	cmdArg >> fromTime >> dChar >> cutoff;
	
	if(cutoff<=0)
	{
		std::cerr << "cutoff must be a floating number greater than zero!" << std::endl;
		return -1;
	}
	
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
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		time=sTime;
		size=sCurrent;
	}
		
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<double> allWidth;
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
		}
		
		if(time>=fromTime)
		{
			threeVector<double> cs;
			threeVector<int> ns;
			ns.x=size.x/cutoff;
			ns.y=size.x/cutoff;
			cs.x=size.x/ns.x;
			cs.y=size.y/ns.y;
			std::vector< std::vector< position<double> > > cellList(ns.x*ns.y);
			
			for(int i=0;i<System.readNParticles();i++)
			{
				bool match=false;
				for(int tIndex=0;!match && tIndex<types.size();tIndex++)
					match=(p[i].type==types[tIndex]);
				
				if(match)
				{
					threeVector<int> ni;
					ni.x=p[i].s[(dimension+1)%3];
					ni.x=p[i].s[(dimension+2)%3];
					cellList[ni.x+ni.y*ns.x].push_back(p[i]);
				}
			}
			
			double totalAvg=0,totalCount=0;
			for(auto &cell:cellList)
			{
				double avg=0,avg2=0,count=0;
			
				for(auto &particle:cell)
				{
					avg+=particle.s[dimension];
					avg2+=particle.s[dimension]*particle.s[dimension];
					count++;
				}
				if(count!=0) avg/=count;
				if(count!=0) avg2/=count;
				totalAvg+=std::sqrt(avg2-avg*avg);
				if(count!=0) totalCount++;
			}
			
			if(totalCount!=0) totalAvg/=totalCount;
			
			allWidth.push_back(totalAvg);
			std::cerr << time << std::endl;
		}
		time+=System.readStoreInterval();
	}
	double avgWidth=0;
	for(auto width:allWidth)
		avgWidth+=width;
	if(allWidth.size()>0) avgWidth/=allWidth.size();
	
	double widthStd=0;
	for(auto width:allWidth)
		widthStd+=std::pow(width-avgWidth,2);
	if(allWidth.size()>0)
		widthStd/=allWidth.size();
	widthStd=std::sqrt(widthStd);
	
	std::cout << avgWidth << '\t' << widthStd << std::endl;
	
	return 0;
}

