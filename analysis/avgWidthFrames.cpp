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
		std::cerr << argv[0] << " name fromTime dimension type1 type2 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime=0;
	threeVector<int> n;
	std::string dChar;
	cmdArg >> fromTime >> dChar;
	
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
	
	for(int i=4;i<argc;i++)
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
	
	double widthAvg=0;
	std::vector<double> allWidth;
	while(xyzFile.load())
	{
		std::vector< std::vector< position<double> > > cellList;
		
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
		if(time>=fromTime)
		{
			//initialize avgs
			double avg=0,avg2=0,count=0;
			
			for(int i=0;i<System.readNParticles();i++)
			{
				bool match=false;
				for(int tIndex=0;!match && tIndex<types.size();tIndex++)
					match=(p[i].type==types[tIndex]);
				
				if(match)
				{
					avg+=p[i].s[dimension];
					avg2+=p[i].s[dimension]*p[i].s[dimension];
					count++;
				}
			}
			
			if(count!=0)
				avg/=count;
			if(count!=0)
				avg2/=count;
			
			double width=std::sqrt(avg2-std::pow(avg,2));
			widthAvg+=width;
			allWidth.push_back(width);
			frames++;
			std::cerr << time << std::endl;
			//std::cout << time << ' ' << width << std::endl;
		}
		time+=System.readStoreInterval();
	}
	if(frames>0)
		widthAvg/=frames;
	double widthStd=0;
	for(auto width:allWidth)
		widthStd+=std::pow(width-widthAvg,2);
	if(frames>0)
		widthStd/=frames;
	widthStd=std::sqrt(widthStd);
	std::cout << widthAvg << '\t' << widthStd << std::endl;
	
	return 0;
}

