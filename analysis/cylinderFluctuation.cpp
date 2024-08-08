/** \brief Loads each frame of a frames file and outputs frames in a new frames file.
 */

#include "../include/MD.h"
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=6)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "outputs fluctuation of cylinder composed of flucType about comType along a cylinder through a plane\n";
		std::cerr << "\tusage: " << argv[0] << " framesFile.xyz comType flucType radius plane\n";
		std::cerr << "\tplane: 0=X, 1=Y, 2=Z\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	int comType, flucType, plane;
	double radius;
	cmdArg >> comType >> flucType >> radius >> plane;
	
	std::cerr << comType << '\t' << flucType << '\t' << radius << '\t' << plane << std::endl;
	
	//Get our frames file
	position<double> *p=NULL;
	int nParticles=0;
	
	//This one is now handling allocation of the 'p' pointer.
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(name, std::ios::in);
	
	//Make sure there is something in the file
	
	time_t start=time(NULL);
	int frame=0;
	
	std::vector<int> comIndex;
	
	double histSpacing=0.5;
	
	while(xyzFile.load())
	{
		p=xyzFile.getPositions();
		nParticles=xyzFile.readNParticles();
		std::cerr << "Frame: " << frame << std::endl;
		if(comIndex.size()==0)
		{
			for(int i=0;i<nParticles;i++)
				if(p[i].type==comType)
					comIndex.push_back(i);
			std::cerr << nParticles << '\t' << comIndex.size() << std::endl;
		}
		threeVector<double> com(0.0);
		for(int j=0;j<comIndex.size();j++)
		{
			int i=comIndex[j];
			com.x+=p[i].x;
			com.y+=p[i].y;
			com.z+=p[i].z;
		}
		if(comIndex.size()>0)
			com/=static_cast<double>(comIndex.size());
		double flucCenter=0;
		std::vector<double> flucHist;
		int flucCount=0;
		double radSquared=radius*radius;
		for(int i=0;i<nParticles;i++)
		{
			if(p[i].type==flucType)
			{
				int histIndex=p[i].s[plane]/histSpacing;
				threeVector<double> d;
				d.x=p[i].x-com.x;
				d.y=p[i].y-com.y;
				d.z=p[i].z-com.z;
				double dr=0;
				for(int j=0;j<3;j++)
					if(j!=plane)
						dr+=d.s[j]*d.s[j];
				if(dr<radSquared)
				{
					while(flucHist.size()<=histIndex)
						flucHist.push_back(0);
					flucHist[histIndex]++;
					flucCount++;
					flucCenter+=p[i].s[plane];
				}
			}
		}
		int center=0;
		double sum=0;
		if(flucCount>0)
			flucCenter/=static_cast<double>(flucCount);
		for(int i=0;i<flucHist.size() && flucCount>0;i++)
		{
			flucHist[i]/=static_cast<double>(flucCount)*histSpacing;
			sum+=flucHist[i]*histSpacing;
			if(sum>=0.5 && center==0)
				center=i;
		}
		for(int i=0;i<flucHist.size() && flucCount>0;i++)
			std::cout << static_cast<double>(i)*histSpacing-flucCenter << '\t' << flucHist[i] << std::endl;
		std::cout << std::endl;
		frame++;
	}
	time_t end=time(NULL);
	std::cerr << "Loads at " << static_cast<double>(frame)/static_cast<double>(end-start) << " frames per second ";
	std::cerr << " in " << end-start << " seconds!\n";
	//close our old files
	xyzFile.close();
	
	return 0;
}

