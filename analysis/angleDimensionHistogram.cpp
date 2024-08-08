/** \brief Loads each frame of a frames file and outputs frames in a new frames file.
 */

#include "../include/MD.h"
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=7)
	{
		std::cerr << "\tUsage: " << argv[0] << " name fromTime thickness dimension molIndex dCos\n";
		std::cerr << "Calculates bond cos(angle) for slabs along a dimension\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int mIndex;
	double thickness,dCos,fromTime;
	std::string dimStr;
	cmdArg >> fromTime >> thickness >> dimStr >> mIndex >> dCos;
	int dim=-1;
	if(dimStr=="x" || dimStr=="X" || dimStr=="0") dim=0;
	if(dimStr=="y" || dimStr=="Y" || dimStr=="1") dim=1;
	if(dimStr=="z" || dimStr=="Z" || dimStr=="2") dim=2;
	if(dim==-1)
	{
		std::cerr << "dimension must be one of {x,y,z,X,Y,Z,0,1,2}!" << std::endl;
		return -1;
	}
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	if(System.readNMolecules()<=mIndex)
	{
		std::cerr << "molIndex out of range!" << std::endl;
		std::cerr << "\tSystem: " << System.readNMolecules() << '\t' << "molIndex: " << mIndex << std::endl;
		return -1;
	}
	
	molecule< double, fourVector<int> > *m=System.getMolecule();
	if(!(m[mIndex].readType()==CHAIN || m[mIndex].readType()==BOND))
	{
		std::cerr << "Wrong type!" << std::endl;
		std::cerr << "Should be one of CHAIN or BOND types (" << CHAIN << " or " << BOND << ")!" << std::endl;
		return -1;
	}
	
	std::vector<fourVector<int>> bondPairs;
	if(m[mIndex].readType()==BOND)
		while(bondPairs.size()<m[mIndex].readNBond())
			bondPairs.emplace_back(m[mIndex].getBonds()[bondPairs.size()]);
	else for(int i=0;i<m[mIndex].readNBond();i++)
	{
		int start=m[mIndex].getBonds()[i].s[START];
		int nChains=m[mIndex].getBonds()[i].s[NCHAINS];
		int chainLength=m[mIndex].getBonds()[i].s[CHAINLENGTH];
		int end=start+nChains*chainLength;
		for(int j=start;j<end;j+=chainLength)
		{
			for(int k=0;k<chainLength-1;k++)
			{
				fourVector<int> bond;
				bond.s[0]=k+j;
				bond.s[1]=k+j+1;
				bondPairs.emplace_back(bond);
			}
		}
	}
	
	//Get our size file
	std::string sizeName("size_");
	sizeName+=name;
	sizeName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent=System.readSize();
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		System.setSize(sCurrent);
	}
	
	//Get our frames file
	std::string frameName("frames_");
	frameName+=name;
	frameName+=".xyz";
	position<double> *p=System.getPositions();
	int nParticles=0;
	
	//This one is now handling allocation of the 'p' pointer.
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(frameName.c_str(), std::ios::in);
	
	double time=0;
	int frame=0;
	
	//start reading frames
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=System.readSize();//just in case it is current
			while(sTime<time && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			System.setSize(sCurrent);
		}
		if(time>=fromTime)
		{
		threeVector<double> size=System.readSize();
		std::vector<std::vector<double>> slabs;
		threeVector<double> zHat=0;
		zHat.s[dim]=1.0;
		for(auto bP:bondPairs)
		{
			threeVector<double> d;
			d.x=p[bP.s[1]].x-p[bP.s[0]].x;
			d.y=p[bP.s[1]].y-p[bP.s[0]].y;
			d.z=p[bP.s[1]].z-p[bP.s[0]].z;
			threeVector<double> offset=0;
			if(d.x>size.x/2.0)offset.x=-size.x;
			if(d.y>size.y/2.0)offset.y=-size.y;
			if(d.z>size.z/2.0)offset.z=-size.z;
			if(d.x<-size.x/2.0)offset.x=size.x;
			if(d.y<-size.y/2.0)offset.y=size.y;
			if(d.z<-size.z/2.0)offset.z=size.z;
			d.x+=offset.x;
			d.y+=offset.y;
			d.z+=offset.z;
			d=unitVector(d);
			threeVector<double> com;
			com.x=(p[bP.s[0]].x+p[bP.s[1]].x+offset.x)/2.0;
			com.y=(p[bP.s[0]].y+p[bP.s[1]].y+offset.y)/2.0;
			com.z=(p[bP.s[0]].z+p[bP.s[1]].z+offset.z)/2.0;
			int slabIndex=std::floor(com.s[dim]/thickness);
			//double theta=std::acos(dotProduct(d,zHat));
			double theta=1+dotProduct(d,zHat);
			while(slabs.size()<=slabIndex)
				slabs.push_back(std::vector<double>());
			slabs[slabIndex].push_back(theta);
		}
		//These are between 0 and pi, not 0 and 2*pi, so no need for circular averages
		for(int i=0;i<slabs.size();i++)
		{
			if(slabs[i].size()>0)
			{
				std::vector<double> histogram;
				for(auto a:slabs[i])
				{
					int hIndex=a/dCos;
					while(histogram.size()<hIndex+1)
						histogram.push_back(0);
					histogram[hIndex]+=1.0/(slabs[i].size()*dCos);
				}
				
				int angleIndex=0;
				for(auto h:histogram)
				{
					std::cout << ((0.5+angleIndex)*dCos)-1 << '\t' << h << '\t' << i*thickness << '\t' << time << '\n';
					angleIndex++;
				}
				std::cout << std::endl;
			}
		}
		}
		std::cerr << "Read frame " << frame << " with time " << time << " tau." << std::endl;
		
		//update time for next frame
		time+=System.readStoreInterval();
		
		//update frame counter
		frame++;
	}
	//close our old files
	xyzFile.close();
	
	return 0;
}

