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
		std::cout << argv[0] << " name molId sliceSize\n";
		return 0;
	}

	std::stringstream cmdArg;
	cmdArg << argv[2] << ' ' << argv[3];
	double sliceSize;
	int molId;
	cmdArg >> molId >> sliceSize;

	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//size file
	std::string sizeName("size_");
	sizeName+=argv[1];
	sizeName+=".dat";
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	
	threeVector<double> size=System.readSize();
	int nTypes=System.readNTypes();
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	double time=0;
	while(xyzFile.load())
	{
		if(sizeFile.is_open())
		{
			double sTime=-1;
			threeVector<double> sCurrent;
			while(sTime<time-System.readMeasureInterval()+System.readDeltaT() && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
			time=sTime;
		}
		std::map<int,std::vector<int>> beads;
		std::map<int,double> profile;
		
		for(int i=0;i<nParticles;i++)
		{
			int slice=static_cast<int>(p[i].z/sliceSize);
			beads[slice].push_back(i);
			profile[slice]=0;
		}
		
		molecule<double, fourVector<int> > *m=System.getMolecule();
		
		//this checks to see if it is actually a beed type
		if(m[molId].readType()==BEAD)
		{
			//x=index, y=bead-bead type, z=bead-particle type, t=radius
			std::vector<threeVector<int> > beadInfo;
			std::vector<double> beadRadius;
			
			for(int k=0;k<m[molId].readNBond();k++)
			{
				threeVector<int> info;
				info.x=m[molId].getBonds()[k].x;
				info.y=0;//m[i].getConstants()[1];
				info.z=p[info.x].type;
				//this is actually repeated across all particles
				beadRadius.push_back(m[molId].getConstants()[BEADRADIUS]);
				//std::cout << info.t << std::endl;
				beadInfo.push_back(info);
			}
			
			double *C=m[molId].getConstants();	
			
			double cutoffSquaredRad=m[molId].getConstants()[0];//radius+rc
			cutoffSquaredRad*=cutoffSquaredRad;
			
			//for all other particles
			for(int j=0;j<m[molId].readNBond();j++)
			{
				for(auto slice:beads)
				{
					double nBeadsInSlice=0;
					for(auto k:slice.second)
					{
						threeVector<double> d;
						int first=beadInfo[j].x;
						int second=k;
						if(first!=second)
						{
							d.x=p[first].x-p[second].x;
							if(d.x>size.x/2.0) d.x-=size.x;
							if(d.x<-size.x/2.0) d.x+=size.x;
							d.y=p[first].y-p[second].y;
							if(d.y>size.y/2.0) d.y-=size.y;
							if(d.y<-size.y/2.0) d.y+=size.y;
							d.z=p[first].z-p[second].z;
							if(d.z>size.z/2.0) d.z-=size.z;
							if(d.z<-size.z/2.0) d.z+=size.z;
							
							int cindex=nBEADCONST*((beadInfo[j].z*nTypes)+p[k].type);
							
							double potential=beadPotential(d,&(C[cindex]), cutoffSquaredRad);
							if(potential!=0)
							{
								profile[slice.first]+=potential;
								nBeadsInSlice++;
							}
						}
					}
					//if(nBeadsInSlice!=0)
					//	profile[slice.first]/=nBeadsInSlice;
				}
			}
		}
		else
		{
			std::cerr << molId << " is not bead type!" << std::endl;
			return 1;
		}
		
		for(auto slice:profile)
		{
			std::cout << ((static_cast<double>(slice.first)+0.5)*sliceSize) << '\t' << slice.second << '\n';
		}
		std::cout << '\n';
		time+=System.readStoreInterval();
		
	}
	return 0;
}
