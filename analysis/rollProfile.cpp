//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc!=2 && argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cout << "To see available molecules:\n";
		std::cout << "\tusage: " << argv[0] << " name\n";
		
		std::cout << "To get the roll profile of a molecule at molIndex:\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex\n";
		
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	molecule<double, fourVector<int> > *m=System.getMolecule();
	
	if(argc==2)
	{
		for(int i=0;i<System.readNMolecules();i++)
		{
			if(m[i].readType()==CHAIN)
			{
				int length=m[i].getBonds()[0].s[CHAINLENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Chain " << i << " with " << nChains << " chains of length " << length << ".\n";
			}
			if(m[i].readType()==BOND)
			{
				int length=m[i].getBonds()[0].s[CHAINLENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Bond " << i << " with " << m[i].readNBond() << " bonds.\n";
			}
			if(m[i].readType()==BEND)
			{
				int length=m[i].getBonds()[0].s[CHAINLENGTH];
				int nChains=m[i].getBonds()[0].s[NCHAINS];
				std::cout << "Bend " << i << " with " << m[i].readNBond() << " bends.\n";
			}
		}
	}
	
	if(argc==3)
	{
		xyzFormat<double> xyzFile(p,nParticles);
		xyzFile.open(framesName.c_str(), std::ios::in);
		threeVector<double> size=System.readSize();
		int molIndex;
		std::stringstream cmdArg;
		for(int i=2;i<argc;i++)
			cmdArg << argv[i] << '\t';
		cmdArg >> molIndex;
		
		if(molIndex>=System.readNMolecules() || molIndex<0)
		{
			std::cerr << "Molecule index " << molIndex << " is out of bounds [0, " << System.readNMolecules() << ").\n";
			return 0;
		}
		
		double time=0;
		
		std::vector<int> pIndices;
		
		for(int bond=0;bond<m[molIndex].readNBond();bond++)
		{
			switch(m[molIndex].readType())
			{
				case CHAIN:
				{
					int start=m[molIndex].getBonds()[bond].s[START];
					int length=m[molIndex].getBonds()[bond].s[CHAINLENGTH];
					int nChains=m[molIndex].getBonds()[bond].s[NCHAINS];
					
					for(int i=start;i<start+nChains*length;i+=length)
						for(int j=i;j<i+length-1;j++)
							pIndices.push_back(j);
				}
				break;
				case BOND:
				{
					pIndices.push_back(m[molIndex].getBonds()[bond].s[0]);
					pIndices.push_back(m[molIndex].getBonds()[bond].s[1]);
				}
				break;
				case BEND:
				{
					pIndices.push_back(m[molIndex].getBonds()[bond].s[0]);
					pIndices.push_back(m[molIndex].getBonds()[bond].s[1]);
					pIndices.push_back(m[molIndex].getBonds()[bond].s[2]);
				}
				break;
			}
		}
		std::sort(pIndices.begin(), pIndices.end());
		std::vector<int>::iterator it;
		it=std::unique(pIndices.begin(), pIndices.end());
		pIndices.resize(std::distance(pIndices.begin(),it-1));
		
		std::vector< threeVector<double> > lastPosition;
		threeVector<double> lastComM;
		
		if(xyzFile.load())
		{
			for(int i=0;i<pIndices.size();i++)
			{
				threeVector<double> lP;
				lP.x=p[pIndices[i]].x;
				lP.y=p[pIndices[i]].y;
				lP.z=p[pIndices[i]].z;
				
				if(i!=0)
				{
					threeVector<double> d;
					d.x=lP.x-lastPosition[0].x;
					d.y=lP.y-lastPosition[0].y;
					d.z=lP.z-lastPosition[0].z;
					lP.x-=(d.x>size.x/2.0)?size.x:0;
					lP.y-=(d.y>size.y/2.0)?size.y:0;
					lP.z-=(d.z>size.z/2.0)?size.z:0;
					lP.x+=(d.x<-size.x/2.0)?size.x:0;
					lP.y+=(d.y<-size.y/2.0)?size.y:0;
					lP.z+=(d.z<-size.z/2.0)?size.z:0;
				}
				lastComM.x+=lP.x;
				lastComM.y+=lP.y;
				lastComM.z+=lP.z;
				lastPosition.push_back(lP);
			}
			lastComM.x/=static_cast<double>(pIndices.size());
			lastComM.y/=static_cast<double>(pIndices.size());
			lastComM.z/=static_cast<double>(pIndices.size());
		}
		
		
		std::string dataRollName("roll_");
		dataRollName+=name;
		dataRollName+=".dat";
		
		std::string dataMSDName("MSD_");
		dataMSDName+=name;
		dataMSDName+=".dat";
		
		std::fstream rollFile, MSDFile;
		rollFile.open(dataRollName.c_str(), std::ios::out);
		MSDFile.open(dataMSDName.c_str(), std::ios::out);
		
		threeVector<double> traj=lastComM;
		std::vector< threeVector<double> > trajList;
		
		while(xyzFile.load())
		{
			std::vector< threeVector<double> > translated;
			threeVector<double> comM;
			for(int i=0;i<pIndices.size();i++)
			{
				threeVector<double> lP;
				lP.x=p[pIndices[i]].x;
				lP.y=p[pIndices[i]].y;
				lP.z=p[pIndices[i]].z;
				
				if(i!=0)
				{
					threeVector<double> d;
					d.x=lP.x-translated[0].x;
					d.y=lP.y-translated[0].y;
					d.z=lP.z-translated[0].z;
					lP.x-=(d.x>size.x/2.0)?size.x:0;
					lP.y-=(d.y>size.y/2.0)?size.y:0;
					lP.z-=(d.z>size.z/2.0)?size.z:0;
					lP.x+=(d.x<-size.x/2.0)?size.x:0;
					lP.y+=(d.y<-size.y/2.0)?size.y:0;
					lP.z+=(d.z<-size.z/2.0)?size.z:0;
				}
				comM.x+=lP.x;
				comM.y+=lP.y;
				comM.z+=lP.z;
				translated.push_back(lP);
			}
			comM.x/=static_cast<double>(pIndices.size());
			comM.y/=static_cast<double>(pIndices.size());
			comM.z/=static_cast<double>(pIndices.size());
			
			threeVector<double> theta(0);
			
			for(int i=0;i<pIndices.size();i++)
			{
				threeVector<double> tCenter, lCenter, d0,d1;
				tCenter.x=translated[i].x-comM.x;
				tCenter.y=translated[i].y-comM.y;
				tCenter.z=translated[i].z-comM.z;
				//tCenter=unitVector(tCenter);
				
				lCenter.x=lastPosition[i].x-lastComM.x;
				lCenter.y=lastPosition[i].y-lastComM.y;
				lCenter.z=lastPosition[i].z-lastComM.z;
				//lCenter=unitVector(lCenter);
				
				//d.x=asin((tCenter.z*lCenter.y-lCenter.z*tCenter.y)/(lCenter.y*lCenter.y+lCenter.z*lCenter.z));
				//d.y=asin((lCenter.z*lCenter.x-tCenter.z*lCenter.x)/(lCenter.x*lCenter.x+lCenter.z*lCenter.z));
				//d.z=asin((tCenter.z*lCenter.x-tCenter.x*lCenter.y)/(lCenter.x*lCenter.x+lCenter.y*lCenter.y));
				
				//d.x=(tCenter.z*lCenter.y-lCenter.z*tCenter.y)/(lCenter.y*lCenter.y+lCenter.z*lCenter.z);
				//d.y=(lCenter.z*lCenter.x-tCenter.z*lCenter.x)/(lCenter.x*lCenter.x+lCenter.z*lCenter.z);
				//d.z=(tCenter.z*lCenter.x-tCenter.x*lCenter.y)/(lCenter.x*lCenter.x+lCenter.y*lCenter.y);
				
				//theta.x+=tCenter.x-lCenter.x;
				//theta.y+=tCenter.y-lCenter.y;
				//theta.z+=tCenter.z-lCenter.z;
				
				//std::cout << d.x << '\t' << d.y << '\t' << d.z << std::endl;
				
				d0.x=0;
				d0.y=lCenter.y;
				d0.z=lCenter.z;
				d0=unitVector(d0);
				d1.x=0;
				d1.y=tCenter.y;
				d1.z=tCenter.z;
				d1=unitVector(d1);
				theta.x+=(d1.z*d0.y-d0.z*d1.y)/(d0.y*d0.y+d0.z*d0.z);
				
				d0.x=lCenter.x;
				d0.y=0;
				d0.z=lCenter.z;
				d0=unitVector(d0);
				d1.x=tCenter.x;
				d1.y=0;
				d1.z=tCenter.z;
				d1=unitVector(d1);
				theta.y+=(d0.z*d1.x-d1.z*d0.x)/(d0.x*d0.x+d0.z*d0.z);
				
				d0.x=lCenter.x;
				d0.y=lCenter.y;
				d0.z=0;
				d0=unitVector(d0);
				d1.x=tCenter.x;
				d1.y=tCenter.y;
				d1.z=0;
				d1=unitVector(d1);
				theta.z+=(d1.y*d0.x-d1.x*d0.y)/(d0.x*d0.x+d0.y*d0.y);
				
			}
			
			theta.x/=static_cast<double>(pIndices.size());
			theta.y/=static_cast<double>(pIndices.size());
			theta.z/=static_cast<double>(pIndices.size());
			
			
			threeVector<double> offSetCOM;
			offSetCOM.x=comM.x-lastComM.x;
			offSetCOM.y=comM.y-lastComM.y;
			offSetCOM.z=comM.z-lastComM.z;
			offSetCOM.x-=(offSetCOM.x>size.x/2.0)?size.x:0;
			offSetCOM.y-=(offSetCOM.y>size.y/2.0)?size.y:0;
			offSetCOM.z-=(offSetCOM.z>size.z/2.0)?size.z:0;
			offSetCOM.x+=(offSetCOM.x<-size.x/2.0)?size.x:0;
			offSetCOM.y+=(offSetCOM.y<-size.y/2.0)?size.y:0;
			offSetCOM.z+=(offSetCOM.z<-size.z/2.0)?size.z:0;
			rollFile << time << '\t' << theta.x << '\t' << theta.y << '\t' << theta.z << '\t';
			rollFile << offSetCOM.x << '\t' << offSetCOM.y << '\t' << offSetCOM.z << std::endl;
			traj.x+=offSetCOM.x;
			traj.y+=offSetCOM.y;
			traj.z+=offSetCOM.z;
			trajList.push_back(traj);
			//MSDFile << time << '\t' << MSD.x << '\t' << MSD.y << '\t' << MSD.z << std::endl;
			
			
			lastPosition=translated;
			lastComM=comM;
			time+=System.readStoreInterval();
			std::cerr << time << std::endl;
		}
		//for(int i=0;i<trajList.size();i++)
		//	trajList[i]-=trajList[0];
		//std::vector<double> MSD;
		for(int skip=1;skip<trajList.size();skip++)
		{
			double avg=0, nAvg=0;
			for(int i=0;i<trajList.size()-skip;i++)
			{
				threeVector<double> d;
				d.x=trajList[i+skip].x-trajList[i].x;
				d.y=trajList[i+skip].y-trajList[i].y;
				d.z=trajList[i+skip].z-trajList[i].z;
				avg+=(d.x*d.x+d.y*d.y+d.z*d.z);
				nAvg++;
			}
			avg/=nAvg;
			MSDFile << static_cast<double>(skip)*System.readStoreInterval() << '\t' << avg << std::endl;
		}
		MSDFile.close();
		rollFile.close();
	}
	/*
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	while(xyzFile.load())
	{
		//Do data measurements in here
		int a=0;
		for(int i=0;i<nParticles;i++)
			if(p[i].x>System.readSize().x/2.0)
				a++;
		std::cout << a << std::endl;
	}
	*/
	return 0;
}

