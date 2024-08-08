//Loading an MD simulation "using systemMD.h".
//Everything needed to do a simulation is already included in systemMD.h.
//Even standard libraries like iostream and fstream are inlcuded.
//But, if any extra libraries are needed you can include them here or in systemMD.h.
//Permenent libraries (included with standard library or with this package) should
// go in systemMD.h, and anything you need that is external should go in this file.
#include <fenv.h>
//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"
#include "../include/algorithms/threeHash.h"
#include "../include/algorithms/blockKey.h"

#include <deque>
#include <list>
#include <chrono>
#include <unordered_map>

#define sThreshold 20
#define pSig 0.3

struct edgeType {
	int x,y;
};

//for sorting keys, logarithmic search
bool compByX(edgeType a, edgeType b)
{
	//if(a.x<b.x)
	//	return true;
	//return false;
	return a.x<b.x;
}

//for sorting keys, logarithmic search
bool compByT(fourVector<double> a, fourVector<double> b)
{
	//if(a.x<b.x)
	//	return true;
	//return false;
	return a.t<b.t;
}

//for sorting keys, logarithmic search
bool equByXY(edgeType a, edgeType b)
{
	return (a.x==b.x) && (a.y==b.y);
}


//for simple linear search
class equByY
{
	public:
		bool operator () (edgeType a)
		{
			return a.y==y;
		};
		int x;
		int y;
};

//for simple linear search
class equByX
{
	public:
		bool operator () (edgeType a)
		{
			return a.x==x;
		};
		int x;
		int y;
};

template<typename T>
struct lipid {
	threeVector<T> orient;
	position<T> pos;
	int index;
	
	bool operator== (lipid& other)
	{
		return index==other.index;
	}
	bool operator!= (lipid& other)
	{
		return index!=other.index;
	}
};




void cmdHelp()
{
	std::cerr << "Options:\n";
	std::cerr << "--surfSpacing [float]\n\tMinimum spacing between surfaces. Default 0!\n";
	std::cerr << "--nSurfaces [int]\n\tDefault -1 (basically traverse everything)!\n";
	std::cerr << "--minSurfaceSize [int]\n\tMinimum surface size. Default 20!\n";
	std::cerr << "--maxVariance [float]\n\tMaximum variance of surfaces size. Default 0.15!\n";
	std::cerr << "--resurf [int]\n\tRecompute minimum surface spacing every nth frame. Default 20!\n";
	std::cerr << "--frame [float]\n\tWhich frame. Default all!\n";
	std::cerr << "--vertSpacing [float]\n\tSpacing for rotation mesh, patchSize=vertSpacing*sqrt(2.0)/3.0. Required!\n";
	std::cerr << "--nVertSpacing [int]\n\tNumber of rotation mesh tests. Default 1!\n";
	std::cerr << "--type [int]\n\tParticle type. Required!\n";
	std::cerr << "--zProfile [double]\n\tProfile the curvature along the z axis using [double] for slice size.\n";
	std::cerr << "--mIndex [int]\n\tMolecule index. defaults to all CHAIN types!\n";
	std::cerr << "--mCIndex [int] [float]\n\tMolecule contact index and cutoff. defaults none!\n";
	std::cerr << "--cType [int] [float]\n\tContact type and cutoff. defaults none!\n";
	std::cerr << "--initialTime [float]\n\tInitial simulation time (tau). Default 0!\n";
	std::cerr << "\tname is a simulation name not including the extension (.mpd)\n";
	std::cerr << "\tNote that --type defines the end to end vector head (for calculating normal).\n";
	std::cerr << "\tPossible files (# is patchSize):" << std::endl;
	std::cerr << "\t\tpatchProfile_S#_name.dat for multiple patch sizes." << std::endl;
	std::cerr << "\t\tgradSquare_SN#_name.dat for one patch size without contact." << std::endl;
	std::cerr << "\t\tgradSquare_SC#_name.dat for one patch size in contact." << std::endl;
	std::cerr << "\t\tcMap_S#_name.xyz for mapping." << std::endl;
}

int main(int argc, char **argv)
{
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW  | FE_UNDERFLOW);
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double surfSpacing=0, nSurfaces=-1, vertSpacing=-1, time=0,maxVariance=0.15,zSlice=0;
	int frame=-1, nVertSpacing=1, resurf=20, minSurfaceSize=20;
	std::string name;
	std::vector<int> types, mIndices, mCIndices, cTypes;
	std::vector<double> mCCutoffs,cCutoffs;
	
	//options
	while(!cmdArg.eof())
	{
		std::string cmd;
		cmdArg >> cmd;
		if(cmd.c_str()[0]=='-' && cmd.c_str()[1]=='-')
		{
			if((cmd=="--h") || (cmd=="--help"))
			{
				std::cerr << "Usage: " << argv[0] << " name ... Options\n";
				cmdHelp();
				return 0;
			}
			else if(cmd=="--surfSpacing")
				cmdArg >> surfSpacing;
			else if(cmd=="--maxVariance")
				cmdArg >> maxVariance;
			else if(cmd=="--nSurfaces")
				cmdArg >> nSurfaces;
			else if(cmd=="--minSurfaceSize")
				cmdArg >> minSurfaceSize;
			else if(cmd=="--resurf")
				cmdArg >> resurf;
			else if(cmd=="--frame")
				cmdArg >> frame;
			else if(cmd=="--zProfile")
				cmdArg >> zSlice;
			else if(cmd=="--vertSpacing")
				cmdArg >> vertSpacing;
			else if(cmd=="--nVertSpacing")
				cmdArg >> nVertSpacing;
			else if(cmd=="--initialTime")
				cmdArg >> time;
			else if(cmd=="--type")
			{
				int type;
				cmdArg >> type;
				types.push_back(type);
			}
			else if(cmd=="--mIndex")
			{
				int mIndex;
				cmdArg >> mIndex;
				mIndices.push_back(mIndex);
			}
			else if(cmd=="--cType")
			{
				int cType;
				double cCutoff;
				cmdArg >> cType;
				cmdArg >> cCutoff;
				cTypes.push_back(cType);
				if(cmdArg.eof())
				{
					std::cerr << cmd << " is expecting 2 arguments!" << std::endl;
					cmdHelp();
					return -1;
				}
				cCutoffs.push_back(cCutoff);
			}
			else if(cmd=="--mCIndex")
			{
				int mCIndex;
				double mCCutoff;
				cmdArg >> mCIndex;
				cmdArg >> mCCutoff;
				mCIndices.push_back(mCIndex);
				if(cmdArg.eof())
				{
					std::cerr << cmd << " is expecting 2 arguments!" << std::endl;
					cmdHelp();
					return -1;
				}
				mCCutoffs.push_back(mCCutoff);
			}
			else 
			{
				std::cerr << "Unrecognized option \"" << cmd << "\".\n\n";
				cmdHelp();
				return -1;
			}
		}
		else if(cmd.length()>0)
		{
			
			name=cmd;
			std::cerr << "Reading " << name << std::endl;
		}
	}
	
	if(vertSpacing<=0)
	{
		std::cerr << "--vertSpacing is required!" << std::endl;
		return -1;
	}
	double patchSize=vertSpacing*sqrt(2.0)/3.0;
	
	if(types.size()==0)
	{
		std::cerr << "--type is resuired!" << std::endl;
		return -1;
	}
	
	double surfSpacingSqr=surfSpacing*surfSpacing;
	double vertSpacingSqr=vertSpacing*vertSpacing;
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (name in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(name.c_str(),std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string sizeName("size_");
	sizeName+=name;
	sizeName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
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
	
	double beadPotential=-1;
	
	std::string newName("beadPotential_");
	newName+=name;
	newName+=".dat";
	
	std::fstream beadPotentialFile;
	beadPotentialFile.open(newName.c_str(), std::ios::in);
	if(beadPotentialFile.is_open())
	{
		double sTime=-1;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !beadPotentialFile.eof())
			beadPotentialFile >> sTime >> beadPotential;
		time=sTime;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//grab particles, by type, for surfaces
	std::vector<twoVector<int> > pIndex;
	
	for(int i=0;i<System.readNParticles();i++)
	{
		bool match=false;
		for(int tIndex=0;!match && tIndex<types.size();tIndex++)
			match=(p[i].type==types[tIndex]);
		
		if(match)
		{
			twoVector<int> buf;
			buf.x=i;
			buf.y=-1;
			pIndex.push_back(buf);
		}
	}
	
	//molecule chains end to end if needed
	for(int i=0;i<mIndices.size();i++)
	{
		if(mIndices[i]<System.readNMolecules() && mIndices[i]>=0)
		{
			molecule<double,fourVector<int> > m=System.getMolecule()[mIndices[i]];
			if(m.readType()==CHAIN)
			{
				fourVector<int> *bond=m.getBonds();
				
				for(int j=0;j<m.readNBond();j++)
				{
					int start=bond[j].s[START];
					int nChains=bond[j].s[NCHAINS];
					int length=bond[j].s[CHAINLENGTH];
					
					for(int k=0;k<pIndex.size();k++)
					{
						if(pIndex[k].x>=start && pIndex[k].x<start+nChains*length)
						{
							int chain=(pIndex[k].x-start)/length;
							int seg=(pIndex[k].x-start)%length;
							if(seg==0)//its at the 0 end, use the length-1 end
								pIndex[k].y=start+(length-1)+chain*length;
							else if(seg==length-1)//its at the length-1 end, use the 0 end
								pIndex[k].y=start+chain*length;
							else//its in the middle, just use the 0 end
								pIndex[k].y=start+chain*length;
						}
					}
				}
			}
			else
			{
				std::cerr << "Molecule wrong type! Must be CHAIN (" << CHAIN << ") type! #" 
				<< mIndices[i] << " != [0," 
				<< System.readNMolecules() << "]" << std::endl;
				return -1;
			}
		}
		else
		{
			std::cerr << "Molecule out of bounds! Must be CHAIN (" << CHAIN << ") type!"
			<< mIndices[i] << " != [0," 
			<< System.readNMolecules() << "]" << std::endl;
			throw 0;
		}
	}
	
	for(int i=0;mIndices.size()==0 && i<System.readNMolecules();i++)
	{
		molecule<double,fourVector<int> > m=System.getMolecule()[i];
		if(m.readType()==CHAIN)
		{
			fourVector<int> *bond=m.getBonds();
			
			for(int j=0;j<m.readNBond();j++)
			{
				int start=bond[j].s[START];
				int nChains=bond[j].s[NCHAINS];
				int length=bond[j].s[CHAINLENGTH];
				
				for(int k=0;k<pIndex.size();k++)
				{
					if(pIndex[k].x>=start && pIndex[k].x<start+nChains*length)
					{
						int chain=(pIndex[k].x-start)/length;
						int seg=(pIndex[k].x-start)%length;
						if(seg==0)//its at the 0 end, use the length-1 end
							pIndex[k].y=start+(length-1)+chain*length;
						else if(seg==length-1)//its at the length-1 end, use the 0 end
							pIndex[k].y=start+chain*length;
						else//its in the middle, just use the 0 end
							pIndex[k].y=start+chain*length;
					}
				}
			}
		}
	}
	
	//molecule beads
	std::vector< std::vector<int> > pCIndices;
	for(int i=0;i<mCIndices.size();i++)
	{
		std::vector<int> pCIndex;
		if(mCIndices[i]<System.readNMolecules() && mCIndices[i]>=0)
		{
			molecule<double,fourVector<int> > m=System.getMolecule()[mCIndices[i]];
			if(m.readType()==BEAD)
			{
				fourVector<int> *bond=m.getBonds();
				mCCutoffs[i]+=m.getConstants()[BEADRADIUS];
				
				for(int j=0;j<m.readNBond();j++)
					pCIndex.push_back(m.getBonds()[j].x);
			}
			else
			{
				std::cerr << "Molecule wrong type! Must be BEAD (" << BEAD << ") type! #" 
				<< mCIndices[i] << " != [0," 
				<< System.readNMolecules() << "]" << std::endl;
				return -1;
			}
		}
		else
		{
			std::cerr << "Molecule out of bounds! Must be BEAD (" << BEAD << ") type! #" 
			<< mCIndices[i] << " != [0," 
			<< System.readNMolecules() << "]" << std::endl;
			return -1;
		}
		pCIndices.push_back(pCIndex);
	}
	
	//types as beads with cutoff
	for(int i=0;i<cTypes.size();i++)
	{
		mCCutoffs.push_back(cCutoffs[i]);
		std::vector<int> pCIndex;
		for(int j=0;j<nParticles;j++)
			if(p[j].type==cTypes[i])
				pCIndex.push_back(j);
		pCIndices.push_back(pCIndex);
	}
	
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<int> sFlagOrigin;
	int currentFrame=0;
	double avgNeighbors=0;
	
	double trialSpacing=surfSpacing;
	
	while(xyzFile.load())
	{
		std::cerr << time << '\t' << currentFrame << std::endl;
		
		if(beadPotentialFile.is_open())
		{
			double sTime=-1;
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !beadPotentialFile.eof())
				beadPotentialFile >> sTime >> beadPotential;
		}
		
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
		
		/*
			Cell List Creation
		*/
		if(surfSpacing==0)
			surfSpacing=vertSpacing;
		
		threeVector<double> cellSize;
		cellSize.x=surfSpacing;
		cellSize.y=cellSize.x;
		cellSize.z=cellSize.x;
		threeHash<threeVector<double> > cellHash(size,cellSize);
		
		std::unordered_map<int, std::vector<lipid<double> > > vCells;
		//std::unordered_map<int, std::vector<int> > vCells;
		for(int i=pIndex.size()-1;i>=0;--i)
		{
			auto& pI=pIndex[i];
			lipid<double> buf;
			
			//first or last lipid chain position
			buf.pos=p[pI.x];
			
			//orientation from current to opposing element
			buf.orient.x=p[pI.x].x-p[pI.y].x;
			buf.orient.y=p[pI.x].y-p[pI.y].y;
			buf.orient.z=p[pI.x].z-p[pI.y].z;
			buf.orient.x+=(buf.orient.x<-size.x/2.0)?size.x:0;
			buf.orient.y+=(buf.orient.y<-size.y/2.0)?size.y:0;
			buf.orient.z+=(buf.orient.z<-size.z/2.0)?size.z:0;
			buf.orient.x-=(buf.orient.x>size.x/2.0)?size.x:0;
			buf.orient.y-=(buf.orient.y>size.y/2.0)?size.y:0;
			buf.orient.z-=(buf.orient.z>size.z/2.0)?size.z:0;
			buf.orient=unitVector(buf.orient);
			
			//just to keep track of it
			buf.index=i;
			
			vCells[cellHash(p[pI.x])].push_back(buf);
			//vCells[cellHash(pVal)].push_back(pVal.type);
		}
		
		
		/*
			Surface Discrimination
		*/
		
		bool thresholdExceeded=false;
		
		//std::cerr << vCells.size() << std::endl;
		
		std::vector< std::unordered_map<int,std::vector<lipid<double> > > > surfaces;
		
		std::vector<int> sFlag(pIndex.size(),-1);
		
		for(auto& cell:vCells)
		{
			std::deque< lipid<double> > stack;
			stack.push_back(cell.second[0]);
			if(sFlag[cell.second[0].index]==-1)
			{
				sFlag[cell.second[0].index]=surfaces.size();
				
				std::unordered_map<int,std::vector<lipid<double> > > surface;
				while(stack.size()>0)
				{
					//std::cerr << stack.size() << std::endl;
					lipid<double> pC=stack.back();
					stack.pop_back();
					int currentCell=cellHash(pC.pos);
					surface[currentCell].push_back(pC);
					threeVector<int> cellPos=cellHash.unHash(currentCell);
					
					for(blockKey<threeVector<int> > nCellOffset;nCellOffset!=nCellOffset.end();++nCellOffset)
					{
						threeVector<int> nCellPos=cellPos;
						
						nCellPos.x+=(*nCellOffset).x;
						nCellPos.y+=(*nCellOffset).y;
						nCellPos.z+=(*nCellOffset).z;
						
						int nCellHash=cellHash.hash(nCellPos);
						threeVector<double> minImg=cellHash.minImg(nCellPos);
						auto nCell=vCells.find(nCellHash);
						
						if(nCell!=vCells.end())
						{
							for(auto &pN:nCell->second)
							{
								//threeVector<double> pNp=p[pN];
								threeVector<double> d;
								d.x=pN.pos.x-pC.pos.x+minImg.x;
								d.y=pN.pos.y-pC.pos.y+minImg.y;
								d.z=pN.pos.z-pC.pos.z+minImg.z;
								
								double dr=d.x*d.x+d.y*d.y+d.z*d.z;
								if(sFlag[pN.index]==-1 && dr<=trialSpacing*trialSpacing)
								{
									stack.push_back(pN);
									sFlag[pN.index]=surfaces.size();
								}
							}
						}
					}
				}
				surfaces.push_back(surface);
			}
		}
		double mean=0;
		double stdDev=0;
		for(auto& surface:surfaces)
			for(auto& element:surface)
				mean+=element.second.size();
		mean/=surfaces.size();
		for(auto& surface:surfaces)
		{
			int surfaceSize=0;
			for(auto& element:surface)
				surfaceSize+=element.second.size();
			stdDev+=pow(surfaceSize-mean,2.0);
		}
		stdDev/=surfaces.size();
		stdDev=sqrt(stdDev);
		std::cerr << "Number of surfaces found: " << surfaces.size() << " at spacing " << trialSpacing;
		std::cerr << " with mean " << mean << " and standard deviation " << stdDev << std::endl;
	
		/*
			Something went wrong with surface, skipping unless resurf occurs
		*/	
		if(nSurfaces>0 && surfaces.size()!=nSurfaces)
		{
			thresholdExceeded=true;
			if(surfaces.size()>0)
				surfaces.clear();
		}
		
		/*
			Resurf
		*/
		if(currentFrame%resurf==0 && nSurfaces>0 && surfaces.size()!=nSurfaces)
		{
			thresholdExceeded=false;
			std::cerr << "Attempting to obtain " << nSurfaces << " surfaces!" << std::endl;
			
			while(surfaces.size()!=nSurfaces && !thresholdExceeded)
			{
				if(surfaces.size()>0)
					surfaces.clear();
				for(auto& flag:sFlag)
					flag=-1;
				for(auto& cell:vCells)
				{
					std::deque< lipid<double> > stack;
					stack.push_back(cell.second[0]);
					if(sFlag[cell.second[0].index]==-1)
					{
						sFlag[cell.second[0].index]=surfaces.size();
						
						std::unordered_map<int,std::vector<lipid<double> > > surface;
						while(stack.size()>0)
						{
							//std::cerr << stack.size() << std::endl;
							lipid<double> pC=stack.back();
							stack.pop_back();
							int currentCell=cellHash(pC.pos);
							surface[currentCell].push_back(pC);
							threeVector<int> cellPos=cellHash.unHash(currentCell);
							
							for(blockKey<threeVector<int> > nCellOffset;nCellOffset!=nCellOffset.end();++nCellOffset)
							{
								threeVector<int> nCellPos=cellPos;
								
								nCellPos.x+=(*nCellOffset).x;
								nCellPos.y+=(*nCellOffset).y;
								nCellPos.z+=(*nCellOffset).z;
								
								int nCellHash=cellHash.hash(nCellPos);
								threeVector<double> minImg=cellHash.minImg(nCellPos);
								auto nCell=vCells.find(nCellHash);
								
								if(nCell!=vCells.end())
								{
									for(auto &pN:nCell->second)
									{
										//threeVector<double> pNp=p[pN];
										threeVector<double> d;
										d.x=pN.pos.x-pC.pos.x+minImg.x;
										d.y=pN.pos.y-pC.pos.y+minImg.y;
										d.z=pN.pos.z-pC.pos.z+minImg.z;
										
										double dr=d.x*d.x+d.y*d.y+d.z*d.z;
										if(sFlag[pN.index]==-1 && dr<=trialSpacing*trialSpacing)
										{
											stack.push_back(pN);
											sFlag[pN.index]=surfaces.size();
										}
									}
								}
							}
						}
						surfaces.push_back(surface);
					}
				}
				
				//std::cerr << "Number of surfaces found (before pruning): " << surfaces.size() << std::endl;
				
				//for(int i=0;i<sFlag.size();i++)
				//	if(sFlag[i]<50)
				//		std::cout << sFlag[i] << ' ' << p[pIndex[i].x].x << ' ' << p[pIndex[i].x].y << ' ' << p[pIndex[i].x].z << std::endl;
				
				mean=0;
				stdDev=0;
				for(auto& surface:surfaces)
					for(auto& element:surface)
						mean+=element.second.size();
				mean/=surfaces.size();
				for(auto& surface:surfaces)
				{
					int surfaceSize=0;
					for(auto& element:surface)
						surfaceSize+=element.second.size();
					stdDev+=pow(surfaceSize-mean,2.0);
				}
				stdDev/=surfaces.size();
				stdDev=sqrt(stdDev);
				std::cerr << "Number of surfaces found: " << surfaces.size() << " at spacing " << trialSpacing;
				std::cerr << " with mean " << mean << " and standard deviation " << stdDev << std::endl;
				
				while(stdDev>mean*maxVariance && surfaces.size()>1)
				{
					mean=0;
					stdDev=0;
					auto smallestSurface=surfaces.begin();
					for(auto surface=surfaces.end()-1;surface!=surfaces.begin();--surface)
						if((*surface).size()<(*smallestSurface).size())
							smallestSurface=surface;
					surfaces.erase(smallestSurface);
					
					for(auto& surface:surfaces)
						for(auto& element:surface)
							mean+=element.second.size();
					mean/=surfaces.size();
					for(auto& surface:surfaces)
					{
						int surfaceSize=0;
						for(auto& element:surface)
							surfaceSize+=element.second.size();
						stdDev+=pow(surfaceSize-mean,2.0);
					}
					stdDev/=surfaces.size();
					stdDev=sqrt(stdDev);
					std::cerr << "Number of surfaces found: " << surfaces.size() << " at spacing " << trialSpacing;
					std::cerr << " with mean " << mean << " and standard deviation " << stdDev << std::endl;
				}
				
				for(auto& surface:surfaces)
					std::cerr << surface.size() << ' ';
				std::cerr << std::endl;
				
				if(surfaces.size()>nSurfaces)
				{
					trialSpacing+=0.01;
				}
				if(surfaces.size()<nSurfaces && surfaces.size()>0)
				{
					trialSpacing-=0.01;
				}
				if(surfaces.size()==0 || trialSpacing>surfSpacing || trialSpacing==0)
				{
					thresholdExceeded=true;
				}
				if(surfaces.size()>nSurfaces && trialSpacing<surfSpacing)
				{
					trialSpacing=surfSpacing;
				}
				//std::cin.get();
			}
		}
		
		/*
			Contact shells for faster searching
		*/
		std::vector< std::vector< threeVector<int> > > shells;
		
		for(auto& cRad:mCCutoffs)
		{
			double maxCellLength=cellHash.c.x;
			maxCellLength=cellHash.c.x;
			maxCellLength=(cellHash.c.y>maxCellLength)?cellHash.c.y:maxCellLength;
			maxCellLength=(cellHash.c.z>maxCellLength)?cellHash.c.z:maxCellLength;
			std::vector<threeVector<int> > shell;
			int shellLength=static_cast<int>(cRad+2.0)*2+4;
			for(blockKey<threeVector<int> > nCellOffset(shellLength);nCellOffset!=nCellOffset.end();++nCellOffset)
			{
				bool blockInRange=false;
				for(int i=0;i<27;i++)
				{
					threeVector<double> corner;
					corner.x=((*nCellOffset).x+i%3-1)*cellHash.c.x;
					corner.y=((*nCellOffset).y+static_cast<int>(i/3)%3-1)*cellHash.c.y;
					corner.z=((*nCellOffset).z+static_cast<int>(i/9)-1)*cellHash.c.z;
						
					double dr=corner.x*corner.x+corner.y*corner.y+corner.z*corner.z;
					if(dr>=(cRad)*(cRad) && dr<=(cRad+maxCellLength)*(cRad+maxCellLength))
						blockInRange=true;
				}
				if(blockInRange)
					shell.push_back((*nCellOffset));
			}
			shells.push_back(shell);
		}
		
		/*
			Mark contact and non-contact lipids
		*/
		std::vector<bool> cFlag(pIndex.size(),false);
		auto cShell=shells.begin();
		auto mCCutoff=mCCutoffs.begin();
		for(auto& pCIndex:pCIndices)
		{
			double cutoffSqr=(*mCCutoff);
			cutoffSqr*=cutoffSqr;
			for(auto& pCI:pCIndex)
			{
				auto& pC=p[pCI];
				int currentCell=cellHash(pC);
				threeVector<int> cellPos=cellHash.unHash(currentCell);
				
				for(auto& cShellI:(*cShell))
				{
					threeVector<int> nCellPos=cellPos;
					
					nCellPos.x+=cShellI.x;
					nCellPos.y+=cShellI.y;
					nCellPos.z+=cShellI.z;
					
					int nCellHash=cellHash.hash(nCellPos);
					threeVector<double> minImg=cellHash.minImg(nCellPos);
					auto nCell=vCells.find(nCellHash);
					
					if(nCell!=vCells.end())
					{
						for(auto &pN:nCell->second)
						{
							//threeVector<double> pNp=p[pN];
							threeVector<double> d;
							d.x=pN.pos.x-pC.x+minImg.x;
							d.y=pN.pos.y-pC.y+minImg.y;
							d.z=pN.pos.z-pC.z+minImg.z;
							
							double dr=d.x*d.x+d.y*d.y+d.z*d.z;
							if(dr<=cutoffSqr)
								cFlag[pN.index]=true;
						}
					}
				}
			}
			cShell++;
			mCCutoff++;
		}
		
		/*
			Vertex shell for local orientation
		*/
		std::vector<threeVector<int> > vertShell;
		if(vertSpacing!=surfSpacing)
		{
			int shellLength=static_cast<int>(vertSpacing)*2+4;
			for(blockKey<threeVector<int> > nCellOffset(shellLength);nCellOffset!=nCellOffset.end();++nCellOffset)
			{
				bool blockInRange=false;
				for(int i=0;i<27;i++)
				{
					threeVector<double> corner;
					corner.x=((*nCellOffset).x+i%3-1)*cellHash.c.x;
					corner.y=((*nCellOffset).y+static_cast<int>(i/3)%3-1)*cellHash.c.y;
					corner.z=((*nCellOffset).z+static_cast<int>(i/9)-1)*cellHash.c.z;
						
					double dr=corner.x*corner.x+corner.y*corner.y+corner.z*corner.z;
					if(dr<=(vertSpacing)*(vertSpacing))
						blockInRange=true;
				}
				if(blockInRange)
					vertShell.push_back((*nCellOffset));
			}
		}
		else
		{
			for(int i=0;i<27;i++)
			{
				threeVector<int> corner;
				corner.x=(i%3-1);
				corner.y=(static_cast<int>(i/3)%3-1);
				corner.z=(static_cast<int>(i/9)-1);
				vertShell.push_back(corner);
			}
		}
		
		/*
			Determine all local orientations!
		*/
		//for(int nVert=0;nVert<nVertSpacing;nVert++)
		double gradSquaredNonContactTotal=0;
		double laplacianIntNonContactTotal=0;
		double gradSquaredContactTotal=0;
		double laplacianIntContactTotal=0;
		int surfaceIndex=0;
		
		std::vector<double> zProfileContact(1,0);
		std::vector<int> zCountContact(1,0);
		std::vector<double> zProfileNonContact(1,0);
		std::vector<int> zCountNonContact(1,0);
		if(zSlice>0)
		{
			zProfileContact.resize(static_cast<int>(size.z/zSlice)+1,0.0);
			zCountContact.resize(static_cast<int>(size.z/zSlice)+1,0);
			zProfileNonContact.resize(static_cast<int>(size.z/zSlice)+1,0.0);
			zCountNonContact.resize(static_cast<int>(size.z/zSlice)+1,0);
		}
		for(auto& surface:surfaces)
		{
			
			std::vector<double> gradSquareValuesContact;
			std::vector< fourVector<double> > gradSquarePos;
			double gradSquaredContact=0;
			double laplacianIntContact=0;
			double localRmsContact=0;
			
			std::vector<double> gradSquareValuesNonContact;
			double gradSquaredNonContact=0;
			double laplacianIntNonContact=0;
			double localRmsNonContact=0;
			
			int nElements=0;
			double percentCompletion=0;
			for(auto& element:surface)
			{
				if(percentCompletion<=100.0*static_cast<double>(nElements++)/static_cast<double>(surface.size()))
				{
					//std::cerr << percentCompletion << "%" << std::endl;
					percentCompletion+=10;
				}
				
				threeVector<double> oneZ=0;
				oneZ.z=1.0;
				
				int currentHash=element.first;
				threeVector<int> iCoor=cellHash.unHash(currentHash);
				
				
				
				
				std::vector<lipid<double> > nearbyGroups;
				
				for(auto& vSCoor:vertShell)
				{
					threeVector<int> nCoor=iCoor;
					nCoor.x+=vSCoor.x;
					nCoor.y+=vSCoor.y;
					nCoor.z+=vSCoor.z;
					
					int nCellHash=cellHash.hash(nCoor);
					if(nCellHash!=currentHash)
					{
						threeVector<double> minImg=cellHash.minImg(nCoor);
						auto nCell=vCells.find(nCellHash);
						
						if(nCell!=vCells.end())
						{
							for(auto& point:nCell->second)
							{
								nearbyGroups.push_back(point);
								nearbyGroups.back().pos.x+=minImg.x;
								nearbyGroups.back().pos.y+=minImg.y;
								nearbyGroups.back().pos.z+=minImg.z;
							}
						}
					}
				}
				
				auto& aGroup=element.second;
				
				int gradSquareValuesNonContactStart=gradSquareValuesNonContact.size();
				//double gradSquaredLow=gradSquared;
				//do this for every particle in aGroup
				for(auto& aP:aGroup)
				{
					double deltaHdeltaX=1.1;
					//double patchSize=initPatchSize;
					//while(deltaHdeltaX>1.0)
					
					threeVector<double> normalSum=0;
					int nNormals=0;
					for(auto& aG:aGroup)
					{
						threeVector<double> d;
						d.x=aG.pos.x-aP.pos.x;
						d.y=aG.pos.y-aP.pos.y;
						d.z=aG.pos.z-aP.pos.z;
						//closest edge
						if(d.x*d.x+d.y*d.y+d.z*d.z<8.0 && dotProduct(aG.orient,aP.orient)>0)
							normalSum+=aG.orient;
					}
					for(auto& nG:nearbyGroups)
					{
						threeVector<double> d;
						d.x=nG.pos.x-aP.pos.x;
						d.y=nG.pos.y-aP.pos.y;
						d.z=nG.pos.z-aP.pos.z;
						//closest edge
						if(d.x*d.x+d.y*d.y+d.z*d.z<8.0 && dotProduct(nG.orient,aP.orient)>0)
							normalSum+=nG.orient;
					}
					
					//normalSum/=static_cast<double>(nNormals);
					
					//get the projection matrix, ROW_MAJOR
					threeVector<threeVector<double> > M;
					//M=I-(aP*aP^T)/aP^2
					//position<double> aPu=unitVector(aP);
					threeVector<double> aPu=unitVector(normalSum);
							
					oneZ.z=1.0;
					if(dotProduct(aPu,oneZ)<0)
						oneZ.z=-oneZ.z;
					
					threeVector<double> aPuZ=crossProduct(aPu,oneZ);
					
					double w=sqrt(1+dotProduct(aPu,oneZ));
					double b = aPuZ.x;
					double c = aPuZ.y;
					double d = aPuZ.z;
					double a = w;
					
					//stolen from boost, boosted from boost?
					double aa = a*a;
					double ab = a*b;
					double ac = a*c;
					double ad = a*d;
					double bb = b*b;
					double bc = b*c;
					double bd = b*d;
					double cc = c*c;
					double cd = c*d;
					double dd = d*d;
					
					double norme_carre = aa+bb+cc+dd;
					
					M.x.x=(aa + bb - cc - dd)/norme_carre;
					M.x.y=2 * (-ad + bc)/norme_carre;
					M.x.z=2 * (ac + bd)/norme_carre;
					M.y.x= 2 * (ad + bc)/norme_carre;
					M.y.y=(aa - bb + cc - dd)/norme_carre;
					M.y.z=2 * (-ab + cd)/norme_carre;
					M.z.x=2 * (-ac + bd)/norme_carre;
					M.z.y=2 * (ab + cd)/norme_carre;
					M.z.z=(aa - bb - cc + dd)/norme_carre;
							
					//our z oriented patch (heights)
					//number of points used in this patch
					threeVector<threeVector<double> > nZPatch;
					//our z oriented patch (heights)
					threeVector<threeVector<double> > zPatch;
						
					for(int a=0;a<3;a++)
						for(int b=0;b<3;b++)
							zPatch.s[a].s[b]=0;
					for(int a=0;a<3;a++)
						for(int b=0;b<3;b++)
							nZPatch.s[a].s[b]=0;
					
					//transform the neighbors by the projection matrix and place on grid
					for(auto &nG:nearbyGroups)
					{
						threeVector<double> sG,qG;
						
						qG.x=nG.pos.x-aP.pos.x;
						qG.y=nG.pos.y-aP.pos.y;
						qG.z=nG.pos.z-aP.pos.z;
								
						sG.x=qG.x*M.x.x+qG.y*M.x.y+qG.z*M.x.z;
						sG.y=qG.x*M.y.x+qG.y*M.y.y+qG.z*M.y.z;
						sG.z=qG.x*M.z.x+qG.y*M.z.y+qG.z*M.z.z;
								
						sG.z*=oneZ.z;
								
						threeVector<int> patchCoor=-1.0;
						if(sG.x>-1.5*patchSize && sG.x<1.5*patchSize)
							patchCoor.x=static_cast<int>((sG.x+patchSize*1.5)/(patchSize));
						if(sG.y>-1.5*patchSize && sG.y<1.5*patchSize)
							patchCoor.y=static_cast<int>((sG.y+patchSize*1.5)/(patchSize));
						
						//farthest corner
						//if(sG.x*sG.x+sG.y*sG.y+sG.z*sG.z<4.5*patchSize*patchSize)
						if(sG.z<patchSize && sG.z>-patchSize)
						{
							if(dotProduct(nG.orient,aPu)>0)
							{
								if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
								{
									zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
									nZPatch.s[patchCoor.x].s[patchCoor.y]++;
								}
							}
						}
					}
							
					for(auto &nG:aGroup)
					{
						threeVector<double> sG,qG;
						
						qG.x=nG.pos.x-aP.pos.x;
						qG.y=nG.pos.y-aP.pos.y;
						qG.z=nG.pos.z-aP.pos.z;
								
						sG.x=qG.x*M.x.x+qG.y*M.x.y+qG.z*M.x.z;
						sG.y=qG.x*M.y.x+qG.y*M.y.y+qG.z*M.y.z;
						sG.z=qG.x*M.z.x+qG.y*M.z.y+qG.z*M.z.z;
						
						sG.z*=oneZ.z;
						//sG=qG;
						
						threeVector<int> patchCoor=-1.0;
						if(sG.x>-1.5*patchSize && sG.x<1.5*patchSize)
							patchCoor.x=static_cast<int>((sG.x+patchSize*1.5)/(patchSize));
						if(sG.y>-1.5*patchSize && sG.y<1.5*patchSize)
							patchCoor.y=static_cast<int>((sG.y+patchSize*1.5)/(patchSize));
						//if(sG.z>0 && sG.z<1.5*vertSpacing)
						//	patchCoor.z=static_cast<int>(sG.z/vertSpacing);
						//farthest corner
						//if(sG.x*sG.x+sG.y*sG.y+sG.z*sG.z<4.5*patchSize*patchSize)// && dotProduct(*nearbyNormal,aPu)>0)
						if(sG.z<patchSize && sG.z>-patchSize)
						{
							
							if(dotProduct(nG.orient,aPu)>0)
							{
								
								if(patchCoor.x>-1 && patchCoor.y>-1 && patchCoor.x<3 && patchCoor.y<3)
								{
									zPatch.s[patchCoor.x].s[patchCoor.y]+=sG.z;
									nZPatch.s[patchCoor.x].s[patchCoor.y]++;
								}
							}
						}
					}
					//throw 0;
					bool badPatch=false;
					int nL=0;
					for(int a=0;a<3;a++)
					{
						for(int b=0;b<3;b++)
						{
							nL+=nZPatch.s[a].s[b];
							//if(nZPatch.s[a].s[b]==0 && a!=b && !(a==0 && b==2) && !(a==2 && b==0))
							//{
							//	nBadPatches++;
								//badPatch=true;
							//}
						}
					}
					nL+=1;
							
					//average of every patch
					if(nZPatch.s[0].s[0]!=0)
						zPatch.s[0].s[0]/=nZPatch.s[0].s[0];
					if(nZPatch.s[0].s[1]!=0)
						zPatch.s[0].s[1]/=nZPatch.s[0].s[1];
					if(nZPatch.s[0].s[2]!=0)
						zPatch.s[0].s[2]/=nZPatch.s[0].s[2];
					if(nZPatch.s[1].s[0]!=0)
						zPatch.s[1].s[0]/=nZPatch.s[1].s[0];
					if(nZPatch.s[1].s[1]!=0)
						zPatch.s[1].s[1]/=(nZPatch.s[1].s[1]+1.0);//aP is here, but contributes 0
					if(nZPatch.s[1].s[2]!=0)
						zPatch.s[1].s[2]/=nZPatch.s[1].s[2];
					if(nZPatch.s[2].s[0]!=0)
						zPatch.s[2].s[0]/=nZPatch.s[2].s[0];
					if(nZPatch.s[2].s[1]!=0)
						zPatch.s[2].s[1]/=nZPatch.s[2].s[1];
					if(nZPatch.s[2].s[2]!=0)
						zPatch.s[2].s[2]/=nZPatch.s[2].s[2];
					
					zPatch.s[0].s[1]-=zPatch.s[1].s[1];
					zPatch.s[1].s[0]-=zPatch.s[1].s[1];
					zPatch.s[1].s[2]-=zPatch.s[1].s[1];
					zPatch.s[2].s[1]-=zPatch.s[1].s[1];
					
					zPatch.s[0].s[0]-=zPatch.s[1].s[1];
					zPatch.s[0].s[2]-=zPatch.s[1].s[1];
					zPatch.s[2].s[0]-=zPatch.s[1].s[1];
					zPatch.s[2].s[2]-=zPatch.s[1].s[1];
					
					
					zPatch.s[1].s[1]=0.0;
							
					
					
					if(!badPatch)
					{
						//double myGradSquared=0.66689*pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
						//		zPatch.s[2].s[1]+zPatch.s[0].s[1])/2.0+
						//		(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
						//		zPatch.s[2].s[2]+zPatch.s[0].s[0])/4.0-
						//		3.0*zPatch.s[1].s[1])/(patchSize*patchSize),2.0);//average 0.66689 nm^2 per lipid
						
						double patchArea=patchSize*patchSize;
						
						twoVector<double> dZ;
						//dZ.x=(zPatch.s[0].s[1]+zPatch.s[2].s[1])/2.0;
						//dZ.y=(zPatch.s[1].s[0]+zPatch.s[1].s[2])/2.0;
						dZ.x=(zPatch.s[0].s[1]-zPatch.s[1].s[1]);
						dZ.y=(zPatch.s[1].s[0]-zPatch.s[1].s[1]);
						
						double hx=(zPatch.s[2].s[1]-zPatch.s[0].s[1])/(patchSize*2.0);
						double hy=(zPatch.s[1].s[2]-zPatch.s[1].s[0])/(patchSize*2.0);
						
						double hyy=(zPatch.s[1].s[2]+zPatch.s[1].s[0]-zPatch.s[1].s[1]*2.0)/patchArea;
						double hxx=(zPatch.s[2].s[1]+zPatch.s[0].s[1]-zPatch.s[1].s[1]*2.0)/patchArea;
						double hxy=(zPatch.s[2].s[2]+zPatch.s[0].s[0]-zPatch.s[2].s[0]-zPatch.s[0].s[2])/(4.0*patchArea);
							
							
						double areaPerLipid=patchArea/(nZPatch.s[1].s[1]+1.0);
						//double area=sqrt(1.0+(zPatch.s[0].s[1]*zPatch.s[0].s[1]+zPatch.s[1].s[0]*zPatch.s[1].s[0])/areaPerLipid);
						double areaScaler=1.0+(dZ.x*dZ.x+dZ.y*dZ.y)/patchArea;
								
						double myLaplacian=((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
						zPatch.s[2].s[1]+zPatch.s[0].s[1])*2.0/3.0+
						(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
						zPatch.s[2].s[2]+zPatch.s[0].s[0])/6.0-
						(10.0/3.0)*zPatch.s[1].s[1])/(patchArea);
						
						//if(myLaplacian>0)
						//	std::cout << 1 << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << '\n';
						//else
						//	std::cout << 0 << ' ' << aP.x << ' ' << aP.y << ' ' << aP.z << '\n';
						
						//\sum_{j}^{N_p} \sum_{k}^{N_s} \left( (\bigtriangledown^2 h_k)^2 \bigtriangleup a_k \right)
						//double myGradSquared=pow(myLaplacian,2.0)*areaPerLipid;
								
						//double myGradSquared=(1.0+hyy)*hxx+(1.0+hxx)*hyy+2.0*hx*hy*hxy;
						//myGradSquared*=myGradSquared*areaPerLipid;
						//myGradSquared/=pow(1.0+(hx*hx+hy*hy),2.5);
						//Correction in August2022
						double myGradSquared=(1.0+hy*hy)*hxx+(1.0+hx*hx)*hyy-2.0*hx*hy*hxy;
						myGradSquared*=myGradSquared*areaPerLipid;
						myGradSquared/=pow(1.0+(hx*hx+hy*hy),2.5);
						//double myGradSquared=pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
						//zPatch.s[2].s[1]+zPatch.s[0].s[1])*2.0/3.0+
						//(zPatch.s[0].s[2]+zPatch.s[2].s[0]+
						//zPatch.s[2].s[2]+zPatch.s[0].s[0])/6.0-
						//(10.0/3.0)*zPatch.s[1].s[1])/(patchArea),2.0);//average 0.66689 nm^2 per lipid
						//double myGradSquared=pow(((zPatch.s[1].s[2]+zPatch.s[1].s[0]+
						//zPatch.s[2].s[1]+zPatch.s[0].s[1])-
						//4.0*zPatch.s[1].s[1])/(patchArea),2.0);//average 0.66689 nm^2 per lipid
						
						if(cFlag[aP.index])
						{
							gradSquaredContact+=myGradSquared;
							laplacianIntContact+=myLaplacian;
							gradSquareValuesContact.push_back(myGradSquared);
							if(zSlice>0)
							{
								zProfileContact[static_cast<int>(aP.pos.z/zSlice)]+=myGradSquared;
								zCountContact[static_cast<int>(aP.pos.z/zSlice)]++;
							}
						} else {
							gradSquaredNonContact+=myGradSquared;
							laplacianIntNonContact+=myLaplacian;
							gradSquareValuesNonContact.push_back(myGradSquared);
							if(zSlice>0)
							{
								zProfileNonContact[static_cast<int>(aP.pos.z/zSlice)]+=myGradSquared;
								zCountNonContact[static_cast<int>(aP.pos.z/zSlice)]++;
							}
						}
						
						if(nVertSpacing<0)
						{
							fourVector<double> aV;
							aV.x=aP.pos.x;
							aV.y=aP.pos.y;
							aV.z=aP.pos.z;
							aV.t=myGradSquared;
							gradSquarePos.push_back(aV);
							
						}
					}
				}
				
		
			}
			//std::cerr << percentCompletion << "%" << std::endl;
			if(nVertSpacing<0)
			{
				/*
				std::string buf;
				buf="heatMap_SN";
				std::string whatever;
				std::stringstream parseNumber;
				parseNumber << s << "_" << patchSize;
				parseNumber >> whatever;
				
				buf+=whatever;
				buf+="_";
				//buf="patchMin_";
				buf+=name;
				buf+=".xyz";
				std::fstream dataFile;
				
				dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
				//auto gradSquareValue=gradSquareValuesNonContact.begin();
				
				double avgGradSquare=gradSquared/static_cast<double>(gradSquareValuesNonContact.size());
				double minGrad=avgGradSquare;
				double maxGrad=avgGradSquare;
				for(auto &value:gradSquareValuesNonContact)
				{
					if(value<minGrad)
						minGrad=value;
					if(value>maxGrad)
						maxGrad=value;
				}
				double dGrad=(maxGrad-minGrad)/20.0;
				
				
				
				dataFile << gradSquareValuesNonContact.size() << "\nalksdfj\n";
				
				auto cPos=gradSquarePos.begin();
				for(auto &value:gradSquareValuesNonContact)
				{
					dataFile << floor((value-minGrad)/dGrad) << ' ' << cPos->x << ' ' << cPos->y << ' ' << cPos->z << '\n';
					cPos++;
				}
				dataFile << std::endl;
				dataFile.close();
				
				
				buf.clear();
				buf="gradSorted_SN";
				
				buf+=whatever;
				buf+="_";
				//buf="patchMin_";
				buf+=name;
				buf+=".dat";
				
				dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
				
				std::sort(gradSquareValuesNonContact.begin(), gradSquareValuesNonContact.end());
				for(auto &value:gradSquareValuesNonContact)
				{
					dataFile << value << '\n';
				}
				dataFile << std::endl;
				dataFile.close();
				*/
			}
			
			if(surfaces.size()>0)
			{
				if(nVertSpacing>1)
				{
					/*
				std::string buf;
					buf="patchProfile_SN";
						std::string whatever;
					std::stringstream parseNumber;
					parseNumber << s;
					parseNumber >> whatever;
					
					buf+=whatever;
					buf+="_";
					//buf="patchMin_";
					buf+=name;
				buf+=".dat";
						std::fstream dataFile;
				
					dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
					//auto gradSquareValue=gradSquareValuesNonContact.begin();
						
				
						double avgGradSquare=gradSquared/static_cast<double>(gradSquareValuesNonContact.size());
				double rms=0;
						for(auto &value:gradSquareValuesNonContact)
						rms+=pow(value-avgGradSquare,2.0);
					rms=sqrt(rms/static_cast<double>(gradSquareValuesNonContact.size()));
				localRms/=static_cast<double>(gradSquareValuesNonContact.size());
						if(gradSquareValuesNonContact.size()!=0)
					{
					dataFile << patchSize << ' ' << gradSquared << ' ' << rms << ' '; 
							dataFile << localRms << ' ' << gradSquareValuesNonContact.size() << ' ' << laplacianInt << std::endl;
					}
					
					dataFile.close();
					*/
				}
				else if (nVertSpacing==1 || nVertSpacing==0)
				{
					/*
					std::string buf;
					buf="gradSquare_SNC";
					std::string whatever;
					std::stringstream parseNumber;
					parseNumber << surfaceIndex << "_" << patchSize;
					parseNumber >> whatever;
					
					buf+=whatever;
					buf+="_";
					buf+=name;
					buf+=".dat";
					std::fstream dataFile;
					
					dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
					//auto gradSquareValue=gradSquareValuesNonContact.begin();
					
					//int nNonContact=0;
					//for(int i=0;i<cFlag.size();i++)
					//	if(!cFlag[i])
					//		nNonContact++;
					//int nContact=cFlag.size()-nNonContact;
					*/
					if(gradSquareValuesNonContact.size()!=0)
					{
						//dataFile << time << ' ' << gradSquaredNonContact << ' ' << gradSquareValuesNonContact.size() << ' ' << laplacianIntNonContact << std::endl;
						gradSquaredNonContactTotal+=gradSquaredNonContact;
						laplacianIntNonContactTotal+=laplacianIntNonContact;
					}
					/*
					dataFile.close();
					
					buf.clear();
					buf="gradSquare_SC";
					whatever.clear();
					std::stringstream parseNumber2;
					parseNumber2 << surfaceIndex << "_" << patchSize;
//					parseNumber2 << patchSize;
					parseNumber2 >> whatever;
					
					buf+=whatever;
					buf+="_";
					//buf="patchMin_";
					buf+=name;
					buf+=".dat";
					
					dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
					//auto gradSquareValue=gradSquareValuesNonContact.begin();
					*/
					if(gradSquareValuesContact.size()!=0)
					{
						//dataFile << time << ' ' << gradSquaredContact << ' ' << gradSquareValuesContact.size() << ' ' << laplacianIntContact << ' ' << beadPotential << std::endl;
						gradSquaredContactTotal+=gradSquaredContact;
						laplacianIntContactTotal+=laplacianIntContact;
					}
					/*
					dataFile.close();
					*/
				}
			}
			surfaceIndex++;
		}
		int nNonContact=0;
		for(int i=0;i<cFlag.size();i++)
			if(!cFlag[i])
				nNonContact++;
		int nContact=cFlag.size()-nNonContact;
		if(zSlice>0)
		{
			
			std::stringstream zCurvaturesContactFileName;
			zCurvaturesContactFileName << "zCurvatureContact_" << name << ".dat";
			
			std::fstream zCurvaturesContactFile(zCurvaturesContactFileName.str(),std::ios::out|std::ios::app);
			for(int i=0;i<zProfileContact.size();i++)
			{
				if(zCountContact[i]>0)
					zCurvaturesContactFile << zSlice*(static_cast<double>(i)+0.5) << '\t' << (zProfileContact[i]-0.004806*zCountContact[i]) << '\t' << zCountContact[i] << '\n';
				else
					zCurvaturesContactFile << zSlice*(static_cast<double>(i)+0.5) << '\t' << 0 << '\t' << 0 << '\n';
				//zCurvaturesContactFile << zSlice*(static_cast<double>(i)+0.5) << '\t' << zProfileContact[i] << '\t' << zCountContact[i] << '\n';
			}
			zCurvaturesContactFile << std::endl;
			zCurvaturesContactFile.close();
			
			std::stringstream zCurvaturesNonContactFileName;
			zCurvaturesNonContactFileName << "zCurvatureNonContact_" << name << ".dat";
			
			double avgGSNCperLipid=gradSquaredNonContactTotal/static_cast<double>(nNonContact);
			std::fstream zCurvaturesNonContactFile(zCurvaturesNonContactFileName.str(),std::ios::out|std::ios::app);
			for(int i=0;i<zProfileNonContact.size();i++)
			{
				if(zCountNonContact[i]>0)
					zCurvaturesNonContactFile << zSlice*(static_cast<double>(i)+0.5) << '\t' << (zProfileNonContact[i]-0.004806*zCountNonContact[i]) << '\t' << zCountNonContact[i] << '\n';
				else
					zCurvaturesNonContactFile << zSlice*(static_cast<double>(i)+0.5) << '\t' << 0 << '\t' << 0 << '\n';
				//zCurvaturesNonContactFile << zSlice*(static_cast<double>(i)+0.5) << '\t' << zProfileNonContact[i] << '\t' << zCountNonContact[i] << '\n';
			}
			zCurvaturesNonContactFile << std::endl;
			zCurvaturesNonContactFile.close();
		}
		if(gradSquaredNonContactTotal>0 || gradSquaredContactTotal>0)
		{
			std::string buf;
			buf="gradSquare_SNC";
			std::string whatever;
			std::stringstream parseNumber;
			parseNumber << "total_" << patchSize;
			parseNumber >> whatever;
			
			buf+=whatever;
			buf+="_";
			//buf="patchMin_";
			buf+=name;
			buf+=".dat";
			std::fstream dataFile;
			
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			//auto gradSquareValue=gradSquareValuesNonContact.begin();
			
			if(nNonContact!=0)
				dataFile << time << ' ' << gradSquaredNonContactTotal << ' ' << nNonContact << ' ' << laplacianIntNonContactTotal << std::endl;
			dataFile.close();
			
			buf.clear();
			buf="gradSquare_SC";
			whatever.clear();
			std::stringstream parseNumber2;
			parseNumber2 << "total_" << patchSize;
			parseNumber2 >> whatever;
			
			buf+=whatever;
			buf+="_";
			//buf="patchMin_";
			buf+=name;
			buf+=".dat";
			
			dataFile.open(buf.c_str(), std::ios::app | std::ios::out);
			//auto gradSquareValue=gradSquareValuesNonContact.begin();
			
			if(nContact!=0)
			{
				dataFile << time << ' ' << gradSquaredContactTotal << ' ' << nContact << ' ' << laplacianIntContactTotal << ' ' << beadPotential << std::endl;
			}
			dataFile.close();
		}
		//}
		currentFrame++;
		if(currentFrame>frame && frame>0)
			break;
		time+=System.readStoreInterval();
	}
	return 0;
}
