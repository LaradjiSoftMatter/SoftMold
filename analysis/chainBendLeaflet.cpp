//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

bool compareAngles(fourVector<double> a, fourVector<double> b)
{
	return a.t>b.t;
}

class bendRange {
	public:
		bendRange(std::vector<position<double> > positions, std::vector<double> bendMap, double cutoff):
			p(positions),bM(bendMap),cutoff(cutoff)
		{
			cutoffSqr=cutoff*cutoff;
			for(int i=0;i<bM.size();i++)
			{
				avgBend.push_back(bM[i]);
				nNeighbors.push_back(1);
			}
		};
		
		~bendRange()
		{};
		
		const std::vector<double> getBendAverages()
		{
			return avgBend;
		};
		
		void finalize()
		{
			for(int i=0;i<avgBend.size();i++)
				avgBend[i]/=static_cast<double>(nNeighbors[i]);
		};
		
		void operator () (int &i, int &j, threeVector<double> minImg)
		{
			threeVector<double> d;
			d.x=p[i].x-p[j].x-minImg.x;
			d.y=p[i].y-p[j].y-minImg.y;
			d.z=p[i].z-p[j].z-minImg.z;
			if(d.x>cutoff*3.0 || d.y > cutoff*3.0 || d.z>cutoff*3.0)
			{
				std::cerr << d.x << ' ' << d.y << ' ' << d.z << std::endl;
				std::cerr << p[i].x << ' ' << p[i].y << ' ' << p[i].z << std::endl;
				std::cerr << p[j].x << ' ' << p[j].y << ' ' << p[j].z << std::endl;
				std::cerr << minImg.x << ' ' << minImg.y << ' ' << minImg.z << std::endl;
				//std::cin.get();
			}
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				nNeighbors[i]++;
				nNeighbors[j]++;
				avgBend[i]+=bM[j];
				avgBend[j]+=bM[i];
			}	
		};
		
		void operator () (int &i, int &j)
		{
			threeVector<double> d;
			d.x=p[i].x-p[j].x;
			d.y=p[i].y-p[j].y;
			d.z=p[i].z-p[j].z;
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				nNeighbors[i]++;
				nNeighbors[j]++;
				avgBend[i]+=bM[j];
				avgBend[j]+=bM[i];
			}	
		};
	private:
		std::vector<double> bM;
		std::vector<position<double> > p;
		std::vector<double> avgBend;
		std::vector<int> nNeighbors;
		double cutoff;
		double cutoffSqr;
};

class orientRange {
	public:
		orientRange(std::vector<position<double> > positions, std::vector<threeVector<double> > orientMap, double cutoff):
			p(positions),aHat(orientMap),cutoff(cutoff)
		{
			cutoffSqr=cutoff*cutoff;
			for(int i=0;i<orientMap.size();i++)
			{
				nHat.push_back(orientMap[i]);
				//nNeighbors.push_back(1);
			}
		};
		
		~orientRange()
		{};
		
		const std::vector<threeVector<double> > getNHat()
		{
			return nHat;
		};
		
		void finalize()
		{
			for(auto& nH:nHat)
				nH=unitVector(nH);
		};
		
		void operator () (int &i, int &j, threeVector<double> minImg)
		{
			threeVector<double> d;
			d.x=p[i].x-p[j].x-minImg.x;
			d.y=p[i].y-p[j].y-minImg.y;
			d.z=p[i].z-p[j].z-minImg.z;
			if(d.x>cutoff*3.0 || d.y > cutoff*3.0 || d.z>cutoff*3.0 || d.x<-cutoff*3.0 || d.y<-cutoff*3.0 || d.z<-cutoff*3.0)
			{
				std::cerr << cutoff << std::endl;
				std::cerr << d.x << ' ' << d.y << ' ' << d.z << std::endl;
				std::cerr << i << ": " << p[i].x << ' ' << p[i].y << ' ' << p[i].z << std::endl;
				std::cerr << j << ": " << p[j].x << ' ' << p[j].y << ' ' << p[j].z << std::endl;
				std::cerr << minImg.x << ' ' << minImg.y << ' ' << minImg.z << std::endl;
				std::cin.get();
			}
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				//nNeighbors[i]++;
				//nNeighbors[j]++;
				nHat[i]+=aHat[j];
				nHat[j]+=aHat[i];
			}	
		};
		
		void operator () (int &i, int &j)
		{
			threeVector<double> d;
			d.x=p[i].x-p[j].x;
			d.y=p[i].y-p[j].y;
			d.z=p[i].z-p[j].z;
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				//nNeighbors[i]++;
				//nNeighbors[j]++;
				nHat[i]+=aHat[j];
				nHat[j]+=aHat[i];
			}	
		};
	private:
		std::vector<threeVector<double> > aHat;
		std::vector<position<double> > p;
		std::vector<threeVector<double> > nHat;
		//std::vector<int> nNeighbors;
		double cutoff;
		double cutoffSqr;
};

template <typename T, typename V>
class sortByZ 
{
	public:
		sortByZ(T &map, V &zHat, V &center, V &size):m(map),z(zHat),c(center),s(size){}
		bool operator () (const int &i, const int &j)
		{
			double mzi=m[i].z-c.z;
			mzi-=mzi>s.z/2.0?s.z:0;
			mzi+=mzi<-s.z/2.0?s.z:0;
			double mzj=m[j].z-c.z;
			mzj-=mzj>s.z/2.0?s.z:0;
			mzj+=mzj<-s.z/2.0?s.z:0;
			if(z.z>0)
				return mzi>mzj;
			else
				return mzi<mzj;
				
		}
	private:
		const T &m;
		const V &z;
		const V &c;
		const V &s;
};

void showSwitches()
{
	std::cerr << "--name [name]\n\tLoad the [name].mpd simulation." << std::endl;
	std::cerr << "\tOutputs all molecules to terminal. Required for all other options!" << std::endl;
	
	std::cerr << "--molIndex [integer]\n\tLoad a particular CHAIN type molecule." << std::endl;
	std::cerr << "\tOutputs chainOrder_[name].dat for CHAIN type molecule." << std::endl;
	
	std::cerr << "--vesicle\n\tAssume membrane is a vesicle (spherical shape) using center of mass for radial distance." << std::endl;
	std::cerr << "\tOutputs proximal*_[name].dat for CHAIN type molecule." << std::endl;
	std::cerr << "\tOutputs distal*_[name].dat for CHAIN type molecule." << std::endl;
	
	std::cerr << "--NPmolIndex [integer]\n\tLoad a particular BEAD type molecule." << std::endl;
	std::cerr << "\tOutputs chainOrderDist_[name].dat. Requires --molIndex!" << std::endl;
	
	std::cerr << "--NPcutoff [float] [float]\n\tUse this cutoff for determining wrapping angle and orientation, range is (0, inf)." << std::endl;
	std::cerr << "\tLimits all other output to wrapped region." << std::endl;
	std::cerr << "\tOutputs wrapCosAngle_[name].dat. Requires --NPmolIndex!" << std::endl;
	std::cerr << "\tOutputs proximalOrder_[name].dat. Requires --NPmolIndex!" << std::endl;
	std::cerr << "\tOutputs distalOrder_[name].dat. Requires --NPmolIndex!" << std::endl;
	
	std::cerr << "--shellThickness [float]\n\tUse this shell thickness for a profile, range is unlimited." << std::endl;
	std::cerr << "\tOutputs chainOrderShell_[name].dat of chain order parameter and " <<
			"chainProfile_[name].dat from NP center of mass. Requires --NPmolIndex!" << 
			"proximalShell_[chainIndex]_[name].dat from NP center of mass. Requires --type!" << 
			"distalShell_[chainIndex]_[name].dat from NP center of mass. Requires --type!" << std::endl;
	
	std::cerr << "--avgShellThickness [float]\n\tUse this time for a profile, range is unlimited." << std::endl;
	std::cerr << "\tUses [float] to generate avgChainProfile_[float]_[name].dat from NP center of mass. Requires --shellThickness!" << std::endl;
	std::cerr << "\tUses [float] to generate avgProximalShell_[chainIndex]_[float]_[name].dat from NP center of mass. Requires --shellThickness!" << std::endl;
	std::cerr << "\tUses [float] to generate avgDistalShell_[chainIndex]_[float]_[name].dat from NP center of mass. Requires --shellThickness!" << std::endl;
	
	std::cerr << "--bendAngleSteps [integer]\n\tUse this stepping for a distribution of order, range is in (0, inf)." << std::endl;
	std::cerr << "\tOutputs chainOrderHist_[name].dat of chain order parameter. Requires --molIndex!" << std::endl;
	
	std::cerr << "--bendAngleStepSize [float]\n\tUse this stepping for a profile, range is in (0, 2]." << std::endl;
	std::cerr << "\tOutputs chainOrderHist_[name].dat of chain order parameter. Requires --molIndex!" << std::endl;
	
	std::cerr << "--type [float]\n\tUse this type for orientation and position of molIndex." << std::endl;
	std::cerr << "\tJust a flag. Requires --molIndex!" << std::endl;
	
	std::cerr << "--wrapDistance [float]\n\tUses [float] for size of annulus of neck." << std::endl;
	std::cerr << "\tOutputs wrapAngle_[name].dat. Requires --NPcutoff!" << std::endl;
	
	std::cerr << "--wrapOrientAngle [float]\n" << std::endl;
	std::cerr << "\tRestricts output to wrapOrientAngle in radians, range is in [0,pi/2.0]. Requires --NPcutoff!" << std::endl;
	
	std::cerr << "--map\n\tOutput a map of bending angles." << std::endl;
	std::cerr << "\tOutputs chainOrder_[name].xyz. Requires --bendAngleSteps or --bendAngleStepSize!" << std::endl;
	
	std::cerr << "--smoothMap [float]\n\tOutput a smoothed map of bending angles using a cutoff for average." << std::endl;
	std::cerr << "\tOutputs chainOrder_[name].xyz. Requires --bendAngleSteps or --bendAngleStepSize!" << std::endl;
	
	std::cerr << "--frequency [float]\n\tAdjusts outputs of map and smoothed map to occur every [float] tau." << std::endl;
	
	std::cerr << "--normalCutoff [float]\n\tCutoff for calculating normal vector in local order parameter." << std::endl;
	std::cerr << "\tOutputs localOrder_[name]_[float].dat." << std::endl;
	std::cerr << "\tOutputs localProximalOrder_[name]_[float].dat if --NPcutoff is used." << std::endl;
	std::cerr << "\tOutputs localDistalOrder_[name]_[float].dat if --NPcutoff is used." << std::endl;
	
	std::cerr << "--frames\n\tRun program over all frames in a system. Assumes frames_[name].xyz exists!" << std::endl;
	std::cerr << "--potentialWrap\n\tCalculate potential with wrapping angle. Assumes frames_[name].xyz exists!" << std::endl;
	
	std::cerr << "--help\n\tDisplay this list of options." << std::endl;
}

int main(int argc, char* argv[])
{
	if(argc==1)
	{
		showSwitches();
		return 0;
	}
	
	//variables
	std::string name;
	std::vector<int> molIndex, NPmolIndex, type;
	double NPcutoffProximal=0, NPcutoffDistal=0, shellThickness=0, bendAngleStepSize=0;
	double smoothCutoff=0, wrapDistance=0, normalCutoff=0, avgShellThickness=0;
	double wrapOrientAngle=M_PI/2.0;
	double frequency=100000;
	int bendAngleSteps=0;
	bool map=false, frames=false, vesicle=false, potentialWrap=false;
	
	//read in options
	std::stringstream cmdArg;
	for(int i=1;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	enum SWITCHESS {NAME_S, MOLINDEX_S, VESICLE_S, NPMOLINDEX_S,
			NPCUTOFF_S, SHELLTHICKNESS_S, BENDANGLESTEPS_S,
			BENDANGLESTEPSIZE_S, TYPE_S, WRAPDISTANCE_S,
			WRAPORIENTANGLE_S,MAP_S, SMOOTHMAP_S, NORMALCUTOFF_S,
			FRAMES_S, AVGSHELLTHICKNESS_S, POTENTIALWRAP_S, FREQUENCY_S, HELP_S};
	std::map<std::string, SWITCHESS> switches;
	switches[std::string("--name")]=NAME_S;
	switches[std::string("--molIndex")]=MOLINDEX_S;
	switches[std::string("--vesicle")]=VESICLE_S;
	switches[std::string("--NPmolIndex")]=NPMOLINDEX_S;
	switches[std::string("--NPcutoff")]=NPCUTOFF_S;
	switches[std::string("--shellThickness")]=SHELLTHICKNESS_S;
	switches[std::string("--bendAngleSteps")]=BENDANGLESTEPS_S;
	switches[std::string("--bendAngleStepSize")]=BENDANGLESTEPSIZE_S;
	switches[std::string("--type")]=TYPE_S;
	switches[std::string("--wrapDistance")]=WRAPDISTANCE_S;
	switches[std::string("--wrapOrientAngle")]=WRAPORIENTANGLE_S;
	switches[std::string("--map")]=MAP_S;
	switches[std::string("--smoothMap")]=SMOOTHMAP_S;
	switches[std::string("--normalCutoff")]=NORMALCUTOFF_S;
	switches[std::string("--frames")]=FRAMES_S;
	switches[std::string("--avgShellThickness")]=AVGSHELLTHICKNESS_S;
	switches[std::string("--potentialWrap")]=POTENTIALWRAP_S;
	switches[std::string("--frequency")]=FREQUENCY_S;
	switches[std::string("--help")]=HELP_S;
	
	std::string option;
	while(cmdArg >> option)
	{
		//std::cerr << cmdArg.str() << std::endl;
		//cmdArg >> option;
		std::cerr << option << std::endl;
		if(switches.end()==switches.find(option))
		{
			std::cerr << option << ": Unknown switch!" << std::endl;
			showSwitches();
			return -1;
		}
		switch(switches[option])
		{
			case NAME_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a name!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> name;
				//std::cerr << "Found: " << name << std::endl;
				break;
			case MOLINDEX_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs an integer!" << std::endl;
					showSwitches();
					return -1;
				}
				int mol;
				cmdArg >> mol;
				molIndex.push_back(mol);
				break;
			case VESICLE_S:
				vesicle=true;
				break;
			case NPMOLINDEX_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs an integer!" << std::endl;
					showSwitches();
					return -1;
				}
				int npmol;
				cmdArg >> npmol;
				NPmolIndex.push_back(npmol);
				break;
			case NPCUTOFF_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> NPcutoffProximal >> NPcutoffDistal;
				break;
			case SHELLTHICKNESS_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> shellThickness;
				break;
			case BENDANGLESTEPS_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs an integer value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> bendAngleSteps;
				break;
			case BENDANGLESTEPSIZE_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> bendAngleStepSize;
				break;
			case TYPE_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs an integer value!" << std::endl;
					showSwitches();
					return -1;
				}
				int typeValue;
				cmdArg >> typeValue;
				type.push_back(typeValue);
				break;
			case WRAPDISTANCE_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> wrapDistance;
				break;
			case WRAPORIENTANGLE_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> wrapOrientAngle;
				break;
			case MAP_S:
				map=true;
				break;
			case SMOOTHMAP_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				map=true;
				cmdArg >> smoothCutoff;
				break;
			case NORMALCUTOFF_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> normalCutoff;
				break;
			case FREQUENCY_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				cmdArg >> frequency;
				break;
			case FRAMES_S:
				frames=true;
				break;
			case AVGSHELLTHICKNESS_S:
				if(cmdArg.eof())
				{
					std::cerr << option << " needs a floating value!" << std::endl;
					showSwitches();
					return -1;
				}
				std::cerr << cmdArg.eof() << cmdArg.bad() << cmdArg.fail() << std::endl;
				cmdArg >> avgShellThickness;
				std::cerr << avgShellThickness << std::endl;
				std::cerr << cmdArg.eof() << cmdArg.bad() << cmdArg.fail() << std::endl;
				break;
			case POTENTIALWRAP_S:
				potentialWrap=true;
				break;
			case HELP_S:
				showSwitches();
				return 0;
		}
		if(cmdArg.fail())
		{
			std::cerr << "At " << option << ": Bad input!" << std::endl;
			showSwitches();
			return -1;
		}
	}
	
	if(wrapOrientAngle>M_PI/2.0 || wrapOrientAngle<0)
	{
		std::cerr << "wrapOrientAngle should be between [0,pi/2.0]!" << std::endl;
		showSwitches();
		return -1;
	}
	double wrapOrientCosTheta=cos(wrapOrientAngle);
	
	if(name.length()==0)
	{
		std::cerr << "No system name provided!" << std::endl;
		showSwitches();
		return -1;
	}
	
	if(bendAngleSteps>0 && bendAngleStepSize>0)
	{
		std::cerr << "Cannot use --bendAngleStepSize and --bendAngleSteps at the same time!";
		showSwitches();
		return -1;
	}
	
//	if(map && smoothMap>0)
//	{
//		std::cerr << "Cannot use --map and --smoothMap at the same time!";
//		showSwitches();
//		return -1;
//	}
	
	if(molIndex.size()==0 && NPmolIndex.size()>0)
	{
		std::cerr << "--NPmolIndex requires --molIndex!";
		showSwitches();
		return -1;
	}
	
	//if(shellThickness>0 && NPmolIndex.size()==0)
	//{
	//	std::cerr << "--shellThickness requires --NPmolIndex!";
	//	showSwitches();
	//	return -1;
	//}
	
	if(molIndex.size()==0 && type.size()>0)
	{
		std::cerr << "--type requires --molIndex!";
		showSwitches();
		return -1;
	}
	
	if((NPcutoffProximal>0 || NPcutoffDistal>0) && NPmolIndex.size()==0)
	{
		std::cerr << "--NPcutoff requires --NPmolIndex!";
		showSwitches();
		return -1;
	}
	
	if((NPcutoffProximal==0 && NPcutoffDistal==0) && wrapDistance>0)
	{
		std::cerr << "--wrapDistance requires --NPcutoff!";
		showSwitches();
		return -1;
	}
	
	if(normalCutoff>0 && molIndex.size()==0)
	{
		std::cerr << "--normalCutoff requires --molIndex!";
		showSwitches();
		return -1;
	}
	
	if(molIndex.size()==0 && bendAngleSteps>0)
	{
		std::cerr << "--bendAngleSteps requires --molIndex!";
		showSwitches();
		return -1;
	}
	
	if(molIndex.size()==0 && bendAngleStepSize>0)
	{
		std::cerr << "--bendAngleStepSize requires --molIndex!";
		showSwitches();
		return -1;
	}
	
	if(molIndex.size()==0 && NPmolIndex.size()>0)
	{
		std::cerr << "--NPmolIndex requires --molIndex!";
		showSwitches();
		return -1;
	}
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name.c_str(),std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//size file
	std::string sizeName("size_");
	sizeName+=name;
	sizeName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	
	//shortened names
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	molecule<double, fourVector<int> > *m=System.getMolecule();
	int nTypes=System.readNTypes();
	
	//frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	xyzFormat<double> xyzFile(p,nParticles);
	
	if(frames)
		xyzFile.open(framesName.c_str(), std::ios::in);
	
	int maxChainLength=0;
	//dump molecule information
	for(int mol=0;mol<System.readNMolecules();mol++)
	{
		if(m[mol].readType()==CHAIN)
		{
			//go through each group of chains
			for(int bond=0;bond<m[mol].readNBond();bond++)
			{
				int start=m[mol].getBonds()[bond].s[START];
				int length=m[mol].getBonds()[bond].s[CHAINLENGTH];
				if(length>maxChainLength)
					maxChainLength=length;
				int nChains=m[mol].getBonds()[bond].s[NCHAINS];
				std::cerr << "Chain " << mol << " with " << nChains << " chains of length " << length << ".\n";
				std::vector<int> nType(System.readNTypes(),0);
				for(int i=start;i<start+length*nChains;i++)
					nType[p[i].type]++;
				std::cerr << "Number by type (type|nType): ";
				for(int i=0;i<nType.size()-1;i++)
					std::cerr << i << "|" << nType[i] << ", ";
				std::cerr << nType.size()-1 << "|" << nType.back() << std::endl;
			}
		}
		if(m[mol].readType()==BEAD)
		{
			double radius=m[mol].getConstants()[4];
			int nBeads=m[mol].readNBond();
			std::cerr << "Bead " << mol << " with " << nBeads << " beads of radius " << radius << ".\n";
		}
	}
	
	
	//get CHAIN indices
	std::vector<std::vector<int> > chainIndices;
	for(auto& mol:molIndex)
	{
		//bounds and error checking
		if(mol>=System.readNMolecules() || mol<0)
		{
			std::cerr << "Molecule index --molIndex=" << mol << " is out of bounds [0, " << System.readNMolecules() << ").\n";
			return -1;
		}
		if(m[mol].readType()!=CHAIN)
		{
			std::cerr << "Molecule CHAIN (" << CHAIN << ") type, from --molIndex, is wrong!\n";
			return -1;
		}
		//go through each group of chains
		for(int bond=0;bond<m[mol].readNBond();bond++)
		{
			int start=m[mol].getBonds()[bond].s[START];
			int length=m[mol].getBonds()[bond].s[CHAINLENGTH];
			int nChains=m[mol].getBonds()[bond].s[NCHAINS];
			
			//go through every chain in group
			for(int j=start;j<start+nChains*length;j+=length)
			{
				//find an end that matches our type
				int orientOffset=0;
				for(int k=j;k<j+length;k++)
				{
					if(std::find(type.begin(), type.end(),p[k].type)!=type.end())
						orientOffset=k-j;
				}
				//make sure the end is actually at the end
				if(!(orientOffset==0 || orientOffset==length-1))
				{
					std::cerr << "Warning: Molecule CHAIN (" << CHAIN << ") type, from --molIndex, is not flippable with --type!\n";
					//return -1;
				}
				//now orient and group the indices into a chain
				std::vector<int> chain;
				if(orientOffset==0)
					for(int k=j;k<j+length;k++)
						chain.push_back(k);
				else
					for(int k=j+length-1;k>=j;k--)
						chain.push_back(k);
				
				//add the chain to our data structure
				chainIndices.push_back(chain);
			}
		}
	}
	
	//get BEAD or nanoparticle indices
	std::vector<twoVector<int> > npIndices;
	double npRadius=0;
	//std::vector<double> radii;
	for(auto& NPmol:NPmolIndex)
	{
		//bounds and error checking
		if(NPmol>=System.readNMolecules() || NPmol<0)
		{
			std::cerr << "Molecule index --NPmolIndex=" << NPmol << " is out of bounds [0, " << System.readNMolecules() << ").\n";
			return -1;
		}
		if(m[NPmol].readType()!=BEAD)
		{
			std::cerr << "Molecule BEAD (" << BEAD << ") type, from --NPmolIndex, is wrong!\n";
			return -1;
		}
		
		//go through each BEAD and put it in our data structures
		for(int i=0;i<m[NPmol].readNBond();i++)
		{
			twoVector<int> np;
			np.x=m[NPmol].getBonds()[i].s[0];
			np.y=NPmol;
			if(npIndices.size()!=1)
				npIndices.push_back(np);
			if(npRadius==0)
			{
				npRadius=m[NPmol].getConstants()[4];
			}
			else if(npRadius!=m[NPmol].getConstants()[4])
			{
				std::cerr << "BEAD Molecule " << NPmol << " is not the same radius as " << NPmolIndex[0] << "!\n";
				return -1;
			}
			//radii.push_back(m[NPmol].getConstants()[4]);
		}
	}
	
	threeVector<double> size=System.readSize();
	double time=0;
	
	//we should use size file to determine what the start time is
	if(frames)
	{
		std::cerr << "Using frames!" << std::endl;
		if(sizeFile.is_open())
		{
			double sTime=-1;
			threeVector<double> sCurrent;
			while(sTime<0 && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			size=sCurrent;
			time=sTime;
			std::cerr << "Got size (" << size.x << ", " << size.y << ", " << size.z << ") at time " << time << std::endl;
		}
	}
	
	//current time averaged shell profile
	std::vector<std::vector<double> > proximalShellProfile, distalShellProfile;
	int nProfiles=0;
	double profileTimeSpan=0;
	
	int frameIndex=0;
	double lastSmoothMapOutput=time;
	double lastMapOutput=time;
	
	//start going through frames, or assume last frame is loaded
	while(xyzFile.load())
	{	
		std::cerr << "frame " << frameIndex << "!" << std::endl;
		
		//pMap=position map where = com position or "head" position
		std::vector< position<double> > pMap;//=position
		
		//get positions based on type
		if(type.size()>0)
		{
			//for every "head", index [0] in a chain
			for(const auto& chain:chainIndices)
			{
				position<double> mapPosition;
				mapPosition.x=p[chain[0]].x;
				mapPosition.y=p[chain[0]].y;
				mapPosition.z=p[chain[0]].z;
				if(mapPosition.x==size.x) mapPosition.x-=size.x;
				if(mapPosition.y==size.y) mapPosition.y-=size.y;
				if(mapPosition.z==size.z) mapPosition.z-=size.z;
				pMap.push_back(mapPosition);
			}
		}
		//get positions based on center of mass
		else
		{
			//for every chain
			for(const auto& chain:chainIndices)
			{
				//we are summing up all the positions in each chain index
				position<double> mapPosition;
				mapPosition.x=0;
				mapPosition.y=0;
				mapPosition.z=0;
				for(const auto& index:chain)
				{
					threeVector<double> d, minImg=0;
					//difference between the first particle and current
					d.x=p[index].x-p[chain[0]].x;
					d.y=p[index].y-p[chain[0]].y;
					d.z=p[index].z-p[chain[0]].z;
					//minimum image
					if(d.x>size.x/2.0) minImg.x=size.x;
					if(d.y>size.y/2.0) minImg.y=size.y;
					if(d.z>size.z/2.0) minImg.z=size.z;
					if(d.x<-size.x/2.0) minImg.x=-size.x;
					if(d.y<-size.y/2.0) minImg.y=-size.y;
					if(d.z<-size.z/2.0) minImg.z=-size.z;
					//sum them
					mapPosition.x+=(p[index].x-minImg.x)/static_cast<double>(chain.size());
					mapPosition.y+=(p[index].y-minImg.y)/static_cast<double>(chain.size());
					mapPosition.z+=(p[index].z-minImg.z)/static_cast<double>(chain.size());
				}
				//now make sure it is in our box
				if(mapPosition.x>=size.x) mapPosition.x-=size.x;
				if(mapPosition.y>=size.y) mapPosition.y-=size.y;
				if(mapPosition.z>=size.z) mapPosition.z-=size.z;
				if(mapPosition.x<0) mapPosition.x+=size.x;
				if(mapPosition.y<0) mapPosition.y+=size.y;
				if(mapPosition.z<0) mapPosition.z+=size.z;
				//put it in our data structure
				pMap.push_back(mapPosition);
			}
		}
		
		//oMap=bend order parameter map where 
		// {x,y,z,t}={endToEnd-com={x,y,z},ord. para.={t}}
		std::vector< threeVector<double> > oMap;//={endToEnd-com={x,y,z},ord. para.={t}}
		
		//find orientation of every chain
		for(const auto& chain:chainIndices)
		{
			//our orientation
			threeVector<double> mapOrient=0;
			//alternatively, we could find orientation along the chain rather than end-to-end
			//for(const auto& index:chain)
			{
				threeVector<double> d, minImg=0;
				//difference between end-to-end vector
				mapOrient.x=p[chain.front()].x-p[chain.back()].x;
				mapOrient.y=p[chain.front()].y-p[chain.back()].y;
				mapOrient.z=p[chain.front()].z-p[chain.back()].z;
				//make sure they are next to one another in the box
				if(mapOrient.x>size.x/2.0) mapOrient.x-=size.x;
				if(mapOrient.y>size.y/2.0) mapOrient.y-=size.y;
				if(mapOrient.z>size.z/2.0) mapOrient.z-=size.z;
				if(mapOrient.x<-size.x/2.0) mapOrient.x+=size.x;
				if(mapOrient.y<-size.y/2.0) mapOrient.y+=size.y;
				if(mapOrient.z<-size.z/2.0) mapOrient.z+=size.z;
			}
			//make it a unit vector so we don't have to compute this each time
			mapOrient=unitVector(mapOrient);
			//put it in our data structure
			oMap.push_back(mapOrient);
		}
		
		//now we need the order parameter
		//bendMap=bend order parameter where =ord. para.=cos(theta)
		std::vector<double> bendMap;//=costheta
		
		//find order parameter of every chain
		for(const auto& chain:chainIndices)
		{
			//we are averaging it over the chain
			double bendValue=0;
			for(auto i=chain.begin();i<chain.end()-2;i++)
			{
				threeVector<double> da,db;
				//vector A->B=DA
				da.x=p[*(i)].x-p[*(i+1)].x;
				da.y=p[*(i)].y-p[*(i+1)].y;
				da.z=p[*(i)].z-p[*(i+1)].z;
				if(da.x>size.x/2.0) da.x-=size.x;
				if(da.x<-size.x/2.0) da.x+=size.x;
				if(da.y>size.y/2.0) da.y-=size.y;
				if(da.y<-size.y/2.0) da.y+=size.y;
				if(da.z>size.z/2.0) da.z-=size.z;
				if(da.z<-size.z/2.0) da.z+=size.z;
				
				//vector B->C=DB
				db.x=p[*(i+1)].x-p[*(i+2)].x;
				db.y=p[*(i+1)].y-p[*(i+2)].y;
				db.z=p[*(i+1)].z-p[*(i+2)].z;
				if(db.x>size.x/2.0) db.x-=size.x;
				if(db.x<-size.x/2.0) db.x+=size.x;
				if(db.y>size.y/2.0) db.y-=size.y;
				if(db.y<-size.y/2.0) db.y+=size.y;
				if(db.z>size.z/2.0) db.z-=size.z;
				if(db.z<-size.z/2.0) db.z+=size.z;
				
				//make them unit vectors
				da=unitVector(da);
				db=unitVector(db);
				
				//should be a value from -1 to 1, from <DA|DB>=|DA||DB|cos(theta)
				double costheta=(da.x*db.x)+(da.y*db.y)+(da.z*db.z);
				bendValue+=costheta;
			}
			//average the order parameter
			bendValue/=static_cast<double>(chain.size()-2);
			//put it in our data structure
			bendMap.push_back(bendValue);
		}
		
		//the mean (average) value of the order parameter
		double mean=0;
		for(const auto& bendValue:bendMap)
			mean+=bendValue;
		if(bendMap.size()>0)
			mean/=static_cast<double>(bendMap.size());
		
		//the standard deviation of the order parameter
		double stdDev=0;
		for(const auto& bendValue:bendMap)
			stdDev+=(bendValue-mean)*(bendValue-mean);
		if(bendMap.size()-1>0)
			stdDev/=static_cast<double>(bendMap.size()-1);
		stdDev=sqrt(stdDev);
		
		//output the values to a file
		if(molIndex.size()>0)
		{
			std::string chainOrderName("chainOrder_");
			chainOrderName+=name;
			chainOrderName+=".dat";
			std::cerr << "Writing " << chainOrderName << " at time " << time << "..." << std::endl;
			std::fstream chainOrder(chainOrderName.c_str(), std::ios::out | std::ios::app);
			if(chainOrder.is_open())
				chainOrder << time << '\t' << mean << '\t' << stdDev << std::endl;
			//std::cout << "Mean cos(theta)=" << mean << " and standard deviation=" << stdDev << std::endl;
		}
		
		//find the minima and maxima values of order parameter for variable step lengths
		double min=1.0, max=-1.0, stepSize=0.0;
		int steps=0;
		if(bendAngleSteps>0)
		{
			for(const auto& bendValue:bendMap)
			{
				min=(bendValue<min)?bendValue:min;
				max=(bendValue>max)?bendValue:max;
			}
			stepSize=(max-min)/static_cast<double>(bendAngleSteps);
			steps=bendAngleSteps;
		}
		
		//used fixed steps and entire range [-1,1]
		if(bendAngleStepSize>0)
		{
			stepSize=bendAngleStepSize;
			min=-1.0;
			max=1.0;
			steps=floor((max-min)/stepSize);
		}
		
		//output a histogram of order parameter
		if(bendAngleSteps>0 || bendAngleStepSize>0)
		{	
			std::string chainOrderHistName("chainOrderHist_");
			chainOrderHistName+=name;
			chainOrderHistName+=".dat";
			
			std::cerr << "Writing " << chainOrderHistName << " at time " << time << "..." << std::endl;
			std::fstream chainOrderHistFile;
			chainOrderHistFile.open(chainOrderHistName.c_str(), std::ios::out | std::ios::app);
			if(chainOrderHistFile.is_open())
			{
				std::vector<double> chainOrderHist(steps,0);
				for(const auto& bendValue:bendMap)
					chainOrderHist[floor((bendValue-min)/stepSize)]++;
				for(int i=0;i<chainOrderHist.size();i++)
				{
					chainOrderHist[i]/=static_cast<double>(bendMap.size())*stepSize;
					chainOrderHistFile << ((static_cast<double>(i)+1)*stepSize)+min << '\t' << chainOrderHist[i] << '\n';
				}
				chainOrderHistFile << std::endl;
			}
		}
		
		threeVector<double> NPcom=0;
		for(auto &NPindex:npIndices)
		{
			NPcom.x+=p[NPindex.x].x;
			NPcom.y+=p[NPindex.x].y;
			NPcom.z+=p[NPindex.x].z;
		}
		if(npIndices.size()>0)
		{
			NPcom/=static_cast<double>(npIndices.size());
			//shift the coordinates to the center of the simulation box
			NPcom.x+=size.x/2.0;
			NPcom.y+=size.y/2.0;
			NPcom.z+=size.z/2.0;
		}
		
		
			
		/*Top layer of proximal leaflet, wrapDistance is right bracket ']':
		       
		       ====   ] <---------This part, in contact with NP
		      //  \\                              -
		     || NP ||
		      \\__//
		       ----
		  ____      _____Membrane
		      \    /
		       |__|
		      //  \\  ] <---------This part, in contact with NP
		     || NP ||
		      \\__//
		       ----
		
		  ___   __   _____Membrane
		     \ /  \ /  ] <---------This part, in contact with NP
		     || NP ||
		      \\__//
		       ----
		
			__		
		       /  \   Membrane
		  ___ | NP |_____  ] <---This part, in contact with NP
		      \\__//
		       ----                            z
		                                       ^
		                                       |
		       /  \                            *--> x or y
		      | NP |   Membrane
		       \__/ 
		 ---------------- ] <---This part, in contact with NP
		*/
		//indices for the proximal leaflet
		std::vector<int> proximal,distal;
		std::vector<int> proximal90,distal90;
		double wrapCostheta=0;
			
		//unit vector along z direction
		threeVector<double> zHat=0;
			
		//wrapping cos(angle) with cutoff
		if(NPcutoffProximal>0 && NPcutoffDistal>0 && wrapDistance>0)
		{
			//cutoffSqr=cutoff*cutoff
			double NPcutoffSqr=npRadius+NPcutoffProximal;
			NPcutoffSqr*=NPcutoffSqr;
			
			//go through every particle in the map to get the maximum and minimum phi component of
			// the proximal leaflet
			//i is a consistent index between pMap, bendMap, chainIndices, and oMap
			for(int i=0;i<pMap.size() && pMap.size()==oMap.size() && pMap.size()==bendMap.size();i++)
			{
				//go through each nanoparticle index
				for(auto& npIndex:npIndices)
				{
					//get the distance between nanoparticle and each chain center of mass
					threeVector<double> d;
					d.x=p[npIndex.x].x-pMap[i].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[npIndex.x].y-pMap[i].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[npIndex.x].z-pMap[i].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					double r=d.x*d.x+d.y*d.y+d.z*d.z;
					
					//find the distance from  (P=proximal and D=distal leaflets) {
					/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
					         z______PDDDDDDDDDDDDDDDDDDDD  /
					npCutoff ^     P\D  } wrapDistance
					         |       P\D
					         |         P\D
					npRadius |__        P|D
					         |  \       P|D
					         |   |       P\D
					         |   |       P|D
					       NP*-------> x or y
					
					*/
					//if it is nearby
					if(r<NPcutoffSqr)
					{
						threeVector<double> uD=unitVector(d);
						double costheta=dotProduct(oMap[i],uD);
						d=unitVector(d);
						//is the lipid facing the NP? i.e. proximal 
						//(HP and TP are head and tail proximal beads)
						//(HD and TD are head and tail distal beads)
						//vectors: T*-->H = oMap
						//vectors: NP*-->H = d
						/*{{                  
						                    ^
						                 HD/ 
						                TD/ } <-- costheta= -1
						         z     TD/ 
						npCutoff ^      *  oMap 
						         |   TP/ 
						         |  TP/   } <-- costheta= 1
						npRadius |_HP/      
						         |  u\                 -
						         |    \                -
						         | ^   \               -
						         |/ <-d |
						       NP*-------> x or y
						*/
						if(costheta>wrapOrientCosTheta)
							proximal.push_back(i);
					}
				}
			}
			
			
			if(proximal.size()==0)
			{
				std::cerr << "No proximal lipids Found! No files written!" << std::endl;
			}
			
			//go through all chains in proximal leaflet to find orientation
			if(proximal.size()>0)
			for(auto &i:proximal)
			{
				//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets) {
				/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
				         z      PDDDDDDDDDDDDDDDDDDDDD /
				npCutoff ^     PD  } wrapDistance
				         |    / / PD
				         |   y y<--PD
				npRadius |_      <--PD
				         | \      <--PD
				         |  |     <--PD
				         |   |     <--PD                         ^
				       NP*-------> x or y  ----> overall zHat {  |
				             |     <--PD                         *
				            |     <--PD
				           /      <--PD
				         -       <--PD
				         ^^^^^^^^<-PD
				         |||||||\\PD
				         PPPPPPPPPD
				         DDDDDDDDD
				*/
				zHat.x+=oMap[i].x;
				zHat.y+=oMap[i].y;
				zHat.z+=oMap[i].z;
			}
			zHat=unitVector(zHat);
			std::cerr << "zHat: (" << zHat.x << ',' << zHat.y << ',' << zHat.z << ')' << std::endl;

			//sort proximal leaflet along z direction
			threeVector<double> center;
			center.x=p[npIndices[0].x].x;
			center.y=p[npIndices[0].x].y;
			center.z=p[npIndices[0].x].z;
			std::sort(proximal.begin(), proximal.end(), sortByZ< std::vector< position<double> >, threeVector<double> >(pMap,zHat,center,size));
			
			//for average cos(angle)
			int nWrapCostheta=0;
			int nAnulus=0;
			threeVector<double> anulusCom=0;
			double anulusDiameter=0;
			//std::fstream anulusDump("anulusDump.xyz",std::ios::app | std::ios::out);
			//anulusDump << 2000 << "\nasdf\n";
			//go through all chains in proximal leaflet
			if(proximal.size()>0)
			for(auto &i:proximal)
			{
				//get distance to top particle in proximal leaflet
				double dz;
				dz=pMap[proximal.front()].z-pMap[i].z;
				if(dz>size.z/2.0) dz-=size.z;
				if(dz<-size.z/2.0) dz+=size.z;
				//is it within the layer thickness?
				if(abs(dz)<wrapDistance)
				{
					threeVector<double> d;
					//get the angle it forms with the nanoparticles
					//go through each nanoparticle index
					for(auto& npIndex:npIndices)
					{
						//get the distance between nanoparticle and each chain center of mass
						//threeVector<double> d;
						d.x=p[npIndex.x].x-pMap[i].x;
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						d.y=p[npIndex.x].y-pMap[i].y;
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						d.z=p[npIndex.x].z-pMap[i].z;
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						
						//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets)
						/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
						         z      PDDDDDDDDDDDDDDDDDDDDD /
						npCutoff ^     PD  } wrapDistance
						         |     ^ PD
						         |    /    PD
						npRadius |__ /      PD
						         |  \       PD
						         | / |       PD
						         |/  |<- costheta
						       NP*-------> x or y
						
						*/
						d=unitVector(d);
						double costheta=dotProduct(d,zHat);
						wrapCostheta+=costheta;
						nWrapCostheta++;
					}
					//anulusDump << 9 << '\t' << pMap[i].x << '\t' << pMap[i].y << '\t' << pMap[i].z << std::endl;
					threeVector<double> recenter=0;
					d.x=pMap[i].x-p[npIndices[0].x].x;
					if(d.x<-size.x/2.0) recenter.x=size.x;
					if(d.x>size.x/2.0) recenter.x=-size.x;
					d.y=pMap[i].y-p[npIndices[0].x].y;
					if(d.y<-size.y/2.0) recenter.y=size.y;
					if(d.y>size.y/2.0) recenter.y=-size.y;
					d.z=pMap[i].z-p[npIndices[0].x].z;
					if(d.z<-size.z/2.0) recenter.z=size.z;
					if(d.z>size.z/2.0) recenter.z=-size.z;
					anulusCom.x+=pMap[i].x+recenter.x;
					anulusCom.y+=pMap[i].y+recenter.y;
					anulusCom.z+=pMap[i].z+recenter.z;
					nAnulus+=1;
				}
			}
			if(nAnulus>0)
			{
				anulusCom.x/=static_cast<double>(nAnulus);
				anulusCom.y/=static_cast<double>(nAnulus);
				anulusCom.z/=static_cast<double>(nAnulus);
			}
			anulusCom.x=p[npIndices[0].x].x;
			anulusCom.y=p[npIndices[0].x].y;
			//for(int asdf=nAnulus;asdf<2000;asdf++)
			//	anulusDump << 9 << '\t' << p[npIndices[0].x].x << '\t' << p[npIndices[0].x].y << '\t' << p[npIndices[0].x].z << std::endl;
			//anulusCom.x=(npRadius+4.0)*zHat.x+p[npIndices[0].x].x;
			//anulusCom.y=(npRadius+4.0)*zHat.y+p[npIndices[0].x].y;
			//anulusCom.z=(npRadius+4.0)*zHat.z+p[npIndices[0].x].z;
			std::cerr << "nWrapCostheta: " << nWrapCostheta << std::endl;
			std::cerr << "nAnulus: " << nAnulus << std::endl;
			if(nWrapCostheta>0)
				wrapCostheta/=static_cast<double>(nWrapCostheta);
				
			double potentialNear=0, potentialNeck=0;
			//go through each nanoparticle index for potential
			if(potentialWrap)
			for(auto& npIndex:npIndices)
			{
				double *C=m[npIndex.y].getConstants();
				double cutoffRadSqr=C[0];//radius+rc
				cutoffRadSqr*=cutoffRadSqr;//(radius+rc)^2
				
				//go through every particle in the map to get the maximum and minimum phi component of
				// the proximal leaflet
				#pragma omp parallel for reduction(+:potentialNear) reduction(+:potentialNeck)
				for(int i=0;i<chainIndices.size();i++)
				{
					//get the distance between NP and chain head
					threeVector<double> d;
					d.x=p[npIndex.x].x-pMap[i].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[npIndex.x].y-pMap[i].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[npIndex.x].z-pMap[i].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					double r=d.x*d.x+d.y*d.y+d.z*d.z;
					
					d=unitVector(d);
					double costheta=dotProduct(d,zHat);
					//wrapCostheta+=costheta;
					
					//for outputting potential at NP surface boundary and "neck" anulus
					//---------------------------------------------------------------
					// r and d are modified here!
					//---------------------------------------------------------------
					if(r<cutoffRadSqr)
					{
						for(auto &ci:chainIndices[i])
						{
							d.x=p[npIndex.x].x-p[ci].x;
							if(d.x>size.x/2.0) d.x-=size.x;
							if(d.x<-size.x/2.0) d.x+=size.x;
							d.y=p[npIndex.x].y-p[ci].y;
							if(d.y>size.y/2.0) d.y-=size.y;
							if(d.y<-size.y/2.0) d.y+=size.y;
							d.z=p[npIndex.x].z-p[ci].z;
							if(d.z>size.z/2.0) d.z-=size.z;
							if(d.z<-size.z/2.0) d.z+=size.z;
							double r=d.x*d.x+d.y*d.y+d.z*d.z;
							
							int cindex=nBEADCONST*((p[ci].type*nTypes)+p[npIndex.x].type);
							if(costheta>wrapCostheta)
								potentialNear+=beadPotential(d,&(C[cindex]), cutoffRadSqr);
							else
								potentialNeck+=beadPotential(d,&(C[cindex]), cutoffRadSqr);
						}
					} 
				}
			}
			if(potentialWrap)
			{
				std::fstream(std::string("potentialWrap_")+name+".dat", std::ios::app | std::ios::out) 
					<< time << '\t' << potentialNear << '\t' << potentialNeck << std::endl;
			}
			//for standard deviation of cos(angle)
			double wrapCosthetaStd=0;
			nWrapCostheta=0;
			
			//go through all chains in proximal leaflet
			if(proximal.size()>0)
			for(auto &i:proximal)
			{
				//get distance to top particle in proximal leaflet
				double dz;
				dz=pMap[proximal.front()].z-pMap[i].z;
				if(dz>size.z/2.0) dz-=size.z;
				if(dz<-size.z/2.0) dz+=size.z;
				//is it within the layer thickness?
				if(abs(dz)<wrapDistance)
				{
					threeVector<double> d;
					//get the angle it forms with the nanoparticles
					//go through each nanoparticle index
					for(auto& npIndex:npIndices)
					{
						//get the distance between nanoparticle and each chain center of mass
						//threeVector<double> d;
						d.x=p[npIndex.x].x-pMap[i].x;
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						d.y=p[npIndex.x].y-pMap[i].y;
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						d.z=p[npIndex.x].z-pMap[i].z;
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						
						//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets)
						/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
						         z      PDDDDDDDDDDDDDDDDDDDDD /
						npCutoff ^     PD  } wrapDistance
						         |     ^ PD
						         |    /    PD
						npRadius |__ /      PD
						         |  \       PD
						         | / |       PD
						         |/  |<- costheta
						       NP*-------> x or y
						
						*/
						d=unitVector(d);
						double costheta=dotProduct(d,zHat);
						wrapCosthetaStd+=pow(wrapCostheta-costheta,2);
						nWrapCostheta++;
					}
					d.x=anulusCom.x-pMap[i].x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					if(d.x>size.x/2.0) d.x-=size.x;
					d.y=anulusCom.y-pMap[i].y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					if(d.y>size.y/2.0) d.y-=size.y;
					d.z=anulusCom.z-pMap[i].z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					if(d.z>size.z/2.0) d.z-=size.z;
					anulusDiameter+=sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
				}
			}
			
			if(nAnulus>0)
			{
				anulusDiameter/=static_cast<double>(nAnulus);
				std::string anulusName("anulusDiameter_");
				anulusName+=name;
				anulusName+=".dat";
				std::cerr << "Writing " << anulusName << " at time " << time << "..." << std::endl;
				std::fstream anulusFile(anulusName.c_str(), std::ios::out | std::ios::app);
				if(anulusFile.is_open())
					anulusFile << time << '\t' << anulusDiameter << std::endl;
			}
			std::string anulusComName("anulusCOM_");
			anulusComName+=name;
			anulusComName+=".xyz";
			std::cerr << "Writing " << anulusComName << " at time " << time << "..." << std::endl;
			std::fstream anulusComFile(anulusComName.c_str(), std::ios::out | std::ios::app);
			if(anulusComFile.is_open())
			{
				anulusComFile << 1 << "\nneckAnulusCom\n";
				anulusComFile << 10 << '\t' << anulusCom.x << '\t' << anulusCom.y << '\t' << anulusCom.z << std::endl;
			}
			
			if(nWrapCostheta-1>0)
				wrapCosthetaStd/=static_cast<double>(nWrapCostheta-1);
			wrapCosthetaStd=sqrt(wrapCosthetaStd);
			
			//output wrapping angle to file
			if(proximal.size()>0)
			{
				std::string wrapName("wrapCosAngle_");
				wrapName+=name;
				wrapName+=".dat";
				std::cerr << "Writing " << wrapName << " at time " << time << "..." << std::endl;
				std::fstream wrapCosAngle(wrapName.c_str(), std::ios::out | std::ios::app);
				if(wrapCosAngle.is_open())
					wrapCosAngle << time << '\t' << wrapCostheta << '\t' << wrapCosthetaStd << std::endl;
			}
			
			//cutoffSqr=cutoff*cutoff
			double NPcutoffSqr2=npRadius+NPcutoffDistal;
			NPcutoffSqr2*=NPcutoffSqr2;
			
			//go through every particle in the map to get the distal leaflet
			//i is a consistent index between pMap, bendMap, chainIndices, and oMap
			if(proximal.size()>0)
			for(int i=0;i<pMap.size() && pMap.size()==oMap.size() && pMap.size()==bendMap.size();i++)
			{
				//go through each nanoparticle index
				for(auto& npIndex:npIndices)
				{
					//get the distance between nanoparticle and each chain center of mass
					threeVector<double> d;
					d.x=p[npIndex.x].x-pMap[i].x;
					if(d.x>size.x/2.0) d.x-=size.x;
					if(d.x<-size.x/2.0) d.x+=size.x;
					d.y=p[npIndex.x].y-pMap[i].y;
					if(d.y>size.y/2.0) d.y-=size.y;
					if(d.y<-size.y/2.0) d.y+=size.y;
					d.z=p[npIndex.x].z-pMap[i].z;
					if(d.z>size.z/2.0) d.z-=size.z;
					if(d.z<-size.z/2.0) d.z+=size.z;
					double r=d.x*d.x+d.y*d.y+d.z*d.z;
					
					//find the cos(angle) it makes with the z component
					d=unitVector(d);
					double costheta=dotProduct(d,zHat);
					
					//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets)
					/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
					         z      PDDDDDDDDDDDDDDDDDDDDD /
					npCutoff ^     PD  } wrapDistance
					         |     ^ PD
					         |    /    PD
					npRadius |__ /      PD
					         |  \       PD
					         | / |       PD
					         |/  |<- costheta
					       NP*-------> x or y
					
					*/
					//if it is nearby and intersects the proximal wrapping angle
					if(r<NPcutoffSqr2 && costheta>wrapCostheta)
					{
						threeVector<double> uD=unitVector(d);
						costheta=dotProduct(oMap[i],uD);
						//is the lipid facing away from the NP? i.e. distal 
						//(HP and TP are head and tail proximal beads)
						//(HD and TD are head and tail distal beads)
						//vectors: T*-->H = oMap
						//vectors: NP*-->H = d
						/*{{                  
						npCutoff2           ^
						                 HD/ 
						                TD/ } <-- costheta= -1
						         z     TD/ 
						npCutoff ^      *  oMap 
						         |   TP/ 
						         |  TP/   } <-- costheta= 1
						npRadius |_HP/      
						         |  u\                 -
						         |    \                -
						         | ^   \               -
						         |/ <-d |
						       NP*-------> x or y
						*/
						if(costheta<=-wrapOrientCosTheta)
						{
							distal.push_back(i);
						}
					}
				}
			}
			
			double proximalBend=0;
			if(proximal.size()>0)
			for(auto &i:proximal)
				proximalBend+=bendMap[i];
			if(proximal.size()>0)
				proximalBend/=static_cast<double>(proximal.size());
			double proximalBendStd=0;
			if(proximal.size()>0)
			for(auto &i:proximal)
				proximalBendStd+=pow(proximalBend-bendMap[i],2);
			if(proximal.size()>0)
				proximalBendStd/=static_cast<double>(proximal.size()-1);
			proximalBendStd=sqrt(proximalBendStd);
			
			//output proximal order parameter to file
			if(proximal.size()>0)
			{
				std::string proximalName("proximalOrder_");
				proximalName+=name;
				proximalName+=".dat";
				std::cerr << "Writing " << proximalName << " at time " << time << "..." << std::endl;
				std::fstream proximalOrder(proximalName.c_str(), std::ios::out | std::ios::app);
				if(proximalOrder.is_open())
					proximalOrder << time << '\t' << proximalBend << '\t' << proximalBendStd << std::endl;
			}
			
			double distalBend=0;
			for(auto &i:distal)
				distalBend+=bendMap[i];
			if(distal.size()>0)
				distalBend/=static_cast<double>(distal.size());
			double distalBendStd=0;
			if(distal.size()>0)
			for(auto &i:distal)
				distalBendStd+=pow(distalBend-bendMap[i],2);
			if(distal.size()>0)
				distalBendStd/=static_cast<double>(distal.size()-1);
			distalBendStd=sqrt(distalBendStd);
			
			//output distal order parameter to file
			if(distal.size()>0)
			{
				std::string distalName("distalOrder_");
				distalName+=name;
				distalName+=".dat";
				std::cerr << "Writing " << distalName << " at time " << time << "..." << std::endl;
				std::fstream distalOrder(distalName.c_str(), std::ios::out | std::ios::app);
				if(distalOrder.is_open())
					distalOrder << time << '\t' << distalBend << '\t' << distalBendStd << std::endl;
			}
			
			if(proximal.size()>0 && distal.size()>0 && map)
			{
				std::string proximalXyzName("proximal_");
				proximalXyzName+=name;
				proximalXyzName+=".xyz";
				std::cerr << "Writing " << proximalXyzName << " at time " << time << "..." << std::endl;
				std::fstream proximalXyz(proximalXyzName.c_str(), std::ios::out | std::ios::app);
				if(proximalXyz.is_open())
				{
					proximalXyz << proximal.size() << "\ntest\n";
					for(auto &i:proximal)
						proximalXyz << 1 << '\t' << pMap[i].x << '\t' << pMap[i].y << '\t' << pMap[i].z << '\n';
				}
				std::string distalXyzName("distal_");
				distalXyzName+=name;
				distalXyzName+=".xyz";
				std::cerr << "Writing " << distalXyzName << " at time " << time << "..." << std::endl;
				std::fstream distalXyz(distalXyzName.c_str(), std::ios::out | std::ios::app);
				if(distalXyz.is_open())
				{
					distalXyz << distal.size() << "\ntest\n";
					for(auto &i:distal)
						distalXyz << 1 << '\t' << pMap[i].x << '\t' << pMap[i].y << '\t' << pMap[i].z << '\n';
				}
			}
		}
		
		/*smoothed map in a vmd "xyz" file {
			z
			^ 
			|
			*--> x or y
		     Normal VMD:                Colors as "type":
		  ____      ____Membrane         0001      1000
		      \    /                         5    5
		       |  |                           5  5
		      /    \              ---->      5    5
		     |  NP  |                       5  NP  5
		      \    /                         5    5
		       ----                           5555
		
		Smothing (* is for center lipid, + is for nearby lipid values):
			y
			^ 
			|
			*--> x
		    ___    
		   / + \   } smoothCutoff length
		  |+ * +|
		   \___/
		
		*/
		if(smoothCutoff>0 && stepSize>0 && time-lastSmoothMapOutput>frequency)
		{
			bendRange bR(pMap,bendMap,smoothCutoff);
			Cell<double> pairInteractions(&pMap[0], pMap.size(), smoothCutoff, size);
			pairInteractions.build();
			pairInteractions.twoWayComparePeriodicIndex(bR);
			bR.finalize();
			
			//double sMin=1.0,sMax=-1.0;
			//for(int i=0;i<bendAvg.size();i++)
			//{
			//	sMin=(sMin>bendAvg[i])?bendAvg[i]:sMin;
			//	sMax=(sMax<bendAvg[i])?bendAvg[i]:sMax;
			//}
			std::stringstream chainOrderMapName;
			chainOrderMapName << "chainOrderSmoothMap_" << name << "_" << time << ".csv";
			std::string cOMN;
			chainOrderMapName >> cOMN;
			std::cerr << "Writing " << cOMN << " at time " << time << "..." << std::endl;
			
			std::fstream chainOrderMapFile;
			chainOrderMapFile.open(cOMN.c_str(), std::ios::out);
			if(chainOrderMapFile.is_open())
			{
				chainOrderMapFile << "cosThetaAvg,cosTheta,x,y,z,type2,nX,nY,nZ,type,proximal\n";
				std::vector<double> bendAvg=bR.getBendAverages();
				for(auto& i:proximal)
				{
					chainOrderMapFile << bendAvg[i] << ',';
					chainOrderMapFile << bendMap[i] << ',';
					position<double> pOut=pMap[i];
					pOut.x-=NPcom.x;
					pOut.y-=NPcom.y;
					pOut.z-=NPcom.z;
					while(pOut.x<0)pOut.x+=size.x;
					while(pOut.y<0)pOut.y+=size.y;
					while(pOut.z<0)pOut.z+=size.z;
					chainOrderMapFile << pOut.x << ',' << pOut.y << ',' << pOut.z << ',' << pOut.type << ',';
					chainOrderMapFile << oMap[i].x << ',' << oMap[i].y << ',' << oMap[i].z << ',' << pOut.type << ",1\n";
				}
				for(auto& i:distal)
				{
					chainOrderMapFile << bendAvg[i] << ',';
					chainOrderMapFile << bendMap[i] << ',';
					position<double> pOut=pMap[i];
					pOut.x-=NPcom.x;
					pOut.y-=NPcom.y;
					pOut.z-=NPcom.z;
					while(pOut.x<0)pOut.x+=size.x;
					while(pOut.y<0)pOut.y+=size.y;
					while(pOut.z<0)pOut.z+=size.z;
					chainOrderMapFile << pOut.x << ',' << pOut.y << ',' << pOut.z << ',' << pOut.type << ',';
					chainOrderMapFile << oMap[i].x << ',' << oMap[i].y << ',' << oMap[i].z << ',' << pOut.type << ",0\n";
				}
			}
			lastSmoothMapOutput=time;
		}
		
		/*map in a vmd "xyz" file
			z
			^ 
			|
			*--> x or y
		     Normal VMD:                Colors as "type":
		  ____      ____Membrane         0000      0000
		      \    /                         5    5
		       |  |                           7  7
		      /    \              ---->      5    5
		     |  NP  |                       3  NP  3
		      \    /                         5    3
		       ----                           4355
		*/
		if(map && stepSize>0 && time-lastMapOutput>frequency && smoothCutoff==0)
		{
			std::stringstream chainOrderMapName;
			chainOrderMapName << "chainOrderMap_" << name << "_" << time << ".csv";
			std::string cOMN;
			chainOrderMapName >> cOMN;
			std::cerr << "Writing " << cOMN << " at time " << time << "..." << std::endl;
			
			std::fstream chainOrderMapFile;
			chainOrderMapFile.open(cOMN.c_str(), std::ios::out);
			if(chainOrderMapFile.is_open())
			{
				chainOrderMapFile << "cosTheta,x,y,z,nX,nY,nZ,type,proximal\n";
				for(auto& i:proximal)
				{
					chainOrderMapFile << bendMap[i] << ',';
					position<double> pOut=pMap[i];
					pOut.x-=NPcom.x;
					pOut.y-=NPcom.y;
					pOut.z-=NPcom.z;
					while(pOut.x<0)pOut.x+=size.x;
					while(pOut.y<0)pOut.y+=size.y;
					while(pOut.z<0)pOut.z+=size.z;
					chainOrderMapFile << pOut.x << ',' << pOut.y << ',' << pOut.z << ',';
					chainOrderMapFile << oMap[i].x << ',' << oMap[i].y << ',' << oMap[i].z << ',' << pOut.type << ",1\n";
				}
				for(auto& i:distal)
				{
					chainOrderMapFile << bendMap[i] << ',';
					position<double> pOut=pMap[i];
					pOut.x-=NPcom.x;
					pOut.y-=NPcom.y;
					pOut.z-=NPcom.z;
					while(pOut.x<0)pOut.x+=size.x;
					while(pOut.y<0)pOut.y+=size.y;
					while(pOut.z<0)pOut.z+=size.z;
					chainOrderMapFile << pOut.x << ',' << pOut.y << ',' << pOut.z << ',';
					chainOrderMapFile << oMap[i].x << ',' << oMap[i].y << ',' << oMap[i].z << ',' << pOut.type << ",1\n";
				}
			}
			lastMapOutput=time;
		}
		//get the proximal and distal leaflets using center of mass assuming a closed vesicle
		if(vesicle)
		{
			//this needs to be outside to use in place of p[npIndex]
			threeVector<double> com=0;
			//This is for the 0.9 costheta search restriction
			zHat.x=0;
			zHat.y=0;
			zHat.z=1.0;
			//vesicle center of mass
			for(auto& pM:pMap)
			{
				//get the distance between first particle and each chain center of mass
				threeVector<double> d;
				d.x=pM.x-pMap[0].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=pM.y-pMap[0].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=pM.z-pMap[0].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				com.x+=pM.x;
				com.y+=pM.y;
				com.z+=pM.z;
			}
			com.x/=static_cast<double>(pMap.size());
			com.y/=static_cast<double>(pMap.size());
			com.z/=static_cast<double>(pMap.size());
			while(com.x<0) com.x+=size.x;
			while(com.y<0) com.y+=size.y;
			while(com.z<0) com.z+=size.z;
			while(com.x>=size.x) com.x-=size.x;
			while(com.y>=size.y) com.y-=size.y;
			while(com.z>=size.z) com.z-=size.z;
			
			//go through every particle in the map to get the maximum and minimum phi component of
			// the proximal leaflet
			//i is a consistent index between pMap, bendMap, chainIndices, and oMap
			for(int i=0;i<pMap.size() && pMap.size()==oMap.size() && pMap.size()==bendMap.size();i++)
			{
				//get the distance between nanoparticle and each chain center of mass
				threeVector<double> d;
				d.x=com.x-pMap[i].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=com.y-pMap[i].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=com.z-pMap[i].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				
				//2 conditions, FIRST:
				//find the distance from  (P=proximal and D=distal leaflets) 
				/*        DDDDDDD    <--- Bilayer
				         z-------D  
					 ^PPPPPPP\D  
				         |       P\D
				         |         P\D
				         |          P|D
				         |          P|D
				         |           P\D
				         |           P|D
				      COM*-------> x or y
				
				*/
				double r=d.x*d.x+d.y*d.y+d.z*d.z;
				
				//d=unitVector(d);
				d/=sqrt(r);
				//if(r<lipoCutoffSqr)
				{
					double costheta=dotProduct(oMap[i],d);
					//SECOND:
					//is the lipid facing the COM? i.e. proximal 
					//(HP and TP are head and tail proximal beads)
					//(HD and TD are head and tail distal beads)
					//vectors: T*-->H = oMap
					//vectors: COM*-->H = d
					/*{{                  
					                    ^
					                 HD/ 
					                TD/ } <-- costheta= -1
					         z     TD/ 
					         ^      *  oMap 
					         |   TP/ 
					         |  TP/   } <-- costheta= 1
					vRadius  |_HP/      
					         |  u\                 -
					         |    \                -
					         | ^   \               -
					         |/ <-d |
					      COM*-------> x or y
					*/
					if(costheta>wrapOrientCosTheta)
						proximal.push_back(i);
					else
						distal.push_back(i);
					//ALTERNATIVE:
					//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets)
					/*        DDDDDD   <--- Bilayer
					         zPPPPPPD
					         ^      PD  
					         |     ^  PD
					         |    /    PD
					         |   /      PD
					         |  /       PD
					         | /\        PD
					         |/  \<- costheta
					      COM*---|---> x or y
					
					*/
					double costheta2=dotProduct(d,zHat);
					//if it is nearby and within the lower costheta range
					if(costheta2>0.9)
					{
						if(costheta>wrapOrientCosTheta)
							proximal90.push_back(i);
						else
							distal90.push_back(i);
					}
				}
			}
			std::cerr << proximal.size() << "\t" << proximal90.size() << std::endl;
			std::cerr << distal.size() << "\t" << distal90.size() << std::endl;
			double proximalBend=0;
			if(proximal.size()>0)
			for(auto &i:proximal)
				proximalBend+=bendMap[i];
			if(proximal.size()>0)
				proximalBend/=static_cast<double>(proximal.size());
			double proximalBendStd=0;
			if(proximal.size()>0)
			for(auto &i:proximal)
				proximalBendStd+=pow(proximalBend-bendMap[i],2);
			if(proximal.size()>0)
				proximalBendStd/=static_cast<double>(proximal.size()-1);
			proximalBendStd=sqrt(proximalBendStd);
			
			//output proximal order parameter to file
			if(proximal.size()>0)
			{
				std::string proximalName("proximalOrder_");
				proximalName+=name;
				proximalName+=".dat";
				std::cerr << "Writing " << proximalName << " at time " << time << "..." << std::endl;
				std::fstream proximalOrder(proximalName.c_str(), std::ios::out | std::ios::app);
				if(proximalOrder.is_open())
					proximalOrder << time << '\t' << proximalBend << '\t' << proximalBendStd << std::endl;
			}
			
			double distalBend=0;
			for(auto &i:distal)
				distalBend+=bendMap[i];
			if(distal.size()>0)
				distalBend/=static_cast<double>(distal.size());
			double distalBendStd=0;
			if(distal.size()>0)
			for(auto &i:distal)
				distalBendStd+=pow(distalBend-bendMap[i],2);
			if(distal.size()>0)
				distalBendStd/=static_cast<double>(distal.size()-1);
			distalBendStd=sqrt(distalBendStd);
			
			//output distal order parameter to file
			if(distal.size()>0)
			{
				std::string distalName("distalOrder_");
				distalName+=name;
				distalName+=".dat";
				std::cerr << "Writing " << distalName << " at time " << time << "..." << std::endl;
				std::fstream distalOrder(distalName.c_str(), std::ios::out | std::ios::app);
				if(distalOrder.is_open())
					distalOrder << time << '\t' << distalBend << '\t' << distalBendStd << std::endl;
			}
			if(proximal.size()>0 && distal.size()>0 && map)
			{
				std::string proximalXyzName("proximal_");
				proximalXyzName+=name;
				proximalXyzName+=".xyz";
				std::cerr << "Writing " << proximalXyzName << " at time " << time << "..." << std::endl;
				std::fstream proximalXyz(proximalXyzName.c_str(), std::ios::out | std::ios::app);
				if(proximalXyz.is_open())
				{
					proximalXyz << proximal.size() << "\ntest\n";
					for(auto &i:proximal)
						proximalXyz << 1 << '\t' << pMap[i].x << '\t' << pMap[i].y << '\t' << pMap[i].z << '\n';
				}
				std::string distalXyzName("distal_");
				distalXyzName+=name;
				distalXyzName+=".xyz";
				std::cerr << "Writing " << distalXyzName << " at time " << time << "..." << std::endl;
				std::fstream distalXyz(distalXyzName.c_str(), std::ios::out | std::ios::app);
				if(distalXyz.is_open())
				{
					distalXyz << distal.size() << "\ntest\n";
					for(auto &i:distal)
						distalXyz << 1 << '\t' << pMap[i].x << '\t' << pMap[i].y << '\t' << pMap[i].z << '\n';
				}
			}
			
			if(shellThickness>0 && proximal90.size()>0 && distal90.size()>0)
			{
				//for averaging over the proximal shell (shellCount[i][shellIndex]=count)
				// where i is the chain offset
				std::vector< std::vector<int> > proximalShellCount(maxChainLength);
				std::vector<double> avgProximalRadius(maxChainLength,0), avgDistalRadius(maxChainLength,0);
				std::vector<double> avgProximalSurfaceDensity(maxChainLength,0), avgDistalSurfaceDensity(maxChainLength,0);
				//go through every particle in the map
				//i is a consistent index between pMap, bendMap, chainIndices, and oMap
				//for(int i=0;i<pMap.size() && pMap.size()==oMap.size() && pMap.size()==bendMap.size();i++)
				//for(const auto& pM:pMap)//no consistent index between pMap, bendMap, and oMap
				for(const auto& i:proximal90)//no consistent index between pMap, bendMap, and oMap
				{
					auto& chain=*(chainIndices.begin()+i);
					//get the distance between nanoparticle and each chain bead
					//chain offset=j
					for(int j=0;j<chain.size();j++)
					{
						//measured distance with NP
						threeVector<double> d;
						d.x=com.x-p[chain[j]].x;
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						d.y=com.y-p[chain[j]].y;
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						d.z=com.z-p[chain[j]].z;
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						double r=d.x*d.x+d.y*d.y+d.z*d.z;
						
						avgProximalRadius[j]+=sqrt(r);
						avgProximalSurfaceDensity[j]+=1;
						
						int shellIndex=floor((sqrt(r))/shellThickness);
						while(shellIndex>=proximalShellCount[j].size())
						{
							for(auto& psc:proximalShellCount)
								psc.push_back(0);
						}
						proximalShellCount[j][shellIndex]++;
					}
				}
				
				std::vector< std::vector<int> > distalShellCount(maxChainLength);
				for(const auto& i:distal90)//no consistent index between pMap, bendMap, and oMap
				{
					auto& chain=*(chainIndices.begin()+i);
					//get the distance between nanoparticle and each chain bead
					//chain offset=j
					for(int j=0;j<chain.size();j++)
					{
						//measured distance with NP
						threeVector<double> d;
						d.x=com.x-p[chain[j]].x;
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						d.y=com.y-p[chain[j]].y;
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						d.z=com.z-p[chain[j]].z;
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						double r=d.x*d.x+d.y*d.y+d.z*d.z;
						
						avgDistalRadius[j]+=sqrt(r);
						avgDistalSurfaceDensity[j]+=1;
						
						int shellIndex=floor((sqrt(r))/shellThickness);
						while(shellIndex>=distalShellCount[j].size())
						{
							for(auto& dsc:distalShellCount)
								dsc.push_back(0);
						}
						distalShellCount[j][shellIndex]++;
					}
				}
				//match the length of the distal shell to the proximal shell length
				while(proximalShellCount[0].size()>distalShellCount[0].size())
					for(auto& dsc:distalShellCount)
						dsc.push_back(0);
					
				//match the length of the proximal shell to the distal shell length
				while(proximalShellCount[0].size()<distalShellCount[0].size())
					for(auto& psc:proximalShellCount)
						psc.push_back(0);
				if(avgShellThickness>0)
				{
					if(profileTimeSpan>avgShellThickness)
					{
						proximalShellProfile.clear();
						distalShellProfile.clear();
						nProfiles=0;
						profileTimeSpan=0;
					}
					if(proximalShellProfile.size()<proximalShellCount.size())
						proximalShellProfile.resize(proximalShellCount.size(),std::vector<double>(proximalShellCount[0].size(),0));
					if(distalShellProfile.size()<distalShellCount.size())
						distalShellProfile.resize(distalShellCount.size(),std::vector<double>(distalShellCount[0].size(),0));
					while(proximalShellProfile[0].size()<proximalShellCount[0].size())
						for(auto& psp:proximalShellProfile)
							psp.push_back(0);
					while(distalShellProfile[0].size()<distalShellCount[0].size())
						for(auto& dsp:distalShellProfile)
							dsp.push_back(0);
				}
				
				//std::vector<std::vector<double> > proximalShellProfile, distalShellProfile;
				//int nProfiles=0;
				//double profileTimeSpan=0;
				nProfiles++;
				for(int i=0;i<proximalShellCount.size();i++)
				{
					std::stringstream shellName;
					shellName << "proximalShell_" << i << "_" << name << ".dat";
					std::cerr << "Writing " << shellName.str() << " at time " << time << "..." << std::endl;
					
					std::fstream shellFile;
					shellFile.open(shellName.str(), std::ios::out | std::ios::app);
					for(int j=0;j<proximalShellCount[i].size() && shellFile.is_open();j++)
					{
						double shellInnerR=static_cast<double>(j)*shellThickness;
						double shellOuterR=static_cast<double>(j+1)*shellThickness;
						//double shellVolume=(4.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0));
						double shellVolume=(2.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0))*(1-0.9);
						//shellVolume*=acos(wrapCostheta)/M_PI;
						shellFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) 
						  << '\t' << static_cast<double>(proximalShellCount[i][j])/shellVolume << std::endl;
						if(j==proximalShellCount[i].size()-1)
							shellFile << std::endl;
						proximalShellProfile[i][j]+=static_cast<double>(proximalShellCount[i][j])/shellVolume;
					}
				}
				if(profileTimeSpan+System.readStoreInterval()>avgShellThickness && avgShellThickness>0)
				for(int i=0;i<proximalShellProfile.size();i++)
				{
					std::stringstream proximalShellProfileName;
					proximalShellProfileName << "avgProximalShell_" << i << "_" << avgShellThickness << "_" << name << ".dat";
					std::cerr << "Writing " << proximalShellProfileName.str() << " at time " << time << "..." << std::endl;
					
					std::fstream proximalShellProfileFile;
					proximalShellProfileFile.open(proximalShellProfileName.str(), std::ios::out | std::ios::app);
					for(int j=0;j<proximalShellProfile[i].size() && proximalShellProfileFile.is_open();j++)
					{
						proximalShellProfileFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) << '\t'
							<< proximalShellProfile[i][j]/static_cast<double>(nProfiles) << std::endl;
						if(j==proximalShellProfile[i].size()-1)
							proximalShellProfileFile << std::endl;
					}
				}
				
				
				for(int i=0;i<distalShellCount.size();i++)
				{
					std::stringstream shellName;
					shellName << "distalShell_" << i << "_" << name << ".dat";
					std::cerr << "Writing " << shellName.str() << " at time " << time << "..." << std::endl;
					
					std::fstream shellFile;
					shellFile.open(shellName.str(), std::ios::out | std::ios::app);
					for(int j=0;j<distalShellCount[i].size() && shellFile.is_open();j++)
					{
						double shellInnerR=static_cast<double>(j)*shellThickness;
						double shellOuterR=static_cast<double>(j+1)*shellThickness;
						//double shellVolume=(4.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0));
						double shellVolume=(2.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0))*(1-0.9);
						//shellVolume*=acos(wrapCostheta)/M_PI;
						shellFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) 
						  << '\t' << static_cast<double>(distalShellCount[i][j])/shellVolume << std::endl;
						if(j==distalShellCount[i].size()-1)
							shellFile << std::endl;
						distalShellProfile[i][j]+=static_cast<double>(distalShellCount[i][j])/shellVolume;
					}
				}
				
				if(profileTimeSpan+System.readStoreInterval()>avgShellThickness && avgShellThickness>0)
				for(int i=0;i<distalShellProfile.size();i++)
				{
					std::stringstream distalShellProfileName;
					distalShellProfileName << "avgDistalShell_" << i << "_" << avgShellThickness << "_" << name << ".dat";
					std::cerr << "Writing " << distalShellProfileName.str() << " at time " << time << "..." << std::endl;
					
					std::fstream distalShellProfileFile;
					distalShellProfileFile.open(distalShellProfileName.str(), std::ios::out | std::ios::app);
					for(int j=0;j<distalShellProfile[i].size() && distalShellProfileFile.is_open();j++)
					{
						distalShellProfileFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) << '\t'
							<< distalShellProfile[i][j]/static_cast<double>(nProfiles) << std::endl;
						if(j==distalShellProfile[i].size()-1)
							distalShellProfileFile << std::endl;
					}
				}
				profileTimeSpan+=System.readStoreInterval();
			
				for(int i=0;i<avgDistalRadius.size();i++)
				{
					//avgDistalSurfaceDensity[i]==distal.size(), which is actually a count here
					if(avgDistalSurfaceDensity[i]>0)
						avgDistalRadius[i]/=avgDistalSurfaceDensity[i];
					//Now it is a density
					if(avgDistalRadius[i]>0)
						avgDistalSurfaceDensity[i]/=avgDistalRadius[i]*avgDistalRadius[i]*2.0*M_PI*(1.0-0.9);
				}
				
				for(int i=0;i<avgProximalRadius.size();i++)
				{
					//avgProximalSurfaceDensity[i]==proximal.size(), which is actually a count here
					if(avgProximalSurfaceDensity[i]>0)
						avgProximalRadius[i]/=avgProximalSurfaceDensity[i];
					//Now it is a density
					if(avgProximalRadius[i]>0)
						avgProximalSurfaceDensity[i]/=avgProximalRadius[i]*avgProximalRadius[i]*2.0*M_PI*(1.0-0.9);
				}
				
				if(avgProximalSurfaceDensity[0]>0 && avgDistalSurfaceDensity[0]>0)
				{
					std::string distalSurfaceDensityName("distalSurfaceDensity_");
					distalSurfaceDensityName+=name;
					distalSurfaceDensityName+=".dat";
					std::cerr << "Writing " << distalSurfaceDensityName << " at time " << time << "..." << std::endl;
					
					std::fstream distalSurfaceDensityFile;
					distalSurfaceDensityFile.open(distalSurfaceDensityName.c_str(), std::ios::out | std::ios::app);
					if(distalSurfaceDensityFile.is_open())
					{
						distalSurfaceDensityFile << time;
						for(auto& dsd:avgDistalSurfaceDensity)
							distalSurfaceDensityFile << '\t' << dsd;
						for(auto& dr:avgDistalRadius)
							distalSurfaceDensityFile << '\t' << dr;
						distalSurfaceDensityFile << '\t' << distal.size() << std::endl;
						distalSurfaceDensityFile.close();
					}
					
					std::string proximalSurfaceDensityName("proximalSurfaceDensity_");
					proximalSurfaceDensityName+=name;
					proximalSurfaceDensityName+=".dat";
					std::cerr << "Writing " << proximalSurfaceDensityName << " at time " << time << "..." << std::endl;
					
					std::fstream proximalSurfaceDensityFile;
					proximalSurfaceDensityFile.open(proximalSurfaceDensityName.c_str(), std::ios::out | std::ios::app);
					if(proximalSurfaceDensityFile.is_open())
					{
						proximalSurfaceDensityFile << time;
						for(auto& psd:avgProximalSurfaceDensity)
							proximalSurfaceDensityFile << '\t' << psd;
						for(auto& pr:avgProximalRadius)
							proximalSurfaceDensityFile << '\t' << pr;
						proximalSurfaceDensityFile << '\t' << proximal.size() << std::endl;
						proximalSurfaceDensityFile.close();
					}
				}
			}
		}
		
		//checking that proximal and distal exist (exclusive or)
		if((proximal.size()>0 && distal.size()==0) || (proximal.size()==0 && distal.size()>0))
		{
			std::cerr << "Warning: A proximal or distal layer exists, but not concurrently!" << std::endl;
			std::cerr << "Result: No shells are output!" << std::endl;
		}
		
		if(normalCutoff>0)
		{
			orientRange oR(pMap,oMap,normalCutoff);
			std::cerr << size.x << ' ' << size.y << ' ' << size.z << std::endl;
			Cell<double> pairInteractions(&pMap[0], pMap.size(), normalCutoff, size);
			pairInteractions.build();
			pairInteractions.twoWayComparePeriodicIndex(oR);
			oR.finalize();
			std::vector<threeVector<double> > nHat=oR.getNHat();
			
			std::string localOrderName("localOrder_");
			localOrderName+=name;
			localOrderName+=".dat";
			std::cerr << "Writing " << localOrderName << " at time " << time << "..." << std::endl;
			
			std::fstream localOrderFile;
			localOrderFile.open(localOrderName.c_str(), std::ios::out | std::ios::app);
			if(localOrderFile.is_open())
			{
				double avgS=0;
				for(int i=0;i<nHat.size() && i<oMap.size();i++)
					avgS+=0.5*(3.0*pow(dotProduct(oMap[i],nHat[i]),2.0)-1);
				if(nHat.size()>0 && oMap.size()>0)
					avgS/=static_cast<double>(nHat.size());
				localOrderFile << time << '\t' << avgS << std::endl;
			}
			
			if(proximal.size()>0)
			{
				std::string localProximalOrderName("localProximalOrder_");
				localProximalOrderName+=name;
				localProximalOrderName+=".dat";
				std::cerr << "Writing " << localProximalOrderName << " at time " << time << "..." << std::endl;
				
				std::fstream localProximalOrderFile;
				localProximalOrderFile.open(localProximalOrderName.c_str(), std::ios::out | std::ios::app);
				if(localProximalOrderFile.is_open())
				{
					double avgS=0;
					for(auto& i:proximal)
						avgS+=0.5*(3.0*pow(dotProduct(oMap[i],nHat[i]),2.0)-1);
					avgS/=static_cast<double>(proximal.size());
					localProximalOrderFile << time << '\t' << avgS << std::endl;
				}
			}
			
			if(distal.size()>0)
			{
				std::string localDistalOrderName("localDistalOrder_");
				localDistalOrderName+=name;
				localDistalOrderName+=".dat";
				std::cerr << "Writing " << localDistalOrderName << " at time " << time << "..." << std::endl;
				
				std::fstream localDistalOrderFile;
				localDistalOrderFile.open(localDistalOrderName.c_str(), std::ios::out | std::ios::app);
				if(localDistalOrderFile.is_open())
				{
					double avgS=0;
					for(auto& i:distal)
						avgS+=0.5*(3.0*pow(dotProduct(oMap[i],nHat[i]),2.0)-1);
					avgS/=static_cast<double>(distal.size());
					localDistalOrderFile << time << '\t' << avgS << std::endl;
				}
			}
		}
		
		//consolidate the proximal and distal leaflets to within cos(theta)<0.9 of zHat
		//go through all chains in proximal leaflet
		if(proximal.size()>0 && !vesicle)
		for(auto &i:proximal)
		{
			//go through each nanoparticle index
			for(auto& npIndex:npIndices)
			{
				//get the distance between nanoparticle and each chain center of mass
				threeVector<double> d;
				d.x=p[npIndex.x].x-pMap[i].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[npIndex.x].y-pMap[i].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[npIndex.x].z-pMap[i].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				double r=d.x*d.x+d.y*d.y+d.z*d.z;
				
				//find the cos(angle) it makes with the z component
				d=unitVector(d);
				double costheta=dotProduct(d,zHat);
				
				//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets)
				/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
				         z      PDDDDDDDDDDDDDDDDDDDDD / {
				npCutoff ^     PD  } wrapDistance
				         |     ^ PD
				         |    /    PD
				npRadius |__ /      PD
				         |  \       PD
				         | / |       PD
				         |/  |<- costheta
				       NP*-------> x or y
				
				*/
				//if it is nearby and intersects the proximal wrapping angle
				if(costheta>0.9)
					proximal90.push_back(i);
			}
		}
		//go through all chains in proximal leaflet
		if(distal.size()>0 && !vesicle)
		for(auto &i:distal)
		{
			//go through each nanoparticle index
			for(auto& npIndex:npIndices)
			{
				//get the distance between nanoparticle and each chain center of mass
				threeVector<double> d;
				d.x=p[npIndex.x].x-pMap[i].x;
				if(d.x>size.x/2.0) d.x-=size.x;
				if(d.x<-size.x/2.0) d.x+=size.x;
				d.y=p[npIndex.x].y-pMap[i].y;
				if(d.y>size.y/2.0) d.y-=size.y;
				if(d.y<-size.y/2.0) d.y+=size.y;
				d.z=p[npIndex.x].z-pMap[i].z;
				if(d.z>size.z/2.0) d.z-=size.z;
				if(d.z<-size.z/2.0) d.z+=size.z;
				double r=d.x*d.x+d.y*d.y+d.z*d.z;
				
				//find the cos(angle) it makes with the z component
				d=unitVector(d);
				double costheta=dotProduct(d,zHat);
				
				//find the cos(angle) it makes with the z component (P=proximal and D=distal leaflets)
				/*               PPPPPPPPPPPPPPPPPPPPP \ <--- Bilayer
				         z      PDDDDDDDDDDDDDDDDDDDDD / {
				npCutoff ^     PD  } wrapDistance
				         |     ^ PD
				         |    /    PD
				npRadius |__ /      PD
				         |  \       PD
				         | / |       PD
				         |/  |<- costheta
				       NP*-------> x or y
				
				*/
				//if it is nearby and intersects the proximal wrapping angle
				if(costheta>0.9)
					distal90.push_back(i);
			}
		}
		proximal=proximal90;
		distal=distal90;
		//Do shells assuming the proximal and distal leaflets are available
		if(shellThickness>0 && proximal.size()>0 && distal.size()>0 && !vesicle)
		{
			//for averaging over the proximal shell (shellCount[i][shellIndex]=count)
			// where i is the chain offset
			std::vector< std::vector<int> > proximalShellCount(maxChainLength);
			
			std::vector<double> avgProximalRadius(maxChainLength,0), avgDistalRadius(maxChainLength,0);
			std::vector<double> avgProximalSurfaceDensity(maxChainLength,0), avgDistalSurfaceDensity(maxChainLength,0);
			
			//go through every particle in the map
			//i is a consistent index between pMap, bendMap, chainIndices, and oMap
			//for(int i=0;i<pMap.size() && pMap.size()==oMap.size() && pMap.size()==bendMap.size();i++)
			//for(const auto& pM:pMap)//no consistent index between pMap, bendMap, and oMap
			for(const auto& i:proximal)//no consistent index between pMap, bendMap, and oMap
			{
				auto& chain=*(chainIndices.begin()+i);
				//go through each nanoparticle index
				for(auto& npIndex:npIndices)
				{
					//get the distance between nanoparticle and each chain bead
					//chain offset=j
					for(int j=0;j<chain.size();j++)
					{
						//measured distance with NP
						threeVector<double> d;
						d.x=p[npIndex.x].x-p[chain[j]].x;
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						d.y=p[npIndex.x].y-p[chain[j]].y;
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						d.z=p[npIndex.x].z-p[chain[j]].z;
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						double r=d.x*d.x+d.y*d.y+d.z*d.z;
						
						avgProximalRadius[j]+=sqrt(r);
						avgProximalSurfaceDensity[j]+=1;
						
						int shellIndex=floor((sqrt(r))/shellThickness);
						while(shellIndex>=proximalShellCount[j].size())
						{
							for(auto& psc:proximalShellCount)
								psc.push_back(0);
						}
						proximalShellCount[j][shellIndex]++;
					}
				}
			}
			
			//for averaging over the distal shell (shellCount[i][shellIndex]=count)
			// where i is the chain offset
			std::vector< std::vector<int> > distalShellCount(maxChainLength);
			
			//go through every particle in the map
			//i is a consistent index between pMap, bendMap, chainIndices, and oMap
			for(const auto& i:distal)//no consistent index between pMap, bendMap, and oMap
			{
				auto& chain=*(chainIndices.begin()+i);
				//go through each nanoparticle index
				for(auto& npIndex:npIndices)
				{
					//get the distance between nanoparticle and each chain bead
					//chain offset=j
					for(int j=0;j<chain.size();j++)
					{
						//measured distance with NP
						threeVector<double> d;
						//std::cerr << chain[j] << std::endl;
						d.x=p[npIndex.x].x-p[chain[j]].x;
						if(d.x>size.x/2.0) d.x-=size.x;
						if(d.x<-size.x/2.0) d.x+=size.x;
						d.y=p[npIndex.x].y-p[chain[j]].y;
						if(d.y>size.y/2.0) d.y-=size.y;
						if(d.y<-size.y/2.0) d.y+=size.y;
						d.z=p[npIndex.x].z-p[chain[j]].z;
						if(d.z>size.z/2.0) d.z-=size.z;
						if(d.z<-size.z/2.0) d.z+=size.z;
						double r=d.x*d.x+d.y*d.y+d.z*d.z;
						
						avgDistalRadius[j]+=sqrt(r);
						avgDistalSurfaceDensity[j]+=1;
						
						int shellIndex=floor((sqrt(r))/shellThickness);
						while(shellIndex>=distalShellCount[j].size())
							for(auto& dsc:distalShellCount)
								dsc.push_back(0);
						distalShellCount[j][shellIndex]++;
					}
				}
			}
			
			//match the length of the distal shell to the proximal shell length
			while(proximalShellCount[0].size()>distalShellCount[0].size())
				for(auto& dsc:distalShellCount)
					dsc.push_back(0);
				
			//match the length of the proximal shell to the distal shell length
			while(proximalShellCount[0].size()<distalShellCount[0].size())
				for(auto& psc:proximalShellCount)
					psc.push_back(0);
			
			for(int i=0;i<avgDistalRadius.size();i++)
			{
				//avgDistalSurfaceDensity[i]==distal.size(), which is actually a count here
				if(avgDistalSurfaceDensity[i]>0)
					avgDistalRadius[i]/=avgDistalSurfaceDensity[i];
				//Now it is a density
				if(avgDistalRadius[i]>0)
					avgDistalSurfaceDensity[i]/=avgDistalRadius[i]*avgDistalRadius[i]*2.0*M_PI*(1.0-0.9);
			}
			
			for(int i=0;i<avgProximalRadius.size();i++)
			{
				//avgProximalSurfaceDensity[i]==proximal.size(), which is actually a count here
				if(avgProximalSurfaceDensity[i]>0)
					avgProximalRadius[i]/=avgProximalSurfaceDensity[i];
				//Now it is a density
				if(avgProximalRadius[i]>0)
					avgProximalSurfaceDensity[i]/=avgProximalRadius[i]*avgProximalRadius[i]*2.0*M_PI*(1.0-0.9);
			}
			
			if(avgProximalSurfaceDensity[0]>0 && avgDistalSurfaceDensity[0]>0)
			{
				std::string distalSurfaceDensityName("distalSurfaceDensity_");
				distalSurfaceDensityName+=name;
				distalSurfaceDensityName+=".dat";
				std::cerr << "Writing " << distalSurfaceDensityName << " at time " << time << "..." << std::endl;
				
				std::fstream distalSurfaceDensityFile;
				distalSurfaceDensityFile.open(distalSurfaceDensityName.c_str(), std::ios::out | std::ios::app);
				if(distalSurfaceDensityFile.is_open())
				{
					distalSurfaceDensityFile << time;
					for(auto& dsd:avgDistalSurfaceDensity)
						distalSurfaceDensityFile << '\t' << dsd;
					for(auto& dr:avgDistalRadius)
						distalSurfaceDensityFile << '\t' << dr;
					distalSurfaceDensityFile << '\t' << distal.size() << std::endl;
					distalSurfaceDensityFile.close();
				}
				
				std::string proximalSurfaceDensityName("proximalSurfaceDensity_");
				proximalSurfaceDensityName+=name;
				proximalSurfaceDensityName+=".dat";
				std::cerr << "Writing " << proximalSurfaceDensityName << " at time " << time << "..." << std::endl;
				
				std::fstream proximalSurfaceDensityFile;
				proximalSurfaceDensityFile.open(proximalSurfaceDensityName.c_str(), std::ios::out | std::ios::app);
				if(proximalSurfaceDensityFile.is_open())
				{
					proximalSurfaceDensityFile << time;
					for(auto& psd:avgProximalSurfaceDensity)
						proximalSurfaceDensityFile << '\t' << psd;
					for(auto& pr:avgProximalRadius)
						proximalSurfaceDensityFile << '\t' << pr;
					proximalSurfaceDensityFile << '\t' << proximal.size() << std::endl;
					proximalSurfaceDensityFile.close();
				}
			}
			
			if(avgShellThickness>0)
			{
				if(profileTimeSpan>avgShellThickness)
				{
					proximalShellProfile.clear();
					distalShellProfile.clear();
					nProfiles=0;
					profileTimeSpan=0;
				}
				if(proximalShellProfile.size()<proximalShellCount.size())
					proximalShellProfile.resize(proximalShellCount.size(),std::vector<double>(proximalShellCount[0].size(),0));
				if(distalShellProfile.size()<distalShellCount.size())
					distalShellProfile.resize(distalShellCount.size(),std::vector<double>(distalShellCount[0].size(),0));
				while(proximalShellProfile[0].size()<proximalShellCount[0].size())
					for(auto& psp:proximalShellProfile)
						psp.push_back(0);
				while(distalShellProfile[0].size()<distalShellCount[0].size())
					for(auto& dsp:distalShellProfile)
						dsp.push_back(0);
			}
			
			//std::vector<std::vector<double> > proximalShellProfile, distalShellProfile;
			//int nProfiles=0;
			//double profileTimeSpan=0;
			nProfiles++;
			for(int i=0;i<proximalShellCount.size();i++)
			{
				std::stringstream shellName;
				shellName << "proximalShell_" << i << "_" << name << ".dat";
				std::cerr << "Writing " << shellName.str() << " at time " << time << "..." << std::endl;
				
				std::fstream shellFile;
				shellFile.open(shellName.str(), std::ios::out | std::ios::app);
				for(int j=0;j<proximalShellCount[i].size() && shellFile.is_open();j++)
				{
					double shellInnerR=static_cast<double>(j)*shellThickness;
					double shellOuterR=static_cast<double>(j+1)*shellThickness;
					//double shellVolume=(4.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0));
					double shellVolume=(2.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0))*(1-0.9);
					//shellVolume*=acos(wrapCostheta)/M_PI;
					shellFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) 
					  << '\t' << static_cast<double>(proximalShellCount[i][j])/shellVolume << std::endl;
					if(j==proximalShellCount[i].size()-1)
						shellFile << std::endl;
					proximalShellProfile[i][j]+=static_cast<double>(proximalShellCount[i][j])/shellVolume;
				}
			}
			if(profileTimeSpan+System.readStoreInterval()>avgShellThickness && avgShellThickness>0)
			for(int i=0;i<proximalShellProfile.size();i++)
			{
				std::stringstream proximalShellProfileName;
				proximalShellProfileName << "avgProximalShell_" << i << "_" << avgShellThickness << "_" << name << ".dat";
				std::cerr << "Writing " << proximalShellProfileName.str() << " at time " << time << "..." << std::endl;
				
				std::fstream proximalShellProfileFile;
				proximalShellProfileFile.open(proximalShellProfileName.str(), std::ios::out | std::ios::app);
				for(int j=0;j<proximalShellProfile[i].size() && proximalShellProfileFile.is_open();j++)
				{
					proximalShellProfileFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) << '\t'
						<< proximalShellProfile[i][j]/static_cast<double>(nProfiles) << std::endl;
					if(j==proximalShellProfile[i].size()-1)
						proximalShellProfileFile << std::endl;
				}
			}
			
			
			for(int i=0;i<distalShellCount.size();i++)
			{
				std::stringstream shellName;
				shellName << "distalShell_" << i << "_" << name << ".dat";
				std::cerr << "Writing " << shellName.str() << " at time " << time << "..." << std::endl;
				
				std::fstream shellFile;
				shellFile.open(shellName.str(), std::ios::out | std::ios::app);
				for(int j=0;j<distalShellCount[i].size() && shellFile.is_open();j++)
				{
					double shellInnerR=static_cast<double>(j)*shellThickness;
					double shellOuterR=static_cast<double>(j+1)*shellThickness;
					//double shellVolume=(4.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0));
					double shellVolume=(2.0/3.0)*M_PI*(pow(shellOuterR,3.0)-pow(shellInnerR,3.0))*(1-0.9);
					//shellVolume*=acos(wrapCostheta)/M_PI;
					shellFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) 
					  << '\t' << static_cast<double>(distalShellCount[i][j])/shellVolume << std::endl;
					if(j==distalShellCount[i].size()-1)
						shellFile << std::endl;
					distalShellProfile[i][j]+=static_cast<double>(distalShellCount[i][j])/shellVolume;
				}
			}
			
			if(profileTimeSpan+System.readStoreInterval()>avgShellThickness && avgShellThickness>0)
			for(int i=0;i<distalShellProfile.size();i++)
			{
				std::stringstream distalShellProfileName;
				distalShellProfileName << "avgDistalShell_" << i << "_" << avgShellThickness << "_" << name << ".dat";
				std::cerr << "Writing " << distalShellProfileName.str() << " at time " << time << "..." << std::endl;
				
				std::fstream distalShellProfileFile;
				distalShellProfileFile.open(distalShellProfileName.str(), std::ios::out | std::ios::app);
				for(int j=0;j<distalShellProfile[i].size() && distalShellProfileFile.is_open();j++)
				{
					distalShellProfileFile << static_cast<double>(j)*shellThickness+(shellThickness/2.0) << '\t'
						<< distalShellProfile[i][j]/static_cast<double>(nProfiles) << std::endl;
					if(j==distalShellProfile[i].size()-1)
						distalShellProfileFile << std::endl;
				}
			}
			profileTimeSpan+=System.readStoreInterval();
		}
		time+=System.readStoreInterval();
		//we should use size file to determine what the start time is
		if(frames)
		{
			if(sizeFile.is_open())
			{
				double sTime=0;
				threeVector<double> sCurrent;
				while(sTime<time-System.readMeasureInterval()+System.readDeltaT() && !sizeFile.eof())
					sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				size=sCurrent;
				time=sTime;
			std::cerr << "Got (" << size.x << ", " << size.y << ", " << size.z << ") at time " << time << std::endl;
			}
		}
		frameIndex++;
	}
	
	return 0;
}

