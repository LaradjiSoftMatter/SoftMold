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

//For nearest neighbor interactions
#include "cellInclude/cell.h"

//For center of mass and radius of gyration
#include "radiusOfGyrationInclude/radiusOfGyration.h"

#include <set>

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		std::cerr << argv[0] << " name fromTime molIndex type cutoff \n";
		std::cerr << "Output Columns:\n";
		std::cerr << "\t0: Index for counts/Number for histograms.\n";
		std::cerr << "\t1: Histogram of lipid contacts.\n";
		std::cerr << "\t2: Number of lipid contacts per polymer.\n";
		std::cerr << "\t3: Histogram of monomer contacts.\n";
		std::cerr << "\t4: Number of monomers in contact per polymer.\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime=0;
	threeVector<int> n;
	cmdArg >> fromTime;
	
	int molIndex,type;
	double cutoff;
	cmdArg >> molIndex >> type >> cutoff;
	double cutoffSqr=cutoff*cutoff;
	
	/*
	std::vector<int> types;
	
	for(int i=4;i<argc;i++)
	{
		int type;
		cmdArg >> type;
		types.push_back(type);
	}*/
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//get the chains
	molecule<double,fourVector<int>> *m=System.getMolecule();
	int nMol=System.readNMolecules();
	if(molIndex<0 || molIndex>=nMol)
	{
		std::cerr << "molIndex " << molIndex << " is out of range [0," << nMol <<")!" << std::endl;
		return -1;
	}
	if(m[molIndex].readType()!=CHAIN)
	{
		std::cerr << "Molecule " << molIndex << " is not CHAIN type! (" 
			<< m[molIndex].readType() << "!=" << CHAIN << ")" << std::endl;
		return -1;
	}
	if(type<0 || type>=System.readNTypes())
	{
		std::cerr << "type " << type << " is out of range [0," << System.readNTypes() <<")!" << std::endl;
		return -1;
	}
	std::vector<std::vector<int>> chains;
	for(int bond=0;bond<m[molIndex].readNBond();bond++)
	{
		int start=m[molIndex].getBonds()[bond].s[START];
		int nChains=m[molIndex].getBonds()[bond].s[NCHAINS];
		int length=m[molIndex].getBonds()[bond].s[CHAINLENGTH];
		for(int i=start;i<nChains*length;i+=length)
		{
			std::vector<int> chain;
			for(int j=i;j<i+length;j++)
				chain.push_back(j);
			chains.push_back(chain);
		}
		std::cerr << "Found chain " << molIndex << " starting at " << start << " with " << nChains << " chains of length " << length << ".\n";
	}
	
	int nParticles=System.readNParticles();
	position<double> *p=System.getPositions();
	std::vector<int> tIndices;
	for(int i=0;i<nParticles;i++)
		if(p[i].type==type)
			tIndices.push_back(i);
	
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
	
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	int frames=0;
	
	
	//histogram of lipid contacts
	std::vector<double> cHist(chains[0].size(),0.0);
	cHist.reserve(chains[0].size()*4);
	//number of lipid contact by monomer index
	std::vector<double> counts(chains[0].size(),0.0);
	//histogram of monomer contacts
	std::vector<double> mHist(chains[0].size(),0.0);
	//number of monomers in contact per polymer by monomer index
	std::vector<double> mCounts(chains[0].size(),0.0);
	
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
			//make a list of positions
			std::vector<position<double>> tI;
			for(auto &i:tIndices)
			{
				//this changes the type to an index for uniqueness
				p[i].type=i;
				tI.push_back(p[i]);
			}
			threeVector<int> nCells;
			nCells.x=size.x/cutoff;
			nCells.y=size.y/cutoff;
			nCells.z=size.z/cutoff;
			threeVector<double> cSize;
			cSize.x=size.x/nCells.x;
			cSize.y=size.y/nCells.y;
			cSize.z=size.z/nCells.z;
			umapfwd<position<double>> cMap=createCMap(&(*tI.begin()), &(*tI.end()),nCells,cSize);
			
			//initialize contact
			std::vector<std::vector<int>> contacts;
			int totalCount=0;
			
			//go through all chains
			for(auto chain:chains)
			{
				std::set<int> mContact;
				std::vector<int> contact;
				//index along chain
				int index=0;
				int lastContact=-1;
				//all indices in chain
				for(auto i:chain)
				{
					//get the neighbors
					threeVector<int> cellCoord=getCell(p[i],cSize);
					std::vector<int> nIndices=neighIndices(hash(cellCoord,nCells),nCells);
					for(auto neigh:nIndices)
					{
						//for all particles in cell map
						for(auto &j:cMap[neigh])
						{
							//check the distance, accounting for periodic boundaries
							threeVector<double> d;
							d.x=p[i].x-j.x;
							d.y=p[i].y-j.y;
							d.z=p[i].z-j.z;
							if(d.x>size.x/2.0) d.x-=size.x;
							if(d.y>size.y/2.0) d.y-=size.y;
							if(d.z>size.z/2.0) d.z-=size.z;
							if(d.x<-size.x/2.0) d.x+=size.x;
							if(d.y<-size.y/2.0) d.y+=size.y;
							if(d.z<-size.z/2.0) d.z+=size.z;
							double dr=d.x*d.x+d.y*d.y+d.z*d.z;
							if(dr<cutoffSqr)
							{
								//this type is really the index of j
								mContact.insert(j.type);
								counts[index]++;
								contact.push_back(index);
								totalCount++;
							}
						}
					}
					//advance to next chain position
					index++;
				}
				auto it=std::unique(contact.begin(),contact.end());
				contact.resize(std::distance(contact.begin(),it));
				contacts.push_back(contact);
				while(cHist.size()<mContact.size())
					cHist.push_back(0);
				cHist[mContact.size()]++;
			}
			//double cHistTotal=0;
			//std::vector<double> cDiff(chains[0].size(),0.0);
			for(auto contact:contacts)
			{
				//if(contact.size()>40)
				//{
				//	for(auto &m:contact)
				//		std::cerr << m << ' ';
				//	std::cout << std::endl;
				//}	
				mHist[contact.size()]++;
				for(auto &m:contact)
					mCounts[m]++;
			}
			
			//for(int index=0;index<counts.size();index++)
			//	if(totalCount!=0)
			//		std::cout << index << ' ' << counts[index]/totalCount << ' ' << cHist[index] << '\n';
			//	else
			//		std::cout << index << ' ' << 0 << ' ' << 0 << '\n';
			//std::cout << std::endl;
			frames++;
			std::cerr << time << ' ' << totalCount << std::endl;
		}
		time+=System.readStoreInterval();
	}
	double total=chains.size()*frames;
	for(int index=0;index<cHist.size();index++)
	{
		std::cout << index << ' ';
		//if(totalCHist!=0)
		if(frames!=0)
			std::cout << cHist[index]/total << ' ';
		else
			std::cout << 0 << ' ';
		//if(totalCount!=0)
		if(frames!=0 && index<counts.size())
			std::cout << counts[index]/total << ' ';
		else
			std::cout << 0 << ' ';
		//if(totalMHist!=0)
		if(frames!=0 && index<mHist.size())
			std::cout << mHist[index]/total << ' ';
		else
			std::cout << 0 << ' ';
		if(frames!=0 && index<mCounts.size())
			std::cout << mCounts[index]/total << '\n';
		else
			std::cout << 0 << '\n';
	}
	
	return 0;
}

