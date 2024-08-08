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

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		std::cerr << argv[0] << " name fromTime molIndex type cutoff \n";
		std::cerr << "Output columns:" << std::endl;
		std::cerr << "1: Average number of contacts per polymer" << std::endl;
		std::cerr << "2: Average length of contacts per polymer" << std::endl;
		std::cerr << "3,4,5: x y and z radius of gyration of polymers not in contact" << std::endl;
		std::cerr << "6,7,8: x y and z radius of gyration of strands not in contact for polymers in contact" << std::endl;
		std::cerr << "9,10,11: x y and z radius of gyration of strands in contact of polymers in contact" << std::endl;
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
	
	double avgNContact=0,nContact=0,nnContact=0;
	double lContact=0;
	threeVector<double> avgCROG=0,avgNROG=0,avgNNROG=0;
	//std::vector<double> cHist(chains[0].size(),0.0);
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
			{
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				//std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
			size=sCurrent;
		}
		if(time>=fromTime)
		{
			//make a list of positions
			std::vector<position<double>> tI;
			for(auto &i:tIndices)
				tI.push_back(p[i]);
			threeVector<int> nCells;
			nCells.x=size.x/cutoff;
			nCells.y=size.y/cutoff;
			nCells.z=size.z/cutoff;
			threeVector<double> cSize;
			cSize.x=size.x/nCells.x;
			cSize.y=size.y/nCells.y;
			cSize.z=size.z/nCells.z;
			umapfwd<position<double>> cMap=createCMap(&(*tI.begin()), &(*tI.end()),nCells,cSize);
			
			//initialize contact and difference counts
			std::vector<double> counts(chains[0].size(),0.0);
			std::vector<std::vector<int>> contacts,noContacts;
			//go through all chains
			for(auto chain:chains)
			{
				std::vector<int> contact;
				//index along chain
				int index=0;
				int lastContact=-1;
				bool isContact=false;
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
								isContact=true;
								counts[index]++;
								contact.push_back(index);
							}
						}
					}
					//advance to next chain position
					index++;
				}
				auto it=std::unique(contact.begin(),contact.end());
				contact.resize(std::distance(contact.begin(),it));
				contacts.push_back(contact);
				if(!isContact)
				{
					std::vector<int> noContact;
					for(auto i:chain)
						noContact.push_back(i);
					noContacts.push_back(noContact);
				}
			}
			int cIndex=0;
			for(auto contact:contacts)
			{
				if(contact.size()>1)
				{
					avgNContact+=contact.size();
					lContact+=contact.back()-contact.front();
					std::vector<int> n,c;
					//The zeroth index to the first contact portion of the polymer
					for(int i=0;i<contact.front();i++)
						n.push_back(chains[cIndex][i]);
					//The first contact portion to the end of the polymer.
					//This assumes the substrate is connected to the first monomer
					for(int i=contact.front();i<chains[0].size();i++)
						c.push_back(chains[cIndex][i]);
					
					avgNROG+=radiusOfGyrationAxisIndex(n.begin(),n.end(),p,size);
					avgCROG+=radiusOfGyrationAxisIndex(c.begin(),c.end(),p,size);
					nContact++;
				}
					
				cIndex++;
			}
			for(auto noContact:noContacts)
			{
				avgNNROG+=radiusOfGyrationAxisIndex(noContact.begin(),noContact.end(),p,size);
				nnContact++;
			}
			//get counts and chain position probability
			double totalCount=0;
			for(auto &c:counts)
				totalCount+=c;
			frames++;
			std::cerr << time << ' ' << totalCount << std::endl;
		}
		time+=System.readStoreInterval();
	}
	if(nContact!=0)
	{
		lContact/=nContact;
		avgNContact/=(frames*chains.size());
		//avgNContact/=nContact;
		avgNROG/=nContact;
		avgCROG/=nContact;
	}
	if(nnContact!=0)
		avgNNROG/=nnContact;
	
	std::cout << avgNContact << ' ' << lContact << ' ' <<
		avgNNROG.x << ' ' << avgNNROG.y << ' ' << avgNNROG.z << ' ' << 
		avgNROG.x << ' ' << avgNROG.y << ' ' << avgNROG.z << ' ' << 
		avgCROG.x << ' ' << avgCROG.y << ' ' << avgCROG.z << std::endl;
	
	return 0;
}

