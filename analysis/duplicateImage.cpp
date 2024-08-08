//This is a very short and simple program that shows you how to modify a variable already present in an mpd script.
//The reason this is used is because some variables (forces, potentials, geometry, etc...) are not easily modified
// in the mpd script. Using a file that opens the mpd script, modifies, and then saves it (as below) is far easier
// than modifying the actual variable. Most of it looks like the run.cpp file, it just replaces the "run loop"
// with the modification.

#define NO_REQUIRED_COMMANDS

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

int main(int argc, char* argv[])
{
	if(argc<4)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "usage: " << argv[0] << " name nDup dim\n";
		
		std::cerr << "usage: " << argv[0] << " name nDup dim mol0 ... molN\n";
		
		std::cerr << "\tNote: duplicates along chosen axii, where dim=0,1,2 cooresponding to X,Y,Z\n";
		
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	//parse the command arguments
	std::stringstream cmdArg;
	cmdArg << argv[2] << ' ' << argv[3] << ' ';
	int nDup, dim;
	cmdArg >> nDup >> dim;
	std::vector<int> mol;
	if(argc>4)
	{
		for(int i=4;i<argc;i++)
		{
			cmdArg << argv[i] << ' ';
			int molIndex;
			cmdArg >> molIndex;
			std::cerr << molIndex << '\t' << argv[i] << std::endl;
			mol.push_back(molIndex);
		}
	}
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//create a new file name
	std::stringstream createName;
	createName << argv[1] << "_" << nDup << "Dup";
	std::string newName;
	createName >> newName;
	
	//what to duplicate
	std::vector< molecule<double,fourVector<int> > > m;
	int nMol=System.readNMolecules();
	if(mol.size()==0)//duplicate them all
	{
		molecule<double,fourVector<int> > *mP=System.getMolecule();
		for(int i=0;i<nMol;i++)
			m.push_back(mP[i]);
	}
	else//just duplicate these
	{
		molecule<double,fourVector<int> > *mP=System.getMolecule();
		for(int i=0;i<mol.size();i++)
			m.push_back(mP[mol[i]]);
	}
	
	threeVector<double> size=System.readSize();
	threeVector<double> oldSize=size;
	
	//these are the new relative positions.
	//when the new system size offset is applied, the positions are corrected
	std::vector< position<double> > p;
	std::vector< threeVector<double> > v,a;
	
	for(int molIndex=0;molIndex<m.size();molIndex++)
	{
		position<double> *pPointer=System.getPositions();
		threeVector<double> *vPointer=System.getVelocities();
		threeVector<double> *aPointer=System.getAccelerations();
		
		//pick a structure by type
		switch(m[molIndex].readType())
		{
			case BOND:
			{
				std::map<int,std::vector<int> > bonded;
				//std::vector<int> bInd;
				for(int l=0;l<m[molIndex].readNBond();l++)
				{
					//These are the first and second particles of the bond
					int firstParticle=m[molIndex].getBonds()[l].s[0];
					int secondParticle=m[molIndex].getBonds()[l].s[1];
					//make sure chains are fixed along boundary
					if(pPointer[firstParticle].s[dim]-pPointer[secondParticle].s[dim]>size.s[dim]/2.0)
						pPointer[firstParticle].s[dim]-=size.s[dim];
					if(pPointer[firstParticle].s[dim]-pPointer[secondParticle].s[dim]<-size.s[dim]/2.0)
						pPointer[firstParticle].s[dim]+=size.s[dim];
					bonded[firstParticle].push_back(l);
					bonded[secondParticle].push_back(l);
					//bInd.push_back(first);
					//bInd.push_back(second);
				}
				std::map<int,std::vector<int> >::iterator cMol=bonded.begin();
				for(int newOffset=0;cMol!=bonded.end();newOffset++,cMol++)
				{
					//these are unique
					p.push_back(pPointer[cMol->first]);
					pPointer[cMol->first].type++;
					//(*p.end()).type++;
					v.push_back(vPointer[cMol->first]);
					a.push_back(aPointer[cMol->first]);
					
					//reindexed bonds
					for(int l=0;l<cMol->second.size();l++)
					{
						if(cMol->first==m[molIndex].getBonds()[cMol->second[l]].s[0])
							m[molIndex].getBonds()[cMol->second[l]].s[0]=newOffset;
						else if(cMol->first==m[molIndex].getBonds()[cMol->second[l]].s[1])
							m[molIndex].getBonds()[cMol->second[l]].s[1]=newOffset;
						else
							std::cerr << "No bond found!";
					}
				}
				std::vector<int> stack;
				stack.push_back(m[molIndex].getBonds()[0].s[0]);
				break;
			}
			case BEND:
			{
				std::map<int,std::vector<int> > bonded;
				//std::vector<int> bInd;
				for(int l=0;l<m[molIndex].readNBond();l++)
				{
					//These are the first and second particles of the bond
					int firstParticle=m[molIndex].getBonds()[l].s[0];
					int secondParticle=m[molIndex].getBonds()[l].s[1];
					int thirdParticle=m[molIndex].getBonds()[l].s[2];
					//make sure chains are fixed along boundary
					if(pPointer[firstParticle].s[dim]-pPointer[secondParticle].s[dim]>size.s[dim]/2.0)
						pPointer[firstParticle].s[dim]-=size.s[dim];
					if(pPointer[firstParticle].s[dim]-pPointer[secondParticle].s[dim]<-size.s[dim]/2.0)
						pPointer[firstParticle].s[dim]+=size.s[dim];
					if(pPointer[thirdParticle].s[dim]-pPointer[secondParticle].s[dim]>size.s[dim]/2.0)
						pPointer[thirdParticle].s[dim]-=size.s[dim];
					if(pPointer[thirdParticle].s[dim]-pPointer[secondParticle].s[dim]<-size.s[dim]/2.0)
						pPointer[thirdParticle].s[dim]+=size.s[dim];
					bonded[firstParticle].push_back(l);
					bonded[secondParticle].push_back(l);
					bonded[thirdParticle].push_back(l);
					//bInd.push_back(first);
					//bInd.push_back(second);
				}
				std::map<int,std::vector<int> >::iterator cMol=bonded.begin();
				int oldOffset=p.size();
				for(int newOffset=oldOffset;cMol!=bonded.end();newOffset++,cMol++)
				{
					//these are unique
					p.push_back(pPointer[cMol->first]);
					v.push_back(vPointer[cMol->first]);
					a.push_back(aPointer[cMol->first]);
					
					//reindexed bonds
					for(int l=0;l<cMol->second.size();l++)
					{
						if(cMol->first==m[molIndex].getBonds()[cMol->second[l]].s[0])
							m[molIndex].getBonds()[cMol->second[l]].s[0]=newOffset;
						else if(cMol->first==m[molIndex].getBonds()[cMol->second[l]].s[1])
							m[molIndex].getBonds()[cMol->second[l]].s[1]=newOffset;
						else if(cMol->first==m[molIndex].getBonds()[cMol->second[l]].s[2])
							m[molIndex].getBonds()[cMol->second[l]].s[2]=newOffset;
						else
							std::cerr << "No bond found!";
					}
				}
				break;
			}
			case CHAIN:
			{
				for(int l=0;l<m[molIndex].readNBond();l++)
				{
					int start=m[molIndex].getBonds()[l].s[START];
					int length=m[molIndex].getBonds()[l].s[CHAINLENGTH];
					int nChains=m[molIndex].getBonds()[l].s[NCHAINS];
					//reindexed bonds
					m[molIndex].getBonds()[l].s[START]=p.size();
					
					//new offset
					for(int j=start;j<start+length*nChains;j++)
					{
						int chainElement=(j-start)%length;
						//make sure chains are fixed along boundary
						if(chainElement>0)
						{
							if(pPointer[j].s[dim]-pPointer[j-1].s[dim]>size.s[dim]/2.0)
								pPointer[j].s[dim]-=size.s[dim];
							if(pPointer[j].s[dim]-pPointer[j-1].s[dim]<-size.s[dim]/2.0)
								pPointer[j].s[dim]+=size.s[dim];
						}
						p.push_back(pPointer[j]);
						v.push_back(vPointer[j]);
						a.push_back(aPointer[j]);
					}
				}
				break;
			}
			case BEAD:
			{
				for(int l=0;l<m[molIndex].readNBond();l++)
				{
					//These are the first and second particles of the bond
					int firstParticle=m[molIndex].getBonds()[l].s[0];
					//this one isn't bound, so we don't fix the boundaries
					//reindexed bonds
					m[molIndex].getBonds()[l].s[0]=p.size();
					//new offset
					p.push_back(pPointer[firstParticle]);
					v.push_back(vPointer[firstParticle]);
					a.push_back(aPointer[firstParticle]);
				}
			}
			default:
			{
				//does nothing
				break;
			}
		}
	}
	std::cerr << "Number to duplicate: " << p.size() << std::endl;
	
	//std::cout << p.size() << "\ntest\n";
	//for(int i=0;i<p.size();i++)
	//	std::cout << p[i].type << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << std::endl;
	
	/****************************************
	 * duplication process (along x, dim=0):
	 * 
	 *    x
	 *    _ 
	 * y |_|
	 *  
	 * first iteration:
	 *  _ _
	 * |_|_|
	 *  
	 * second iteration:
	 *  _ _ _
	 * |_|_|_|
	 * 
	 * etc ...
	 * 
	 * **************************************/
	
	
	size.s[dim]*=static_cast<double>(nDup+1);
	
	System.setSize(size);
	
	//go back through to make sure the current particles are within boundaries
	for(int i=0;i<System.readNParticles();i++)
		System.setParticle(i,System.getPositions()[i], System.getVelocities()[i], System.getAccelerations()[i]);
	//this makes increments quick and contiguous
	System.allocParticle(System.readNParticles()+p.size()*nDup);
	
	for(int dup=0;dup<nDup;dup++)
	{
		int offset=System.readNParticles();
		//duplicate particles
		for(int i=0;i<p.size();i++)
		{
			position<double> newP=p[i];
			newP.s[dim]+=static_cast<double>(dup+1)*oldSize.s[dim];
			System.addParticle(newP, v[i], a[i]);
		}
		
		//duplicate molecules
		for(int molIndex=0;molIndex<m.size();molIndex++)
		{
			//this is the new duplicate, all the bonds are zero offset with std::vector p
			molecule<double,fourVector<int> > newMol=m[molIndex];
			
			//pick a structure by type
			switch(m[molIndex].readType())
			{
				case BOND:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//These are the first and second particles of the bond
						newMol.getBonds()[l].s[0]+=offset;
						newMol.getBonds()[l].s[1]+=offset;
					}
				}
				break;
				
				case BEND:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//These are the first and second particles of the bond
						newMol.getBonds()[l].s[0]+=offset;
						newMol.getBonds()[l].s[1]+=offset;
						newMol.getBonds()[l].s[2]+=offset;
					}
					break;
				}
				case CHAIN:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//reindexed bonds
						newMol.getBonds()[l].s[START]+=offset;
					}
					break;
				}
				case BEAD:
				{
					for(int l=0;l<m[molIndex].readNBond();l++)
					{
						//reindexed bonds
						newMol.getBonds()[l].s[0]+=offset;
					}
				}
				default:
				{
					//does nothing
					break;
				}
			}
			//put new molecule in system
			System.addMolecule(newMol);
		}
	}
	//Save it
	fileIO.open(newName.c_str(),std::ios::out);
	fileIO.write();
	fileIO.close();
	
	return 0;
}

