/** \brief Locates blebs and necks. Outputs the regions (in an .xyz file) and the sizes (in .dat files).
 * 
 * 
 */

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

class lipidCluster {
	public:
		lipidCluster(std::vector< position<double> > &lipids, double cutoff)
		{
			lip=lipids;
			std::cout << lip[0].type << '\t' << lip[0].x << '\t' << lip[0].y << '\t' << lip[0].z << '\n' << cutoff << '\n';
			cutoffSqr=cutoff*cutoff;
		};
		
		inline int nElements()
		{
			return lip.size();
		};
		
		inline bool compare(int &i, int &j)
		{
			bool result=false;
			threeVector<double> d;
			d.x=lip[i].x-lip[j].x;
			d.y=lip[i].y-lip[j].y;
			d.z=lip[i].z-lip[j].z;
			
			if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
			{
				result=true;
			}
			
			return result;
		};
		
		int next(int i);
	private:
		std::vector< position<double> > lip;
		double cutoffSqr;
};

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Output molecule residue data:\n";
		std::cout << "\tUsage: " << argv[0] << " name\n";
		
		std::cout << "Output bleb and neck regional and size data:\n";
		std::cout << "\tUsage: " << argv[0] << " name residue0 residue1 residue2 ...\n";
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
	
	position<double> *p=System.getPositions();
	
	if(argc>2)
	{
		//get the residues
		std::vector<int> residue;
		std::stringstream cmdArg;
		for(int i=2;i<argc;i++)
			cmdArg << argv[i] << ' ';
		
		while(!cmdArg.eof())
		{
			int buf;
			cmdArg >> buf;
			if(!cmdArg.eof())
			{
				residue.push_back(buf);
			}
		}
		
		//gather lipid center of masses
		std::vector< position<double> > lipCOM;
		
		for(int i=0;i<residue.size();i++)
		{
			//shorthand
			molecule<double, fourVector<int> > *m=System.getMolecule();
			
			//combine center of masses by type
			switch(m[residue[i]].readType())
			{
				case BOND:
				{
					//initialize our center of mass reduction to zero
					position<double> comVector;
					comVector.x=0;
					comVector.y=0;
					comVector.z=0;
					comVector.type=residue[i];
					
					for(int l=0;l<m[residue[i]].readNBond();l++)
					{
						//These are the first and second particles of the bond
						int firstParticle=m[residue[i]].getBonds()[l].s[0];
						
						int secondParticle=m[residue[i]].getBonds()[l].s[1];
						
						comVector.x+=p[firstParticle].x;
						comVector.y+=p[firstParticle].y;
						comVector.z+=p[firstParticle].z;
						
						comVector.x+=p[secondParticle].x;
						comVector.y+=p[secondParticle].y;
						comVector.z+=p[secondParticle].z;
						
					}
					comVector.x/=static_cast<double>(m[residue[i]].readNBond()*2);
					comVector.y/=static_cast<double>(m[residue[i]].readNBond()*2);
					comVector.z/=static_cast<double>(m[residue[i]].readNBond()*2);
					
					lipCOM.push_back(comVector);
					break;
				}
				case BEND:
				{
					//initialize our center of mass reduction to zero
					position<double> comVector;
					comVector.x=0;
					comVector.y=0;
					comVector.z=0;
					comVector.type=residue[i];
					
					for(int l=0;l<m[residue[i]].readNBond();l++)
					{
						//These are the first, second, and third particles of the bend
						int firstParticle=m[residue[i]].getBonds()[l].s[0];
						int secondParticle=m[residue[i]].getBonds()[l].s[1];
						int thirdParticle=m[residue[i]].getBonds()[l].s[2];
						
						comVector.x+=p[firstParticle].x;
						comVector.y+=p[firstParticle].y;
						comVector.z+=p[firstParticle].z;
						
						comVector.x+=p[secondParticle].x;
						comVector.y+=p[secondParticle].y;
						comVector.z+=p[secondParticle].z;
						
						comVector.x+=p[thirdParticle].x;
						comVector.y+=p[thirdParticle].y;
						comVector.z+=p[thirdParticle].z;
					}
					comVector.x/=static_cast<double>(m[residue[i]].readNBond()*3);
					comVector.y/=static_cast<double>(m[residue[i]].readNBond()*3);
					comVector.z/=static_cast<double>(m[residue[i]].readNBond()*3);
					
					lipCOM.push_back(comVector);
					break;
				}
				case CHAIN:
				{
					for(int l=0;l<m[residue[i]].readNBond();l++)
					{
						int start=m[residue[i]].getBonds()[l].s[START];
						int length=m[residue[i]].getBonds()[l].s[CHAINLENGTH];
						int nChains=m[residue[i]].getBonds()[l].s[NCHAINS];
						for(int j=start;j<start+length*nChains;j+=length)
						{
							//initialize our center of mass reduction to zero
							position<double> comVector;
							comVector.x=0;
							comVector.y=0;
							comVector.z=0;
							comVector.type=residue[i];
							
							for(int k=j;k<j+length;k++)
							{
								comVector.x+=p[k].x;
								comVector.y+=p[k].y;
								comVector.z+=p[k].z;
							}
							
							comVector.x/=static_cast<double>(length);
							comVector.y/=static_cast<double>(length);
							comVector.z/=static_cast<double>(length);
							
							lipCOM.push_back(comVector);
						}
					}
					break;
				}
				default:
				{
					//does nothing
					break;
				}
			}
		}
		
		std::fstream xyzFile;
		xyzFile.open("com.xyz", std::ios::out);
		if(xyzFile.is_open())
		{
			xyzFile << lipCOM.size() << "\nasdf\n";
			for(int i=0;i<lipCOM.size();i++)
				xyzFile << lipCOM[i].type << '\t' << lipCOM[i].x << '\t' << lipCOM[i].y << '\t' << lipCOM[i].z << '\n';
			xyzFile.close();
		}
		
		lipidCluster superSet(lipCOM, System.readCutoff()*500);
		std::vector< std::vector<int> > setGroups=cluster(superSet);
		
		std::cout << setGroups.size() << std::endl;
		
	}
	if(argc==2)
	{
		for(int i=0;i<System.readNMolecules();i++)
		{
			//shorthand
			molecule<double, fourVector<int> > *m=System.getMolecule();
			
			switch(m[i].readType())
			{
				case BOND:
					std::cout << "Residue " << i << ":\n";
					std::cout << "\tnBonds: " << m[i].readNBond() << '\n';
					break;
				case BEND:
					std::cout << "Residue " << i << ":\n";
					std::cout << "\tnBends: " << m[i].readNBond() << '\n';
					break;
				case CHAIN:
					std::cout << "Residue " << i << ":\n";
					std::cout << "\tnSubResidues: " << m[i].readNBond() << '\n';
					for(int j=0;j<m[i].readNBond();j++)
					{
						int start=m[i].getBonds()[j].s[START];
						int length=m[i].getBonds()[j].s[CHAINLENGTH];
						int nChains=m[i].getBonds()[j].s[NCHAINS];
						std::cout << "\t\tstart offset: " << start << '\n';
						std::cout << "\t\tchainLength: " << length << '\n';
						std::cout << "\t\tnChains: " << nChains << '\n';
					}
					break;
				default:
					break;
			}
		}
	}
	
	return 0;
}
