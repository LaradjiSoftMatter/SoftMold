#include "functions.h"

#ifndef MD_BONDEDMOLECULES
#define MD_BONDEDMOLECULES

//molecule types
#define CHAINOLD 0
#define BENDOLD 1
#define ALLPAIRS 2
#define TORSION 3
#define BONDOLD 4
#define DIHEDRAL 5
//these are extensions (they do what the above do, but better):
#define BOND 6
#define BEND 7
#define CHAIN 8
//different types
#define BEAD 9
#define SOLID 10
#define BOUNDARY 11
#define OFFSET_BOUNDARY 14
#define FLOATING_BASE 15
#define ZTORQUE 16
#define ZPOWERPOTENTIAL 17
//Same as BEAD, but lacks BEAD-BEAD interaction, just non-bonded
#define NANOCORE 18
//Same as harmonic, but includes cutoff
#define BALL 19

//These could actually be used for a user interface, as they induce a 
// rigid axis rotation and pull. RIGIDBEND is not conservative.
// PULLBEAD is actually just a harmonic force applied across a large number of
// particles.
#define RIGIDBEND 12
#define PULLBEAD 13

//for CHAIN type, this is block data, mol[i].bond[whichChain].s[INFO]
#define START 0
#define NCHAINS 1
#define CHAINLENGTH 2

#define NZBOND 1

/**
 * \brief Molecular object for bonding and bonding data.
 * This is a molecule object that can be used as a container for molecules.
 * It has numerous useful functions:
 * 	accepts simple datatypes for bond structures
 * 	add bonds to heap
 * 	set bonds in heap
 * 	add constants to heap
 * 	set constants in heap
 * 	find bonds
 * It needs some useful functions:
 * 	add bonds unless present
 * 	delete bonds
 * 	delete constants
 * 	sort bonds for rapid location
 * It has some issues:
 * 	when sending it too a function:
 * 		molecule mol;
 * 		someFunc(mol);
 * 		definition should be: returnType someFunc(molecule &mol)
 * 		definition will break if: returnType someFunc(molecule mol)
 * 		because allocation/deallocation occurs within the molecule class
 * It should be called a molecular residue structure since it's just a collection
 * of monomers.
 * 
 * Usually a system "blob" would accumulate bonds like this:
 * 	molecule\<double,fourVector\<int\> \> m;
 * 	system.addMolecule(m);
*/
template <typename T, typename bondData>
class molecule {
	public:
		/** \brief Default constructor for 0 or NULL values.
		 * Use it for constructing empty objects like this:
		 * molecule\<double,fourVector\<int\> \> m; 
		 */
		molecule()
		{
			type=-1;
		};
		
		/** \brief Destructor for any memory allocated
		 * Called when the object is destroyed.
		 */
		~molecule()
		{};
		
		
		///Operator to transfer structure information. Copies bonds and bond information, but not allocation information.
		molecule &operator= (molecule<T,bondData> &mol)
		{
			this->setType(mol.readType());
			this->setBonds(mol.readNBond(), mol.getBonds());
			this->setConstants(mol.readNConstant(), mol.getConstants());
			return *this;
		};
		
		///Read the molecule's type.
		int readType()
		{
			#ifdef WARNINGS_ENABLED
				if(type==-1)
					std::cout << "Warning (molecule): type isn't set!\n";
			#endif
			return type;
		};
		
		///Read the number of bonds present in molecule residue.
		int readNBond()
		{
			return bond.size();
		};
		
		///Read the number of bonds allocated in molecule residue.
		int readNBondAlloc()
		{
			return bond.capacity();
		};
		
		///Read the number of constants in molecule residue.
		int readNConstant()
		{
			return constant.size();
		};
		
		///Read the number of constants allocated in molecule residue.
		int readNConstantAlloc()
		{
			return constant.capacity();
		};
		
		///Get a pointer to the constants in the molecule residue.
		T *getConstants()
		{
			#ifdef WARNINGS_ENABLED
				if(constant.size()==0)
					std::cout << "Warning (molecule): constant is empty!\n";
			#endif
			if(constant.size()==0)
				return NULL;
			return &constant[0];
		};
		
		///Get a pointer to the bonds in the molecule residue.
		bondData *getBonds()
		{
			#ifdef WARNINGS_ENABLED
				if(bond.size()==0)
					std::cout << "Warning (molecule): bond is empty, so no bonds are alloted!\n";
			#endif
			if(bond.size()==0)
				return NULL;
			return &bond[0];
		}
		
		///Allocates a number of bonds without adding any to the structure.
		///Mostly for the purpose of preallocation when you are generating a lot of bonds.
		void allocBonds(int value)
		{
			bond.reserve(value);
		};
		
		///Allocates a number of constants without adding any to the structure.
		///Mostly for the purpose of preallocation when you are generating a lot of constants.
		void allocConstants(int value)
		{
			constant.reserve(value);
		};
		
		///Set the molecule type. See the molecule.h file for information about types.
		void setType(int value)
		{
			type=value;
		};
		
		///Set a list of bonds starting with the 0th bond.
		void setBonds(int nElements, bondData *list)
		{
			#ifdef ERRORS_ENABLED
				if(list==NULL)
				{
					std::cerr << "Error (molecule): bondData *list is null in setBonds()!\n";
					throw 0;
				}
			#endif
			for(int i=0;i<nElements;i++)
			{
				if(bond.size()>i)
					bond[i]=list[i];
				else
					bond.push_back(list[i]);
			}
		};
		
		///Set a list of constants starting with the 0th constant.
		void setConstants(int nElements, T *list)
		{
			#ifdef ERRORS_ENABLED
				if(list==NULL)
				{
					std::cerr << "Error (molecule): bondData *list is null in setConstants()!\n";
					throw 0;
				}
			#endif
			for(int i=0;i<nElements;i++)
			{
				if(constant.size()>i)
					constant[i]=list[i];
				else
					constant.push_back(list[i]);
			}
			
		};
		
		///Set an individual bond index with bond information.
		void setBond(int index, bondData value)
		{
			#ifdef ERRORS_ENABLED
				if(index>=bond.size() || index<0)
				{
					std::cout << "Error (molecule): index in setBond() is out of bounds!\n";
					std::cout << "\tindex=" << index << ", nBond=" << bond.size() << '\n';
					throw 0;
				}
			#endif
			bond[index]=value;
		};
		
		///Set an individual constant index to a value.
		void setConstant(int index, T value)
		{
			#ifdef ERRORS_ENABLED
				if(index>=constant.size() || index<0)
				{
					std::cout << "Error (molecule): index in setConstant() is out of bounds!\n";
					std::cout << "\tindex=" << index << ", nConstant=" << constant.size() << '\n';
					throw 0;
				}
			#endif
			constant[index]=value;
		};
		
		///Add a bond to the total number of bonds.
		void addBond(bondData value)
		{
			bond.push_back(value);
		};
		
		///Add a constant to the total number of constants.
		void addConstant(T value)
		{
			constant.push_back(value);
		};
		
		///Add a list of bonds to the total number of bonds.
		void addBonds(bondData *list, int value)
		{
			this->setBonds(value, list);
		};
		
		///Add a list of constants to the total number of constants.
		void addConstants(T *list, int value)
		{
			this->setConstants(value, list);
		};
		/*
		///Find the residence of the bond index of a molecule and return it. It is -1 when the index isn't present.
		///It will iterate as long as value doesn't change!
		int findBond(int value)
		{
			int index=-1;
			if(value!=currentParticle)
				currentBondIndex=0;
			
			switch(type)
			{
				case CHAIN:
					for(int currentBondIndex=currentBondIndex;currentBondIndex<mBond;currentBondIndex++)
					{
						int start=bond[currentBondIndex].s[START];
						int nChains=bond[currentBondIndex].s[NCHAINS];
						int length=bond[currentBondIndex].s[CHAINLENGTH];
						if(value>start && value<start+nChains*length)
							index=currentBondIndex;
					}
					break;
				case BOND:
					//check that it is part of a bond
					for(int currentBondIndex=currentBondIndex;currentBondIndex<mBond;currentBondIndex++)
						if(bond[currentBondIndex].s[0]==value || bond[currentBondIndex].s[1]==value)
							index=currentBondIndex;
					break;
				case BEND:
					//check that it is part of a bendable chain
					for(int currentBondIndex=currentBondIndex;currentBondIndex<mBond;currentBondIndex++)
						if(bond[currentBondIndex].s[0]==value || bond[currentBondIndex].s[1]==value || bond[currentBondIndex].s[2]==value)
							index=currentBondIndex;
					break;
				default:
					#ifdef WARNINGS_ENABLED
						std::cout << "Warning (molecule): Not a recognized bond, returning failure.\m";
					#endif
					break;
				
			};
			return index;
		};
		*/
		///Find the residence of the bond index of a molecule and return it. index.x is -1 when the index isn't present.
		///If it finds a bond, it also returns which part it belongs to (which chain or s[index.y]).
		///Don't you dare do m.findBond(particle).x, because this requires an extra search request each time you call it!
		///Instead, create a twoVector< int >, and set it to the output of findBond. Then use it.
		twoVector<int> findBond(int value)
		{
			twoVector<int> index;
			index.x=-1;
			if(value!=currentParticle)
				currentBondIndex=0;
			else
				currentBondIndex++;
			currentParticle=value;
			
			switch(type)
			{
				case CHAIN:
					for(;currentBondIndex<bond.size();currentBondIndex++)
					{
//std::cout << currentBondIndex << '\t' << bond.size() << '\t' << value << '\n';
						int start=bond[currentBondIndex].s[START];
						int nChains=bond[currentBondIndex].s[NCHAINS];
						int length=bond[currentBondIndex].s[CHAINLENGTH];
						if(value>start && value<start+nChains*length)
						{
							index.x=currentBondIndex;
							index.y=((value-start)/nChains);//which chain
							break;
						}
					}
//std::cout << "What?\n";
					break;
				case BOND:
					//check that it is part of a bond
					for(;currentBondIndex<bond.size();currentBondIndex++)
					{
//std::cout << currentBondIndex << '\n';
						if(bond[currentBondIndex].s[0]==value)
						{
							index.x=currentBondIndex;
							index.y=0;
							break;
						}
						if(bond[currentBondIndex].s[1]==value)
						{
							index.x=currentBondIndex;
							index.y=1;
							break;
						}
					}
					break;
				case BEND:
					//check that it is part of a bendable chain
					for(;currentBondIndex<bond.size();currentBondIndex++)
					{
						if(bond[currentBondIndex].s[0]==value)
						{
							index.x=currentBondIndex;
							index.y=0;
							break;
						}
						if(bond[currentBondIndex].s[1]==value)
						{
							index.x=currentBondIndex;
							index.y=1;
							break;
						}
						if(bond[currentBondIndex].s[2]==value)
						{
							index.x=currentBondIndex;
							index.y=2;
							break;
						}
					}
					break;
				case BEAD:
					for(;currentBondIndex<bond.size();currentBondIndex++)
					{
						if(bond[currentBondIndex].s[0]==value)
						{
							index.x=currentBondIndex;
							index.y=0;
						}
					}
					break;
				default:
					#ifdef WARNINGS_ENABLED
						std::cout << "Warning (molecule): Not a recognized bond, returning failure to find.\n";
					#endif
					break;
				
			};
			return index;
		};
		
		///Reset the current value for find. Every time you want to start a new search, just reset this!
		void resetFind()
		{
			currentParticle=-1;
			currentBondIndex=0;
		};
		
	private:
		int type;
		
		std::vector<T> constant;
		
		std::vector<bondData> bond;
		int currentParticle;
		int currentBondIndex;
};

/*//Can't do this?
template <typename T>
std::vector< fourVector<int> > mergeFourVectorMolecules(molecules<T, fourVector<int> > *list, int nMolecules)
{
	
}

template <typename T>
std::vector< std::vector<T> > mergeFourVectorConstants(molecules<T, fourVector<int> > *list, int nMolecules)
{
	for(int i=0;i<nMolecules;i++)
	{
		
	}
}
*/
//end of MD_BONDEDMOLECULES
#endif
