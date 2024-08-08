//This contains a molecule object that can be used as a container for molecules.
//It has numerous useful functions:
//	accepts simple datatypes for bond structures
//	add bonds to heap
//	set bonds in heap
//	add constants to heap
//	set constants in heap
//It needs some useful functions:
//	add bonds unless present
//	delete bonds
//	delete constants
//	sort bonds for rapid location
//It has some issues:
//	when sending it too a function:
//		molecule mol;
//		someFunc(mol);
//		definition should be: returnType someFunc(molecule &mol)
//		definition will break if: returnType someFunc(molecule mol)
//		because allocation/deallocation occurs within the molecule class

#include "functions.h"

//molecules, soft version of a rigid construct
template <typename T, typename bondData>
class molecule {
	public:
		//Constructor: initializes values
		molecule()
		{
			bond=NULL;
			constant=NULL;
			nBond=0;
			nBondAlloc=0;
			nConstant=0;
			nConstantAlloc=0;
			type=-1;
		};
		
		//Destructor: deallocates memory if allocated
		~molecule()
		{
			if(nBondAlloc)
				delete bond;
			if(nConstantAlloc)
				delete constant;
		};
		
		//Sets molecule to another molecule's values (topical, doesn't transfer amount of memory allocated)
		molecule operator = (molecule &mol)
		{
			setType(mol.getType());
			setBond(mol.getBonds(), mol.getNBond());
			setConstant(mol.getConstants(), mol.getNConstant());
		};
		
		//All of these getSomething functions return the parameter requested
		int getType()
		{
			#ifdef WARNINGS_ENABLED
				if(type==-1)
					std::cout << "Warning (molecule): type isn't set!\n";
			#endif
			return type;
		};
		int getNBond()
		{
			#ifdef WARNINGS_ENABLED
				if(nBond==0)
					std::cout << "Warning (molecule): nBond is 0!\n";
			#endif
			return nBonded;
		};
		int getNConstant()
		{
			#ifdef WARNINGS_ENABLED
				if(nConstant==0)
					std::cout << "Warning (molecule): nConstant is 0!\n";
			#endif
			return nConstant;
		};
		T *getConstants()
		{
			#ifdef WARNINGS_ENABLED
				if(constant==NULL)
					std::cout << "Warning (molecule): constant is null!\n";
			#endif
			return constant;
		};
		bondData *getBonds()
		{
			#ifdef WARNINGS_ENABLED
				if(bond==NULL)
					std::cout << "Warning (molecule): bond is null!\n";
			#endif
			return bond;
		}
		
		//These allocate if you know how much to allocate in advance
		void allocBonds(int value)
		{
			if(nBondAlloc<value)
			{
				if(nBondAlloc)
				{
					bondData *buf=new bondData[value];
					for(int i=0;i<nBondAlloc;i++)
						buf[i]=bond[i];
					delete bond;
					bond=buf;
					nBondAlloc=value;
				}
				else
				{
					bond=new bondData[value];
					nBondAlloc=value;
				}
			}
			else
			{
				#ifdef WARNINGS_ENABLED
					std::cout << "Warning (molecule): value is less than number of allocated bonds!";
					std::cout << "\tNothing will be done!\n";
				#endif
			}
		};
		void allocConstant(int value)
		{
			if(nConstantAlloc<value)
			{
				if(nConstantAlloc)
				{
					T *buf=new T[value];
					for(int i=0;i<nConstantAlloc;i++)
						buf[i]=constant[i];
					delete constant;
					constant=buf;
					nConstantAlloc=value;
				}
				else
				{
					constant=new T[value];
					nConstantAlloc=value;
				}
			}
			else
			{
				#ifdef WARNINGS_ENABLED
					std::cout << "Warning (molecule): value is less than number of allocated constants!";
					std::cout << "\tNothing will be done!\n";
				#endif
			}
		};
		
		//These will set a particular parameter
		void setType(int value)
		{
			type=value;
		};
		void setBond(bondData *list, int value)
		{
			if(nBondAlloc<value)
			{
				if(nBondAlloc)
					delete bond;
				bond=new bondData[value];
				nBondAlloc=value;
			}
			for(int i=0;i<value;i++)
				bond[i]=list[i];
			nBond=value;
		};
		void setConstant(T *list, int value)
		{
			if(nConstantAlloc<value)
			{
				if(nConstantAlloc)
					delete constant;
				constant=new T[value];
				nConstantAlloc=value;
			}
			for(int i=0;i<value;i++)
				constant[i]=list[i];
			nConstant=value;
			
		};
		void setBond(bondData value, int offset)
		{
			if(offset>=nBond || offset<0)
			{
				#ifdef ERRORS_ENABLED
					std::cout << "Error (molecule): offset in setBond is out of bounds!\n";
					std::cout << "\toffset=" << offset << ", nBond=" << nBond << '\n';
				#endif
				throw 0;
			}
			bond[offset]=value;
		};
		void setConstant(T value, int offset)
		{
			if(offset>=nConstant || offset<0)
			{
				#ifdef ERRORS_ENABLED
					std::cout << "Error (molecule): offset in setConstant is out of bounds!\n";
					std::cout << "\toffset=" << offset << ", nConstant=" << nConstant << '\n';
				#endif
				throw 0;
			}
			constant[offset]=value;
		};
		
		//These will add parameters to heap
		void addBond(bondData value)
		{
			if(nBond+1<nBondAlloc)
			{
				bondData *buf=new bondData[nBond+1];
				for(int i=0;i<nBond;i++)
					buf[i]=bond[i];
				delete bond;
				bond=buf;
			}
			bond[nBonded++]=value;
		};
		void addConstant(T value)
		{
			if(nConstant+1<nConstantAlloc)
			{
				T *buf=new T[nConstant+1];
				for(int i=0;i<nConstant;i++)
					buf[i]=constant[i];
				delete constant;
				constant=buf;
			}
			constant[nConstant++]=value;
		};
		void addBonds(bondData *list, int value)
		{
			if(nBond+value<nBondAlloc)
			{
				bondData *buf=new bondData[nBond+value];
				for(int i=0;i<nBond;i++)
					buf[i]=bond[i];
				delete bond;
				bond=buf;
			}
			for(int i=0;i<value;i++)
				bond[nBond++]=list[i];
		};
		void addConstants(T *list, int value)
		{
			if(nConstant+value<nConstantAlloc)
			{
				T *buf=new T[nConstant+value];
				for(int i=0;i<nConstant;i++)
					buf[i]=constant[i];
				delete constant;
				constant=buf;
			}
			for(int i=0;i<value;i++)
				constant[nConstant++]=list[i];
		};
		
	private:
		int type;
		
		T *constant;
		int nConstant;
		int nConstantAlloc;
		
		bondData *bond;
		int nBond;
		int nBondAlloc;
};

#ifndef MD_BONDEDMOLECULES
#define MD_BONDEDMOLECULES


//end of MD_BONDEDMOLECULES
#endif
