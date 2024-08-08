#include "functions.h"

#ifndef MD_BONDEDMOLECULES
#define MD_BONDEDMOLECULES

/** \brief Simple molecule residue container.
 * Simple molecule residue container that holds the molecule type, bonds, and constants.
 * Usually you would do something like this:
 * 	enum molID {BOND,BEND,CHAIN};
 * 	vector\< residue\<molID,fourVector\<int\>,double\> \> residueList;
 * 	residue\<molID,fourVector\<int\>,double\> newResidue;
 * 	newResidue.type=BOND;
 * 	fourVector\<int\> bond;
 * 	bond.s[0]=0;
 * 	bond.s[1]=1;
 * 	newResidue.bonds.push_back(bond);
 * 	newResidue.bonds.push_back(100);
 * 	newResidue.bonds.push_back(0.7);
 * 	residueList.push_back(newResidue);
 * 
 * You should place molID in the system header, the header that handles your system information.
 */
template <typename ID, typename BONDDATA, typename BONDCONSTANTS>
struct residue {
	///Type of residue to identify it in the handler.
	ID type;
	///Bond information.
	std::vector<BONDDATA> bonds;
	///Constants associated with the residue.
	std::vector<BONDCONSTANTS> constants;
};

/*
//maybe I should do something like this:
template <typename ID, typename BONDDATA, typename BONDCONSTANTS>
class simpleMolecule {
	public:
		///construct molecule
		simpleMolecule() {
			myType=-1;
		};
		
		///construct molecule with another molecule
		simpleMolecule(const simpleMolecule &m) {
			this->type(m.type());
			this->bonds(m.bonds());
			this->constants(m.constants());
		};
		
		///copy another molecule to the current molecule, overwrites it
		simpleMolecule& operator= (const simpleMolecule &m) {
			this->type(m.type());
			this->bonds(m.bonds());
			this->constants(m.constants());
		};
		
		///Destructor.
		~simpleMolecule();
		
		
		
		///return mytype
		ID &type() {
			return myType;
		};
		
		///return bond information from index i
		BONDDATA &bond(const int &i) {
			return myBonds[i];
		};
		
		
		///return constant information from index i
		BONDCONSTANTS &constant(const int &i) {
			return myConstants[i];
		};
		
		
		
		///Examine, modify, and pass bond information.
		std::vector<BONDDATA> &bonds() {
			return myBonds;
		};
		
		///Examine, modify, and pass constant information
		std::vector<BONDCONSTANTS> &constants() {
			return myConstants;
		};
		
		///Set type
		void type(const ID &newType) {
			myType=newType;
		};
		
		///Set bonds. We don't mind if you modify or reference our bonds, 
		/// just don't make it impossible for us to reference it.
		void bonds(const std::vector<BONDDATA> &newBonds) {
			myBonds=newBonds;
		};
		
		void constants(const std::vector<BONDCONSTANTS> &newConstants) {
			myConstants=newConstants;
		};
		
		
		
		
		std::vector<BONDDATA> addBond(const BONDDATA &newBond) {
			myBonds.push_back(newBond);
			return myBonds;
		};
		
		std::vector<BONDCONSTANTS> addConstant(const BONDCONSTANTS &newConstant) {
			myConstants.push_back(newConstant);
			return myConstants;
		};
		
		std::vector<BONDDATA> delBond(const int &i) {
			myBonds.erase(myBonds.begin()+i);
			return myBonds;
		};
		
		std::vector<BONDCONSTANTS> delConstant(const int &i) {
			myConstants.erase(myConstants.begin()+i);
			return myConstants;
		};
		
		
		
		std::vector<BONDDATA> addBonds(const std::vector<BONDDATA> &newBonds) {
			myBonds.push_back(newBond);
			return myBonds;
		};
		
		std::vector<BONDCONSTANTS> addConstants(const std::vector<BONDCONSTANTS> &newConstants) {
			myConstants.push_back(newConstant);
			return myConstants;
		};
		
		std::vector<BONDDATA> delBonds(const int &start, const int &end) {
			myBonds.erase(myBonds.begin()+start,myBonds.begin()+end);
			return myBonds;
		};
		
		std::vector<BONDCONSTANTS> delConstants(const int &start, const int &end) {
			myConstants.erase(myConstants.begin()+start, myConstants.begin()+end);
			return myConstants;
		};
		
	private:
		ID myType;
		std::vector<BONDDATA> myBonds;
		std::vector<BONDCONSTANTS> myConstants;
};
*/

//end of MD_BONDEDMOLECULES
#endif
