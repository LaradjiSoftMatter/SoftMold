/**
 * \brief Simple implicit molecular dynamic system. Just a bilayer in various geometries.
 * The steps to make a system are pretty simple, especially given the lipidModels.h header.
 * Just put together the parameters for initial conditions, generate the geometry, and
 * then run the system. That's it! If you want more, less, or something different, just
 * modify this program. There is a section that takes the command line arguments (just
 * after the main function definition) and checks that there are enough. After that,
 * atoi() and atof() (see some cstdlib documentation) convert the text to values. Those
 * are basically your initial conditions. The Blob class template is also useful for this.
 * Just take a parameter you want to modify and either use setParameter, addParameter, or
 * delParameter members in Blob to modify it. I've used the variable definition
 * 'Blob\<double\> System' for convienience. Functions in lipidModels.h accept System as
 * a reference, and can modify parameters as one would in main. For convienience, this
 * program outputs a name.xyz file to quickly view the initial geometry; this is useful
 * for debugging.
 */

#include <iostream>
#include <cstdlib>

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#include "../models/lipidModels.h"
#include "../models/nanoModels.h"
using namespace std;

#define CUTOFF 2.0
#define RMIN 1.0

#define HEAD 2
#define TAIL 3
#define HEAD2 6
#define TAIL2 7

#define NANOTYPE 4
//#define SINGLEBEAD 5

int main(int argc, char **argv)
{
	if(argc!=7)
	{
		std::cout << "nArgs: " << argc << std::endl;
		std::cout << "Usage: " << argv[0] << " name mIndex nanoRadius nanoHeadUmin nanoHeadUmax sigma" << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int mIndex;
	double nanoRadius=0, nanoHeadUmin=0, nanoNanoUmin=0, nanoHeadUmax=0, sigma=0;
	cmdArg >> mIndex >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> sigma;
	
	//the variables for the simulation
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	if(System.readNMolecules()<=mIndex)
	{
		std::cerr << "mIndex out of range!" << std::endl;
		std::cerr << "\tSystem: " << System.readNMolecules() << '\t' << "mIndex: " << mIndex << std::endl;
		return 0;
	}
	
	//constants for bead type
	std::vector<double> C;
	
	//defaults, 200=Umax, 0=Umin, rc=2, rmin=1, sigma=46
	for(int typeA=0;typeA<System.readNTypes();typeA++)
	{
		for(int typeB=0;typeB<System.readNTypes();typeB++)
		{
			//double rmin=1.0;
			//double rc=2.0*rmin;
			double UminNano=0.0;
			double UmaxNano=10.0;
			
			//bead-particle constants
			C.push_back(CUTOFF+nanoRadius);//0
			C.push_back(-7.0/4.0*RMIN);//1
			C.push_back(2.0*RMIN*RMIN);//2
			C.push_back(UminNano*M_PI*nanoRadius*sigma/(RMIN*RMIN*RMIN));//3
			C.push_back(nanoRadius);//4
			C.push_back(RMIN);//5
			double D=sigma*M_PI*nanoRadius;
			double A=UmaxNano-UminNano;
			
			C.push_back(-D*A/(2.0*RMIN*RMIN));//6,0@B^4
			C.push_back(2.0*D*A/(3.0*RMIN));//7,1,@B^3
			C.push_back(-D*UminNano);//8,2,@B^2
			C.push_back(2.0*D*UminNano*RMIN);//9,3,@B^1
			C.push_back(D*1.3*UminNano*RMIN*RMIN);//10,4
			
			//bead-bead constants
			D=pow(2.0*M_PI*sigma*nanoRadius,2.0)/(RMIN*RMIN*RMIN);
			
			C.push_back(2.0*nanoRadius+RMIN);//0
			C.push_back(2.0*nanoRadius+CUTOFF);//1
			C.push_back(UminNano*D*2.0/30.0);//2, x^6
			C.push_back(-UminNano*D*7.0*RMIN/20.0);//3, x^5
			C.push_back(UminNano*D*RMIN*RMIN/2.0);//4, x^4
			
			C.push_back(-D*A*RMIN/20.0);//5,0@B^5
			C.push_back(D*A*RMIN*RMIN/12.0);//6,1,@B^4
			C.push_back(-D*UminNano*pow(RMIN,3.0)/6.0);//7,2,@B^3
			C.push_back(D*UminNano*pow(RMIN,4.0)/2.0);//8,3,@B^2
			C.push_back(D*1.3*UminNano*pow(RMIN,5.0)/2.0);//9,@B
			C.push_back(D*13.0*UminNano*pow(RMIN,6.0)/60.0);//10
		}
	}
	
	molecule< double, fourVector<int> > *m=System.getMolecule();
	double *mConstants=m[mIndex].getConstants();
	
	//exceptions
	//nanoparticle to head
	if(m[mIndex].readNBond()>0)
	{
		int cindex=nBEADCONST*(HEAD+NANOTYPE*System.readNTypes());
		int cindexMirror=nBEADCONST*(NANOTYPE+HEAD*System.readNTypes());
		
		double UminNano=nanoHeadUmin;
		double UmaxNano=nanoHeadUmax;
		
		//bead-particle constants
		C[cindex+0]=(CUTOFF+nanoRadius);//c0
		C[cindex+1]=(-7.0/4.0*RMIN);//c1
		C[cindex+2]=(2.0*RMIN*RMIN);//c2
		C[cindex+3]=(UminNano*M_PI*nanoRadius*sigma/(RMIN*RMIN*RMIN));//c3
		C[cindex+4]=(nanoRadius);//c4
		C[cindex+5]=(RMIN);//c5
		
		double D=sigma*M_PI*nanoRadius;
		double A=UmaxNano-UminNano;
			
		C[cindex+6]=-D*A/(2.0*RMIN*RMIN);//c6
		C[cindex+7]=2.0*D*A/(3.0*RMIN);//c7
		C[cindex+8]=-D*UminNano;
		C[cindex+9]=2.0*D*UminNano*RMIN;
		C[cindex+10]=D*1.3*UminNano*RMIN*RMIN;
		
		//bead-bead constants
		D=pow(2.0*M_PI*sigma*nanoRadius,2.0)/(RMIN*RMIN*RMIN);
		
		C[cindex+11]=(2.0*nanoRadius+RMIN);//0
		C[cindex+12]=(2.0*nanoRadius+CUTOFF);//1
		C[cindex+13]=(UminNano*D*2.0/30.0);//2, x^6
		C[cindex+14]=(-UminNano*D*7.0*RMIN/20.0);//3, x^5
		C[cindex+15]=(UminNano*D*RMIN*RMIN/2.0);//4, x^4
		
		C[cindex+16]=(-D*A*RMIN/20.0);//5,0@B^5
		C[cindex+17]=(D*A*RMIN*RMIN/12.0);//6,1,@B^4
		C[cindex+18]=(-D*UminNano*pow(RMIN,3.0)/6.0);//7,2,@B^3
		C[cindex+19]=(D*UminNano*pow(RMIN,4.0)/2.0);//8,3,@B^2
		C[cindex+20]=(D*1.3*UminNano*pow(RMIN,5.0)/2.0);//9,@B
		C[cindex+21]=(D*13.0*UminNano*pow(RMIN,6.0)/60.0);//10
		
		for(int i=0;i<nBEADCONST;i++)
			C[cindexMirror+i]=C[cindex+i];
	}
	
	
	for(int i=0;i<C.size();i++)
		mConstants[i]=C[i];
	
	std::cerr << "Storing configuration...";
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


