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
	if(argc!=9)
	{
		std::cout << "nArgs: " << argc << std::endl;
		std::cout << "Usage: " << argv[0] << " name mIndex nanoRadius nanoHeadUmin nanoHeadUmax cutoff sigma type" << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int mIndex,type;
	double nanoRadius=0, nanoHeadUmin=0, nanoNanoUmin=0, nanoHeadUmax=0, sigma=0,cutoff=0;
	cmdArg >> mIndex >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> cutoff >> sigma >> type;
	double rmin=cutoff/2.0;
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
	
	molecule< double, fourVector<int> > *m=System.getMolecule();
	double *mConstants=m[mIndex].getConstants();
	int nanoType=System.getPositions()[m[mIndex].getBonds()[0].s[0]].type;
	
	//exceptions
	//nanoparticle to head
	if(m[mIndex].readNBond()>0)
	{
		int cindex=nBEADCONST*(type+nanoType*System.readNTypes());
		int cindexMirror=nBEADCONST*(nanoType+type*System.readNTypes());
		
		double UminNano=nanoHeadUmin;
		double UmaxNano=nanoHeadUmax;
		
		//bead-particle constants
		mConstants[cindex+0]=(cutoff+nanoRadius);//c0
		mConstants[cindex+1]=(-7.0/4.0*rmin);//c1
		mConstants[cindex+2]=(2.0*rmin*rmin);//c2
		mConstants[cindex+3]=(UminNano*M_PI*nanoRadius*sigma/(rmin*rmin*rmin));//c3
		mConstants[cindex+4]=(nanoRadius);//c4
		mConstants[cindex+5]=(rmin);//c5
		
		double D=sigma*M_PI*nanoRadius;
		double A=UmaxNano-UminNano;
			
		mConstants[cindex+6]=-D*A/(2.0*rmin*rmin);//c6
		mConstants[cindex+7]=2.0*D*A/(3.0*rmin);//c7
		mConstants[cindex+8]=-D*UminNano;
		mConstants[cindex+9]=2.0*D*UminNano*rmin;
		mConstants[cindex+10]=D*1.3*UminNano*rmin*rmin;
		
		//bead-bead constants
		D=pow(2.0*M_PI*sigma*nanoRadius,2.0)/(rmin*rmin*rmin);
		
		mConstants[cindex+11]=(2.0*nanoRadius+rmin);//0
		mConstants[cindex+12]=(2.0*nanoRadius+cutoff);//1
		mConstants[cindex+13]=(UminNano*D*2.0/30.0);//2, x^6
		mConstants[cindex+14]=(-UminNano*D*7.0*rmin/20.0);//3, x^5
		mConstants[cindex+15]=(UminNano*D*rmin*rmin/2.0);//4, x^4
		
		mConstants[cindex+16]=(-D*A*rmin/20.0);//5,0@B^5
		mConstants[cindex+17]=(D*A*rmin*rmin/12.0);//6,1,@B^4
		mConstants[cindex+18]=(-D*UminNano*pow(rmin,3.0)/6.0);//7,2,@B^3
		mConstants[cindex+19]=(D*UminNano*pow(rmin,4.0)/2.0);//8,3,@B^2
		mConstants[cindex+20]=(D*1.3*UminNano*pow(rmin,5.0)/2.0);//9,@B
		mConstants[cindex+21]=(D*13.0*UminNano*pow(rmin,6.0)/60.0);//10
		
		for(int i=0;i<nBEADCONST;i++)
			mConstants[cindexMirror+i]=mConstants[cindex+i];
	}
	
	
	std::cerr << "Storing configuration...";
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


