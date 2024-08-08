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

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		std::cout << argv[0] << " name typeA typeB Umin Umax\n";
		return 0;
	}

	int typeA, typeB;
	double Umin, Umax;
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	cmdArg >> typeA >> typeB >> Umin >> Umax;
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	double * uConst=System.getTwoBodyUconst();
	double * fConst=System.getTwoBodyFconst();
	int types=System.readNTypes();
	int offsetA=typeA+typeB*types;
	int offsetB=typeB+typeA*types;
	double rmin=1.0, cutoff=2.0;
	
	fConst[offsetA*nTWOBODYFCONST+0]=rmin;
	fConst[offsetA*nTWOBODYFCONST+1]=((2.0*(Umax-Umin))/(rmin*rmin));
	fConst[offsetA*nTWOBODYFCONST+2]=0;
	fConst[offsetA*nTWOBODYFCONST+3]=cutoff;
	fConst[offsetA*nTWOBODYFCONST+4]=((6.0*Umin)/((cutoff-rmin)*(cutoff-rmin)));
	fConst[offsetA*nTWOBODYFCONST+5]=((6.0*Umin)/((cutoff-rmin)*(cutoff-rmin)*(cutoff-rmin)));
	for(int i=0;i<nTWOBODYFCONST;i++)
		fConst[offsetB*nTWOBODYFCONST+i]=fConst[offsetA*nTWOBODYFCONST+i];
	
	uConst[offsetA*nTWOBODYUCONST+0]=rmin;
	uConst[offsetA*nTWOBODYUCONST+1]=(((Umax-Umin))/(rmin*rmin));
	uConst[offsetA*nTWOBODYUCONST+2]=Umin;
	uConst[offsetA*nTWOBODYUCONST+3]=cutoff;
	uConst[offsetA*nTWOBODYUCONST+4]=((3.0*Umin)/((cutoff-rmin)*(cutoff-rmin)));
	uConst[offsetA*nTWOBODYUCONST+5]=((2.0*Umin)/((cutoff-rmin)*(cutoff-rmin)*(cutoff-rmin)));
	
	for(int i=0;i<nTWOBODYUCONST;i++)
		uConst[offsetB*nTWOBODYUCONST+i]=uConst[offsetA*nTWOBODYUCONST+i];
	
	
	//Save it
	fileIO.open(argv[1],std::ios::out);
	fileIO.write();
	fileIO.close();
	

	return 0;
}
