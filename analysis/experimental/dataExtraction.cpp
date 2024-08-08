//#define ANCHOR_DATA

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

//For extraction of molecular dynamics data (measurements)
#include "../dataExtractionLimited.h"

#define CUTOFF 2.0

int main(int argc, char* argv[])
{
	if(argc!=2 && argc!=3)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "usage: " << argv[0] << " name\n";
		std::cerr << "Extracts data in place of MD.cpp for already executed systems.\n";
		
		//specify a constant
		std::cerr << "usage: " << argv[0] << " name profileCutoff\n";
		std::cerr << "Extracts data in place of MD.cpp for already executed systems.\n";
		
		std::cout << "Usage: " << argv[0] << " name mIndex nanoRadius nanoHeadUmin nanoHeadUmax sigma" << std::endl;
		
		return 0;
	}
	
	///read in options
	char *name=argv[1];
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	double profileCutoff=-1;
	
	if(argc==3)
		cmdArg >> profileCutoff;
	
	int mIndex;
	double nanoRadius=0, nanoHeadUmin=0, nanoNanoUmin=0, nanoHeadUmax=0, sigma=0;
	
	if(argc==7)
		cmdArg >> mIndex >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> sigma;
	
	///Configuration variables
	Blob<double> System;
	
	///load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	
	fileIO.read();
	
	fileIO.close();
	
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	System.setInitialTime(0.0);
	
	//for data
	dataExtraction<double, Blob <double> > dataCollection(&System,name,profileCutoff);
	dataCollection.initialize();
	
	std::string sizeName("size_");
	sizeName+=argv[1];
	sizeName+=".dat";
		
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent=System.readSize();
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		System.setSize(sCurrent);
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	std::string beadName("beadUmin_");
	beadName+=argv[1];
	beadName+=".dat";
		
	std::fstream beadFile;
	beadFile.open(beadName.c_str(), std::ios::in);
	if(beadFile.is_open() && argc==7)
	{
		double sTime=-1;
		threeVector<double> sCurrent=System.readSize();
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !beadFile.eof())
			beadFile >> sTime >> beadUmin;
		
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	if(System.readNMolecules()<=mIndex && argc==7)
	{
		std::cerr << "mIndex out of range!" << std::endl;
		std::cerr << "\tSystem: " << System.readNMolecules() << '\t' << "mIndex: " << mIndex << std::endl;
		return 0;
	}
	
	if(argc==7)
	{
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
	}
	
	double simTime=0;
	for(int i=0;xyzFile.load();i++)
	{
		std::cerr << static_cast<double>(i)*System.readStoreInterval() << std::endl;
		System.setInitialTime(static_cast<double>(i)*System.readStoreInterval());
		
		simTime=System.readInitialTime();
		
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=System.readSize();//just in case it is current
			while(sTime<simTime && !sizeFile.eof())
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			System.setSize(sCurrent);
		}
		
		//make sure the beadPotential hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=System.readSize();//just in case it is current
			while(sTime<simTime && !sizeFile.eof())
				beadUmin >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
			System.setSize(sCurrent);
		}
		
		dataCollection.compute();
	}
	
	return 0;
}
