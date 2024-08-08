#include "include.h"

#define SURFA 5
#define SURFB 6


#define NANOTYPE 4
//#define SINGLEBEAD 5

int main(int argc, char **argv)
{
	if(argc!=14 && argc!=10 && argc!=6 && argc!=11)
	{
		std::cout << "nArgs: " << argc << std::endl;
		//14, bilayer with cytoskeleton and nanoparticle
		std::cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		std::cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles";
		std::cout << "nMonomers nAnchors.x nAnchors.y anchorHeadUmin " << std::endl;
		std::cout << "\tFor bilayer with cytoskeleton and nanoparticles near surface." << std::endl << std::endl;
		
		//6, nanoparticles only
		std::cout << "Usage: " << argv[0] << " name seed nNanoparticles nanoRadius nanoDensity\n";
		std::cout << "\tFor nanoparticles only." << std::endl << std::endl;
		
		//10, bilayer, nanoparticle
		std::cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		std::cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles" << std::endl;
		std::cout << "\tFor bilayer and nanoparticles." << std::endl << std::endl;
		
		//11, bilayer, nanoparticle, specific spacing
		std::cout << "Usage: " << argv[0] << " name seed nLipids lipidArealDensity ";
		std::cout << "nanoRadius nanoHeadUmin nanoHeadUmax sigma nNanoparticles dNanoparticles" << std::endl;
		std::cout << "\tFor bilayer and nanoparticles with a specific spacing (chain configuration)." << std::endl << std::endl;
		
		return 0;
	}
	
	char *name=argv[1];
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	int seed, nLipids=0, nMonomers=0, nNanoparticles;
	double lipidArealDensity=0, anchorHeadUmin=0, nanoRadius=0, nanoDensity=0, dNanoparticles=0, nanoHeadUmin=0;
	double nanoNanoUmin=0, nanoHeadUmax=0, sigma=0;
	twoVector<int> nAnchors;
	
	if(argc>=3)
		cmdArg >> seed;
	if(argc==6)
		cmdArg >> nNanoparticles >> nanoRadius >> nanoDensity;
	else if(argc>=8)
	{
		cmdArg >> nLipids >> lipidArealDensity >> nanoRadius >> nanoHeadUmin >> nanoHeadUmax >> sigma >> nNanoparticles;
		if(argc==11)
			cmdArg >> dNanoparticles;
		if(argc==14)
			cmdArg >> nMonomers >> nAnchors.x >> nAnchors.y >> anchorHeadUmin;
	}
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(7);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(1000);
	System.setDeltaT(0.008);
	System.setStoreInterval(100);
	System.setMeasureInterval(100);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(3.0);
	System.setFinalTemp(3.0);
	//System.setSolventGamma(5.0);
	
	//initialize constants (force and potential constants)
	//generate force constants
        std::vector<double> Umin(System.readNTypes()*System.readNTypes(),0.0);
        std::vector<double> Umax(System.readNTypes()*System.readNTypes(),0.0);
	
	//could be one loop, but it's good for an example
	//of how to set constants
	for(int i=0;i<System.readNTypes();i++)
	{
		for(int j=0;j<System.readNTypes();j++)
		{
			//i=first type, j=second type, indexed grid
			int k=i+j*System.readNTypes();
			Umin[k]=0;
			Umax[k]=100;
		}
	}
	
	//tail types hydrophobicity
	Umin[TAIL+TAIL*System.readNTypes()]=-6;
	Umax[TAIL+TAIL*System.readNTypes()]=200;
        
	Umin[SURFA+SURFA*System.readNTypes()]=0;
	Umax[SURFA+SURFA*System.readNTypes()]=20;
	
	Umin[SURFA+SURFB*System.readNTypes()]=0;
	Umax[SURFA+SURFB*System.readNTypes()]=20;
	
	Umin[SURFB+SURFA*System.readNTypes()]=Umin[SURFA+SURFB*System.readNTypes()];
	Umax[SURFB+SURFA*System.readNTypes()]=Umax[SURFA+SURFB*System.readNTypes()];
	
	Umin[SURFA+HEAD*System.readNTypes()]=0;
	Umax[SURFA+HEAD*System.readNTypes()]=40;
	
	Umin[HEAD+SURFA*System.readNTypes()]=Umin[SURFA+HEAD*System.readNTypes()];
	Umax[HEAD+SURFA*System.readNTypes()]=Umax[SURFA+HEAD*System.readNTypes()];
	
	Umin[SURFB+HEAD*System.readNTypes()]=-2;
	Umax[SURFB+HEAD*System.readNTypes()]=40;
	
	Umin[HEAD+SURFB*System.readNTypes()]=Umin[SURFB+HEAD*System.readNTypes()];
	Umax[HEAD+SURFB*System.readNTypes()]=Umax[SURFB+HEAD*System.readNTypes()];
	
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//double radius=sqrt(((double)nLipids/lipidArealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;//80;
	size.z=80;
	
	if(nanoDensity>0)
	{
		double volume=nNanoparticles/nanoDensity;
		size.x=pow(volume,1.0/3.0);
		size.y=pow(volume,1.0/3.0);
		size.z=pow(volume,1.0/3.0);
		std::cerr << "Setting system size: " << size.x << '\t' << size.y << '\t' << size.z << std::endl;
	}
	
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;//System.readSize().x/2.0;
	pos.y=0;//System.readSize().y/2.0;
	pos.z=System.readSize().z/2.0;
	double bondLength=0.7;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=100;
	
	//bilayer with cytoskeleton, regular hexagonal cytoskeleton has an aspect ratio sqrt(3/2)=y/x
	if(nLipids>0)
	{
		if(nAnchors.x<=0 || nAnchors.y<=0 || nMonomers<=0)
			bilayer<double>(System, nLipids, 3, pos, bondLength, lipidArealDensity, constants, 4, 1.0);
		else
			bilayer<double>(System, nLipids, 3, pos, bondLength, lipidArealDensity, constants, 4, (double(nAnchors.x)/double(nAnchors.y))/HEXAGONAL_ASPECT_RATIO);
	}
	//pos.z-=bondLength;
	if(nAnchors.x>0 && nAnchors.y>0 && nMonomers>0)
		flatHexagonalCyto<double>(System, nAnchors, nMonomers, pos, constants, 4);
	
	//Add nanoparticle constants
	std::vector<double> UminContinuum(System.readNTypes()*System.readNTypes(),0.0);
	std::vector<double> UmaxContinuum(System.readNTypes()*System.readNTypes(),10.0);
	
	//exceptions
	UminContinuum[SURFB+NANOTYPE*System.readNTypes()]=nanoHeadUmin;
	UminContinuum[NANOTYPE+SURFB*System.readNTypes()]=nanoHeadUmin;
	UmaxContinuum[SURFB+NANOTYPE*System.readNTypes()]=nanoHeadUmax;
	UmaxContinuum[NANOTYPE+SURFB*System.readNTypes()]=nanoHeadUmax;
	UminContinuum[SURFA+NANOTYPE*System.readNTypes()]=nanoHeadUmin;
	UminContinuum[NANOTYPE+SURFA*System.readNTypes()]=nanoHeadUmin;
	UmaxContinuum[SURFA+NANOTYPE*System.readNTypes()]=nanoHeadUmax;
	UmaxContinuum[NANOTYPE+SURFA*System.readNTypes()]=nanoHeadUmax;
	
	std::vector<double> C=mpd::laradjiSpanglerFC(UmaxContinuum, UminContinuum, CUTOFF, 
						     RMIN, nanoRadius-1, 5.88);
	//Add a nanoparticle
	std::vector< threeVector<double> > nanoPos;
	MTRand randNum(System.readSeed());
	
	for(int i=0;i<nNanoparticles;i++)
	{
		std::cout << "Placing nanoparticle " << i << "!" << std::endl;
		threeVector<double> toPlace;
		bool overlap;
		do
		{
			overlap=false;
			if(dNanoparticles<=0)//no chain
			{
				toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				if(nanoDensity<=0)//bilayer
					toPlace.z=pos.z+nanoRadius+2.0+3.0*2.0*0.7;
				else//no bilayer
					toPlace.z=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
			}
			else//chain placement
			{
				double xPos, yPos;
				bool outOfBounds;
				if(i!=0)
				{
					do
					{
						outOfBounds=false;
						double theta=2.0*M_PI*randNum.rand53();
						xPos=sin(theta)*(dNanoparticles+2.0*nanoRadius);
						yPos=cos(theta)*(dNanoparticles+2.0*nanoRadius);
						if(nanoPos[i-1].x+xPos>System.readSize().x-nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].y+yPos>System.readSize().y-nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].x+xPos<nanoRadius)
							outOfBounds=true;
						if(nanoPos[i-1].y+yPos<nanoRadius)
							outOfBounds=true;
					} while(outOfBounds);
					
					toPlace.x=xPos+nanoPos[i-1].x;
					toPlace.y=yPos+nanoPos[i-1].y;
					toPlace.z=pos.z+nanoRadius+2.0+3.0*2.0*0.7;
				}
				else
				{
					toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
					toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
					toPlace.z=pos.z+nanoRadius+2.0+3.0*2.0*0.7;
				}
			}
			for(int j=0;j<nanoPos.size();j++)
			{
				threeVector<double> d;
				d.x=nanoPos[j].x-toPlace.x;
				d.y=nanoPos[j].y-toPlace.y;
				d.z=nanoPos[j].z-toPlace.z;
				d.x-=(d.x>System.readSize().x/2.0)?System.readSize().x:0;
				d.x+=(d.x<-System.readSize().x/2.0)?System.readSize().x:0;
				
				d.y-=(d.y>System.readSize().y/2.0)?System.readSize().y:0;
				d.y+=(d.y<-System.readSize().y/2.0)?System.readSize().y:0;
				
				d.z-=(d.z>System.readSize().z/2.0)?System.readSize().z:0;
				d.z+=(d.z<-System.readSize().z/2.0)?System.readSize().z:0;
				
				if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<2.0*nanoRadius+2.0)
					overlap=true;
			}
		} while(overlap);
		std::cerr << i << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
		nanoPos.push_back(toPlace);
	}
	continuumSphere<double>(System, &(nanoPos[0]), nanoPos.size(), nanoRadius, &C[0], 
				System.readNTypes()*System.readNTypes()*nBEADCONST, NANOTYPE);
	
	for(int i=0;i<nanoPos.size();i++)
	{
		int startIndex=System.readNParticles();
		shellNano(System, nanoPos[0], nanoRadius+1.0, 2, SURFA);
		int endIndex=System.readNParticles();
		for(int j=startIndex+1;j<endIndex;j++)
			System.getPositions()[j].type=(j%2==0)?SURFA:SURFB;
	}
	
	std::cerr << "Storing configuration...";
	
	//Write configuration
	Script<double,Blob <double> > output(name,std::ios::out,&System);
	output.write();
	
	//This stores the xyz file to check new configuration
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	std::string newName("");
	newName+=name;
	newName+=".xyz";
	
	xyzFormat<double> xyzFile(p, nParticles);
	xyzFile.open(newName.c_str(), std::ios::out);
	xyzFile.store();
	xyzFile.close();
	
	std::cerr << "Done.\nExiting...\n";
	return 0;
}


