#include "include.h"

#define BB 0
#define BM 1
#define NANOTYPE 2

int main(int argc, char **argv)
{
	if(argc!=14)
	{
		//14, brush + NP
		std::cout << "Usage: " << argv[0] << " name seed nBrushes arealDensity UminBrushNP UmaxBrushNP nNanoparticles nanoRadius UminBrushBrush UmaxBrushBrush nMonomers bondLength temp" << std::endl;
		return 0;
	}
	char *name=argv[1];
	std::stringstream cmdArg;
	
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << ' ';
	
	int seed,nBrushes,nMonomers,nNanoparticles;
	double arealDensity, UminBrushNP, UmaxBrushNP, UminBrushBrush,UmaxBrushBrush,bondLength,temp,nanoRadius;
	
	cmdArg >> seed >> nBrushes >> arealDensity >> UminBrushNP >> UmaxBrushNP >> nNanoparticles >> nanoRadius >> UminBrushBrush >> UmaxBrushBrush >> nMonomers >> bondLength >> temp;
	
	//the variables for the simulation
	Blob<double> System;
	
	System.setGamma(1.0);
	System.setNTypes(4);
	System.setSeed(seed);
	
	threeVector<bool> wrap;//periodic boundaries flag, be forwarned, this doesn't work yet, it is always on
	wrap.x=true;
	wrap.y=true;
	wrap.z=true;
	System.setPeriodic(wrap);
	System.setCutoff(2.0);
	
	System.setInitialTime(0);
	System.setFinalTime(100);
	System.setDeltaT(0.02);
	System.setStoreInterval(10);
	System.setMeasureInterval(1);
	//System.setDeltaLXY(0.01);//if you don't set it, it doesn't activate
	System.setInitialTemp(temp);
	System.setFinalTemp(temp);
	
	//initialize constants (force and potential constants)
	//generate force constants
	std::vector<double> Umax(System.readNTypes()*System.readNTypes());
	std::vector<double> Umin(System.readNTypes()*System.readNTypes());
	
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
	
	//Brush monomers
	Umin[BM+BM*System.readNTypes()]=UminBrushBrush;
	Umax[BM+BM*System.readNTypes()]=UmaxBrushBrush;
	
	//two body constants
	std::vector<double> twoBodyFconst=mpd::laradjiRevaleeFC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyFconst)
		System.addTwoBodyFconst(v);
	
	std::vector<double> twoBodyUconst=mpd::laradjiRevaleePC(Umax,Umin,CUTOFF,RMIN);
	for(auto &v:twoBodyUconst)
		System.addTwoBodyUconst(v);
	
	//double radius=sqrt(((double)nBrushes/arealDensity)/(4.0*M_PI));//density adjusted
	threeVector<double>size;
	size=0;
	size.z=static_cast<double>(nMonomers)*bondLength*2.0;
	System.setSize(size);
	
	threeVector<double> pos;
	pos.x=0;
	pos.y=0;
	pos.z=1.0;
	
	double *constants=new double[4];
	constants[0]=bondLength;
	constants[1]=100;
	constants[2]=1;
	constants[3]=0;
	
	//simple brush
	brush<double>(System, nBrushes, nMonomers, pos, bondLength, arealDensity, constants, 4, 1.0, BB, BM);
	
	//Add a nanoparticle
	MTRand randNum(System.readSeed());
	
	//constants for bead type
	std::vector<double> C;
	//defaults, 10.0=Umax, 0=Umin, rc=2, rmin=1, sigma=5.88
	for(int typeA=0;typeA<System.readNTypes();typeA++)
	{
		for(int typeB=0;typeB<System.readNTypes();typeB++)
		{
			//double rmin=1.0;
			//double rc=2.0*rmin;
			double UminNano=0.0;
			double UmaxNano=40.0;
			double sigma=5.88;
			
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
	
	//exceptions
	if(nNanoparticles>0)
	{
		int cindex=nBEADCONST*(BM+NANOTYPE*System.readNTypes());
		int cindexMirror=nBEADCONST*(NANOTYPE+BM*System.readNTypes());
		
		double UminNano=UminBrushNP;
		double UmaxNano=UmaxBrushNP;
		double sigma=5.88;
		
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
	std::vector< threeVector<double> > nanoPos;
	bool retryZ=true;
	int maxTrials=100000;
	for(auto zHeight=System.readSize().z;retryZ;zHeight+=1.0)
	{
		if(zHeight>1000)
		{ 
			std::cerr << "Too many retries, zHeight>1000!" << std::endl;
			return -1;
		}
		size=System.readSize();
		size.z=zHeight;
		System.setSize(size);
		std::cerr << "zHeight=" << zHeight << std::endl;
		int nTrials=0;
		int nanoOffset=System.readNParticles();
		for(int i=0;nanoPos.size()<nNanoparticles && nTrials<maxTrials;i++)
		{
			//std::cout << System.readSize().z << std::endl;
			//std::cerr << "Placing nanoparticle " << nanoPos.size() << " of " << nNanoparticles << "!" << std::endl;
			threeVector<double> toPlace;
			bool overlap;
			do
			{
				overlap=false;
				toPlace.x=(System.readSize().x-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				toPlace.y=(System.readSize().y-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				toPlace.z=static_cast<double>(nMonomers+2)*bondLength+1.0+
				(System.readSize().z-bondLength*(nMonomers+2)-2.0*nanoRadius)*randNum.rand53()+nanoRadius;
				//std::cout << toPlace.z << ' ' << System.readSize().z << std::endl;
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
					
					if(sqrt(d.x*d.x+d.y*d.y+d.z*d.z)<2.0*nanoRadius+0.5)
						overlap=true;
				}
				nTrials++;
				//std::cout << "Trial: " << nTrials << std::endl;
			} while(overlap && nTrials<maxTrials);
			//std::cerr << i << '\t' << toPlace.x << '\t' << toPlace.y << '\t' << toPlace.z << std::endl;
			nanoPos.push_back(toPlace);
		
		}
		if(nTrials>=maxTrials)
		{
			nanoPos.clear();
			retryZ=true;
		} else retryZ=false;
	}
	if(nNanoparticles>0)
		continuumSphere<double>(System, &(nanoPos[0]), nanoPos.size(), nanoRadius, &C[0], 
					System.readNTypes()*System.readNTypes()*nBEADCONST, NANOTYPE);
	
	molecule< double, fourVector<int> > boundaryZ;
	molecule< double, fourVector<int> > boundaryZ2;
	boundaryZ.setType(BOUNDARY);
	boundaryZ2.setType(OFFSET_BOUNDARY);
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		if(System.getPositions()[i].type==NANOTYPE)
			boundaryZ2.addBond(buf);
		else
			boundaryZ.addBond(buf);
	}
	boundaryZ.addConstant(2);//along Z, direction
	boundaryZ.addConstant(0.0);//at 0=z position, z_0
	boundaryZ.addConstant(4);//new:Doesn't matter old:half width of z, cutoff
	boundaryZ.addConstant(200);//k, epsilon
	
	boundaryZ2.addConstant(2);//along Z, direction
	boundaryZ2.addConstant(0.0);//at 0=z position, z_0
	boundaryZ2.addConstant(nanoRadius*2);//new:size of nanoparticle
	boundaryZ2.addConstant(200);//k, epsilon
	
	if(boundaryZ.readNBond()>0)
		System.addMolecule(boundaryZ);
	if(boundaryZ2.readNBond()>0)
		System.addMolecule(boundaryZ2);
	
	size=System.readSize();
	//This doubles the size, this is why it was always bigger
	//size.z*=2.0;
	//Instead, let us just add a smaller buffer of 2 times the NP radius
	size.z+=2*nanoRadius;
	System.setSize(size);
	
	/*
	boundaryZ.setType(BOUNDARY);
	for(int i=0;i<System.readNParticles();i++)
	{
		fourVector<int> buf;
		buf.s[0]=i;
		boundaryZ.addBond(buf);
	}
	boundaryZ.addConstant(2);//along Z
	boundaryZ.addConstant(size.z-6.0);//at z position
	boundaryZ.addConstant(3);//half width of z, cutoff
	boundaryZ.addConstant(400);//k
	
	System.addMolecule(boundaryZ);
	*/
	std::cout << "Storing configuration...";
	
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
	
	std::cout << "Done.\nExiting...\n";
	
	delete constants;
	return 0;
}


