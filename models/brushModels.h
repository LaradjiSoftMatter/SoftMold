/**
 * \brief Various model geometries associated with lipid systems.
 * We have flat and spherical cytoskeletons, solvent fills, and flat and spherical lipid bilayers. 
 * Many options for connectivity as well. Assumes your system blob has some functions to utilize 
 * with molecules, particles, and temperature.
 */

#include "../include/algorithms/functions.h"

//aspectRatio is x/y
template <typename T>
bool brush(Blob<T> &System, int nBrushes, int brushLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants, T aspectRatio, int bottomType, int chainType)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type (brush)...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type (brush)...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4![bool brush()]\n";
		return true;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=0;//nBrushes;
		bond.s[CHAINLENGTH]=brushLength;
		//m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nBrushes*brushLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
		
	
	//adjust size of system if it is out of bounds
	threeVector<T> s;
	
	s.x=sqrt((T)nBrushes/arealDensity)*sqrt(aspectRatio);
	s.y=sqrt((T)nBrushes/arealDensity)/sqrt(aspectRatio);
	s.z=pos.z;
	std::cerr << "Prefered system size is (brush): " << s.x << '\t' << s.y << '\t' << s.z << '\n';
	
	threeVector<T> size=System.readSize();
	
	std::cerr << "Current system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	//Adjust number of brushes to match grafting density
	
	size.x=(size.x<s.x)?s.x:size.x;
	size.y=(size.y<s.y)?s.y:size.y;
	size.z=(size.z<s.z)?s.z:size.z;
	
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+2*nBrushes*brushLength);
	
	threeVector<int> latticeSize;
	threeVector<T> latticeLength;
	
	latticeSize.x=static_cast<int>(sqrt((arealDensity*size.x*size.x)/aspectRatio));
	latticeSize.y=static_cast<int>(sqrt((arealDensity*size.y*size.y)*aspectRatio));
	
	int mBrushes=latticeSize.x*latticeSize.y;
	
	size.x=sqrt((T)mBrushes/arealDensity)*sqrt(aspectRatio);
	size.y=sqrt((T)mBrushes/arealDensity)/sqrt(aspectRatio);
	
	latticeLength.x=size.x/static_cast<T>(latticeSize.x);
	latticeLength.y=size.y/static_cast<T>(latticeSize.y);
		
	System.setSize(size);
	
	std::cerr << "Adjusted system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';	
	std::cerr << "Adjusted density is (brush): " << static_cast<T>(latticeSize.x*latticeSize.y)/(size.x*size.y) << '=';
	std::cerr << latticeSize.x << '*' << latticeSize.y << "/(" << size.x << '*' << size.y << ')' << std::endl;
	//latticeLength.x=size.x/static_cast<double>(latticeSize.x);
	//latticeLength.y=size.y/static_cast<double>(latticeSize.y);
	
	//initialize brush positions
	for(int i=0;i<latticeSize.x;i++)
	{
		for(int j=0;j<latticeSize.y;j++)
		{
			int brushIndex=(j*latticeSize.x+i);
			//if(brushIndex%1000==0)
			//	std::cerr << brushIndex << std::endl;
			T theta,phi;
			
			position<T> p;
			threeVector<T> v;
			threeVector<T> a;
			//set all acceleration to 0 initially
			a.x=0;
			a.y=0;
			a.z=0;
			
			for(int k=0;k<brushLength;k++)
			{
				//position
				p.x=i*latticeLength.x+0.001+size.x/2.0;
				p.y=j*latticeLength.y+0.001+size.y/2.0;
				p.z=pos.z+bondLength*(k+1);
				p.type=(k==0)?bottomType:chainType;
				
				//velocity
				theta=M_PI*randNum->rand53();
				phi=M_PI*2*randNum->rand53();
				v.x=Vrms*cos(phi)*sin(theta);
				v.y=Vrms*sin(phi)*sin(theta);
				v.z=Vrms*cos(theta);
				
				//save it
				System.addParticle(p,v,a);
				
				//if it is of BOND type and not the last monomer
				if(nConstants==2 && k!=brushLength-1)
				{
					//bond
					bond.s[0]=System.readNParticles()-1;
					bond.s[1]=System.readNParticles();
					m.addBond(bond);
				}
			}
			if(nConstants==4)
			{
				bond.s[NCHAINS]++;
			}
		}
	}	
	if(nConstants==4)
	{
		m.addBond(bond);
	}
	System.addMolecule(m);
	
	return true;
}

template <typename T>
std::vector<threeVector<T>> relaxSubstrate(std::vector<threeVector<T>> p, T radius, 
			threeVector<T> size, int iterations, int seed)
{
	MTRand randNum(seed);
	T radiusSqr=radius*radius;
	for(int it=0;it<iterations;it++)
	{
		for(int i=0;i<p.size();i++)
		{
			T theta=randNum.rand53()*2.0*M_PI;
			T r=randNum.rand53()*radius;
			auto temp=p[i];
			temp.x+=r*std::cos(theta);
			temp.y+=r*std::sin(theta);
			bool overlap=false;
			for(int j=0;j<p.size();j++)
			{
				if(i!=j)
				{
					threeVector<T> d;
					d.x=temp.x-p[j].x;
					d.y=temp.y-p[j].y;
					d.z=temp.z-p[j].z;
					d.x-=(d.x>=size.x/2.0)?size.x:0.0;
					d.y-=(d.y>=size.y/2.0)?size.y:0.0;
					d.z-=(d.z>=size.z/2.0)?size.z:0.0;
					d.x+=(d.x<-size.x/2.0)?size.x:0.0;
					d.y+=(d.y<-size.y/2.0)?size.y:0.0;
					d.z+=(d.z<-size.z/2.0)?size.z:0.0;
					
					if(d.x*d.x+d.y*d.y+d.z*d.z<radiusSqr)
						overlap=true;
				}
			}
			if(!overlap)
			{
				while(temp.x>=size.x)temp.x-=size.x;
				while(temp.y>=size.y)temp.y-=size.y;
				while(temp.x<0.0)temp.x+=size.x;
				while(temp.y<0.0)temp.y+=size.y;
				p[i]=temp;
			}
		}
		if(it%100==0)
			std::cerr << it << std::endl;
	}
	return p;
}

template <typename T>
std::vector<threeVector<T>> relaxSubstrate(std::vector<threeVector<T>> p, T radius, 
			threeVector<T> size, int iterations, int seed, T k, T temperature)
{
	T radiusSqr=radius*radius;
	auto pot= [radiusSqr,radius,k](const threeVector<T> &d)
	{
		double potential=0;
		T dist=d.x*d.x+d.y*d.y+d.z*d.z;
		if(dist<radiusSqr)
		{
			dist=std::sqrt(dist);
			potential+=0.5*k*(dist-radius)*(dist-radius);
		}
		return potential;
	};
	MTRand randNum(seed);
	for(int it=0;it<iterations;it++)
	{
		int accept=0,reject=0;
		for(int i=0;i<p.size();i++)
		{
			double potentialA=0,potentialB=0;
			auto temp=p[i];
			for(int j=0;j<p.size();j++)
			{
				if(i!=j)
				{
					threeVector<T> d;
					d.x=temp.x-p[j].x;
					d.y=temp.y-p[j].y;
					d.z=temp.z-p[j].z;
					d.x-=(d.x>=size.x/2.0)?size.x:0.0;
					d.y-=(d.y>=size.y/2.0)?size.y:0.0;
					d.z-=(d.z>=size.z/2.0)?size.z:0.0;
					d.x+=(d.x<-size.x/2.0)?size.x:0.0;
					d.y+=(d.y<-size.y/2.0)?size.y:0.0;
					d.z+=(d.z<-size.z/2.0)?size.z:0.0;
					
					potentialA+=pot(d);
				}
			}
			T theta=randNum.rand53()*2.0*M_PI;
			T r=randNum.rand53();
			temp.x+=r*std::cos(theta);
			temp.y+=r*std::sin(theta);
			for(int j=0;j<p.size();j++)
			{
				if(i!=j)
				{
					threeVector<T> d;
					d.x=temp.x-p[j].x;
					d.y=temp.y-p[j].y;
					d.z=temp.z-p[j].z;
					d.x-=(d.x>=size.x/2.0)?size.x:0.0;
					d.y-=(d.y>=size.y/2.0)?size.y:0.0;
					d.z-=(d.z>=size.z/2.0)?size.z:0.0;
					d.x+=(d.x<-size.x/2.0)?size.x:0.0;
					d.y+=(d.y<-size.y/2.0)?size.y:0.0;
					d.z+=(d.z<-size.z/2.0)?size.z:0.0;
					
					potentialB+=pot(d);
				}
			}
			//leave temperature around 1.0 for ~85 accept
			double D=std::exp((potentialA-potentialB)/temperature);
			
			if(D>randNum.rand53() || potentialB<=potentialA)
			//if(potentialB<=potentialA)
			{
				while(temp.x>=size.x)temp.x-=size.x;
				while(temp.y>=size.y)temp.y-=size.y;
				while(temp.x<0.0)temp.x+=size.x;
				while(temp.y<0.0)temp.y+=size.y;
				p[i]=temp;
				accept++;
			}
			else
				reject++;
		}
		if(it%100==0)
		{
			std::cerr << it << ' ' << accept << ' ' << reject << std::endl;
			for(int i=0;i<p.size();i++)
				for(int j=i+1;j<p.size();j++)
				{
					threeVector<T> d;
					d.x=p[i].x-p[j].x;
					d.y=p[i].y-p[j].y;
					d.z=p[i].z-p[j].z;
					d.x-=(d.x>=size.x/2.0)?size.x:0.0;
					d.y-=(d.y>=size.y/2.0)?size.y:0.0;
					d.z-=(d.z>=size.z/2.0)?size.z:0.0;
					d.x+=(d.x<-size.x/2.0)?size.x:0.0;
					d.y+=(d.y<-size.y/2.0)?size.y:0.0;
					d.z+=(d.z<-size.z/2.0)?size.z:0.0;
					if(d.x*d.x+d.y*d.y+d.z*d.z<0.00001)
						std::cerr << "Overlap?" << std::endl;
				}	
			//std::cout << p.size() << "\ntest\n";
			//for(auto &pi:p)
			//	std::cout << "1\t" << pi.x << '\t' << pi.y << '\t' << pi.z << '\n';
		}
	}
	return p;
}

//aspectRatio is x/y
template <typename T>
bool brush(Blob<T> &System, int nBrushes, int brushLength, threeVector<T> pos, T bondLength, T arealDensity, T *constants, int nConstants, T aspectRatio, int bottomType, int chainType, T k, T temperature)
{
	//Error checking
	if(nConstants==4)
		std::cerr << "nConstants is 4, assuming CHAIN type (brush)...\n";
	if(nConstants==2)
		std::cerr << "nConstants is 2, assuming BOND type (brush)...\n";
	
	if(nConstants!=2 && nConstants!=4)
	{
		std::cerr << "Error: nConstants is not 2 or 4![bool brush()]\n";
		return true;
	}
	
	//Initialize molecule
	fourVector<int> bond;
	molecule<T,fourVector<int> > m;
	
	if(nConstants==4)
	{
		m.setType(CHAIN);
		
		//bond for a chain type
		bond.s[START]=System.readNParticles();
		bond.s[NCHAINS]=0;//nBrushes;
		bond.s[CHAINLENGTH]=brushLength;
		//m.addBond(bond);
	}
	
	if(nConstants==2)
	{
		m.setType(BOND);
		m.allocBonds(nBrushes*brushLength);
	}
	
	for(int i=0;i<nConstants;i++)
		m.addConstant(constants[i]);
	
		
	
	//adjust size of system if it is out of bounds
	threeVector<T> s;
	
	s.x=sqrt((T)nBrushes/arealDensity)*sqrt(aspectRatio);
	s.y=sqrt((T)nBrushes/arealDensity)/sqrt(aspectRatio);
	s.z=pos.z;
	std::cerr << "Prefered system size is (brush): " << s.x << '\t' << s.y << '\t' << s.z << '\n';
	
	threeVector<T> size=System.readSize();
	
	std::cerr << "Current system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	//Adjust number of brushes to match grafting density
	
	size.x=s.x;//(size.x<s.x)?s.x:size.x;
	size.y=s.y;//(size.y<s.y)?s.y:size.y;
	size.z=(size.z<s.z)?s.z:size.z;
	//size=s;
	
	
	MTRand *randNum=new MTRand(System.readSeed());
	//Velocities
	T Vrms=sqrt(3.0*System.readInitialTemp());
	
	//This is important to reduce excessive reallocations
	System.allocParticle(System.readNParticles()+2*nBrushes*brushLength);
	
	threeVector<int> latticeSize;
	threeVector<T> latticeLength;
	
	latticeSize.x=static_cast<int>(sqrt((arealDensity*size.x*size.x)/aspectRatio))+1;
	latticeSize.y=static_cast<int>(sqrt((arealDensity*size.y*size.y)*aspectRatio))+1;
	
	//int mBrushes=latticeSize.x*latticeSize.y;
	
	//size.x=sqrt((T)mBrushes/arealDensity)*sqrt(aspectRatio);
	//size.y=sqrt((T)mBrushes/arealDensity)/sqrt(aspectRatio);
	
	latticeLength.x=size.x/static_cast<T>(latticeSize.x);
	latticeLength.y=size.y/static_cast<T>(latticeSize.y);
		
	System.setSize(size);
	
	std::cerr << "Adjusted system size is (brush): " << size.x << '\t' << size.y << '\t' << size.z << '\n';	
	std::cerr << "Adjusted density is (brush): " << static_cast<T>(nBrushes)/(size.x*size.y) << '=';
	std::cerr << nBrushes << "/(" << size.x << '*' << size.y << ')' << std::endl;
	//latticeLength.x=size.x/static_cast<double>(latticeSize.x);
	//latticeLength.y=size.y/static_cast<double>(latticeSize.y);
	
	//initialize substrate positions
	std::vector<threeVector<double>> substrate;
	for(double x=0;x<latticeSize.x;x++)
	{
		for(double y=0;y<latticeSize.y;y++)
		{
			if(substrate.size()<nBrushes)
			{
				threeVector<double> temp;
				temp.x=x*latticeLength.x;temp.y=y*latticeLength.y;temp.z=pos.z;
				substrate.push_back(temp);
			}
		}
	}
	
	T rSub=(latticeLength.x<latticeLength.y)?latticeLength.y:latticeLength.x;
	//randomly space substrate using metropolis hasting monte carlo
	substrate=relaxSubstrate(substrate,rSub*1.0,size,1000,System.readSeed(),k, temperature);
	
	//set up brushes with substrate
	for(auto &sub:substrate)
	{		
		T theta,phi;
		
		position<T> p;
		threeVector<T> v;
		threeVector<T> a;
		//set all acceleration to 0 initially
		a.x=0;
		a.y=0;
		a.z=0;
		for(int k=0;k<brushLength;k++)
		{
			//position
			p.x=sub.x;
			p.y=sub.y;
			p.z=sub.z+bondLength*k;
			p.type=(k==0)?bottomType:chainType;
			
			//velocity
			theta=M_PI*randNum->rand53();
			phi=M_PI*2*randNum->rand53();
			v.x=Vrms*cos(phi)*sin(theta);
			v.y=Vrms*sin(phi)*sin(theta);
			v.z=Vrms*cos(theta);
					
			//save it
			System.addParticle(p,v,a);
			
			//if it is of BOND type and not the last monomer
			if(nConstants==2 && k!=brushLength-1)
			{
				//bond
				bond.s[0]=System.readNParticles()-1;
				bond.s[1]=System.readNParticles();
				m.addBond(bond);
			}
		}
		if(nConstants==4)
		{
			bond.s[NCHAINS]++;
		}
	}
	if(nConstants==4)
	{
		m.addBond(bond);
	}
	System.addMolecule(m);
	
	return true;
}
