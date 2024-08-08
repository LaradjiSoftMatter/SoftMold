/** \brief Outputs CHAIN type molecules per frame and sets their type by how many lipids they have.
 */

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#define MIN_NUMBER_SURFACE_LIPIDS 50

int main(int argc, char* argv[])
{
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Takes frame and molecule data to generate maps of density.\n";
		std::cout << "usage: " << argv[0] << " name molIndex cutoff minTheta\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int molIndex;
	double cutoff, minTheta;
	std::stringstream cmdInp;
	cmdInp << argv[2] << '\t' << argv[3] << '\t' << argv[4] << '\n';
	cmdInp >> molIndex >> cutoff >> minTheta;
	
	double minCosTheta=std::cos(minTheta);
	double cutoffSqr=cutoff*cutoff;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	molecule<double, fourVector<int> > *m=System.getMolecule();
	
	//our current size
	threeVector<double> size=System.readSize();
	
	//extra error checking
	if(molIndex<System.readNMolecules() && molIndex>=0)
	{
		if(m[molIndex].readType()!=CHAIN)
		{
			std::cout << "Error (main): molIndex, " << molIndex << ", not CHAIN, " << CHAIN << ", type!\n";
			return 0;
		}
	}
	else
	{
		std::cout << "Error (main): molIndex, " << molIndex << ", out of bounds [0," << System.readNMolecules()-1 << "].\n";
		return 0;
	}
	
	//check if size varies
	std::string newName("size_");
	newName+=name;
	newName+=".dat";
	
	std::fstream sizeFile;
	sizeFile.open(newName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	std::cout << 0 << '\t' << size.x << '\t' << size.y << '\t' << size.z << '\n';
	
	int molStart=m[molIndex].getBonds()[0].s[START];
	int molLength=m[molIndex].getBonds()[0].s[CHAINLENGTH];
	int molNChains=m[molIndex].getBonds()[0].s[NCHAINS];
	
	//Get our frames file
	std::string framesName("frames_");
	framesName+=name;
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	//molecule center of masses and orientations
	position<double> *molCom=new position<double>[molNChains];
	threeVector<double> *orient=new threeVector<double>[molNChains];
	
	double time=0;
	while(xyzFile.load())
	{
		//make sure the size hasn't changed
		if(sizeFile.is_open())
		{
			double sTime=0;
			threeVector<double> sCurrent=size;//just in case it is current
			//for the first iteration, I expect this to exit the loop after one read,
			// otherwise, this is a template for subsequent reads in the loop that is
			// coming up.
			while(sTime<time && !sizeFile.eof())
			{
				sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
			size=sCurrent;
		}
		
		//our next surfaces
		//for each molecule
		for(int i=0;i<molNChains;i++)
		{
			//Get our orientation
			threeVector<double> endToEnd;
			int first=molStart+i*molLength;
			int last=molStart+(i+1)*molLength-1;
			
			if(p[first].type>p[last].type)
			{
				int temp=first;
				first=last;
				last=temp;
			}
			
			endToEnd.x=p[first].x-p[last].x;
			endToEnd.y=p[first].y-p[last].y;
			endToEnd.z=p[first].z-p[last].z;
			
			if(endToEnd.x>size.x/2.0)
				endToEnd.x-=size.x;
			if(endToEnd.y>size.y/2.0)
				endToEnd.y-=size.y;
			if(endToEnd.z>size.z/2.0)
				endToEnd.z-=size.z;
			
			if(endToEnd.x<-size.x/2.0)
				endToEnd.x+=size.x;
			if(endToEnd.y<-size.y/2.0)
				endToEnd.y+=size.y;
			if(endToEnd.z<-size.z/2.0)
				endToEnd.z+=size.z;
			
			orient[i]=unitVector(endToEnd);
			
			
			
			//Get our center of mass
			position<double> comVector;
			comVector.x=0;
			comVector.y=0;
			comVector.z=0;
			comVector.type=-1;//unflagged
			
			//every particle of a particular molecule
			for(int j=molStart+i*molLength;j<molStart+(i+1)*molLength;j++)
			{
				threeVector<double> minImg;
				minImg=0;
				
				//check with the first one for relative frame
				threeVector<double> d;
				
				d.x=p[j].x-p[molStart+i*molLength].x;
				d.y=p[j].y-p[molStart+i*molLength].y;
				d.z=p[j].z-p[molStart+i*molLength].z;
				
				//set minimum image
				//Also, I believe orthoganal ordering of branches should be optimal
				// (or at least the compiler/cpu should see it that way)
				if(d.x>size.x/2.0)
					minImg.x=-size.x;
				
				if(d.y>size.y/2.0)
					minImg.y=-size.y;
				
				if(d.z>size.z/2.0)
					minImg.z=-size.z;
				
				if(d.x<-size.x/2.0)
					minImg.x=size.x;
				
				if(d.y<-size.y/2.0)
					minImg.y=size.y;
				
				if(d.z<-size.z/2.0)
					minImg.z=size.z;
				
				
				//com is the accumulated average of all particles in the same frame
				comVector.x+=p[j].x+minImg.x;
				comVector.y+=p[j].y+minImg.y;
				comVector.z+=p[j].z+minImg.z;
			}
			
			//average center of mass
			comVector.x/=static_cast<double>(molLength);
			comVector.y/=static_cast<double>(molLength);
			comVector.z/=static_cast<double>(molLength);
			
			//make sure it is inside our system bounds 
			// (all particles and subsequent molecules center of masses should be)
			if(comVector.x<0)
				comVector.x+=size.x;
			if(comVector.y<0)
				comVector.y+=size.y;
			if(comVector.z<0)
				comVector.z+=size.z;
			if(comVector.x>=size.x)
				comVector.x-=size.x;
			if(comVector.y>=size.y)
				comVector.y-=size.y;
			if(comVector.z>=size.z)
				comVector.z-=size.z;
			
			//now we have our center of mass
			molCom[i]=comVector;
			//std::cout << 1 << '\t' << molCom[i].x << '\t' << molCom[i].y << '\t' << molCom[i].z << '\n';
		}
		
		//set up our cell lists
		Cell<double> newSearch(molCom, molNChains, cutoff, size);
		
		newSearch.build();
		
		//set up surfaces
		for(int i=0;i<molNChains;i++)
		{
			molCom[i].type=0;
			threeVector<double> minImg;
			
			//This section will be changed once I verify Cell->query()
			for(int j=newSearch.query(i,minImg);j!=-1;j=newSearch.query(i,minImg))
			{
				threeVector<double> d;
				minImg=0;
				
				d.x=molCom[j].x-molCom[i].x;
				d.y=molCom[j].y-molCom[i].y;
				d.z=molCom[j].z-molCom[i].z;
				
				//set minimum image, I forgot to take this out, and now I'm too lazy to do so
				//Also, I believe orthoganal ordering of branches should be optimal
				// (or at least the compiler/cpu should see it that way)
				if(d.x>size.x/2.0)
					minImg.x=-size.x;
				
				if(d.y>size.y/2.0)
					minImg.y=-size.y;
				
				if(d.z>size.z/2.0)
					minImg.z=-size.z;
				
				if(d.x<-size.x/2.0)
					minImg.x=size.x;
				
				if(d.y<-size.y/2.0)
					minImg.y=size.y;
				
				if(d.z<-size.z/2.0)
					minImg.z=size.z;
				
				
				d.x+=minImg.x;
				d.y+=minImg.y;
				d.z+=minImg.z;
				
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)//within cutoff?
					if(dotProduct(orient[i],orient[j])>minCosTheta)//is it oriented?
						molCom[i].type++;
			}
		}
		
		int minNeighbors=molNChains;
		for(int i=0;i<molNChains;i++)
			if(molCom[i].type<minNeighbors)
				minNeighbors=molCom[i].type;
		
		//output surface configuration using current surface values
		std::fstream nextFileFrame;
		std::stringstream nextFileName;
		nextFileName << "density_" << name << "_" << time << ".xyz";
		nextFileFrame.open(nextFileName.str().c_str(),std::ios::out);
		
		nextFileFrame << molNChains*molLength << "\nasdf\n";
		for(int i=0;i<molNChains;i++)
			for(int j=molStart+i*molLength;j<molStart+(i+1)*molLength;j++)
				nextFileFrame << molCom[i].type-minNeighbors << '\t' << p[j].x << '\t' << p[j].y << '\t' << p[j].z << '\n';
		nextFileFrame << '\n';
		nextFileFrame.close();
		time+=System.readStoreInterval();
	}
	
	if(sizeFile.is_open())
		sizeFile.close();
	delete molCom, orient;
	
	return 0;
}

