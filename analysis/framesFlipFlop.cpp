/** \brief Outputs molecules with nearby orientations per frame. Effectively separating lipids per layer.
 */

//For molecular dynamics forces and potentials
#include "../include/MD.h"

//For the molecular dynamics variables
#include "../include/system.h"

#define MIN_NUMBER_SURFACE_LIPIDS 50

int main(int argc, char* argv[])
{
	if(argc!=6 && argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cout << "Takes frame and molecule data to generate maps of flipping events.\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex cutoff minTheta frameSkip\n";
		std::cout << "\tusage: " << argv[0] << " name molIndex cutoff minTheta\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int molIndex,frameSkip=1;
	double cutoff, minTheta;
	std::stringstream cmdInp;
	if(argc==6)
	{
		cmdInp << argv[2] << '\t' << argv[3] << '\t' << argv[4] << '\t' << argv[5] << '\n';
		cmdInp >> molIndex >> cutoff >> minTheta >> frameSkip;
	}
	if(argc==5)
	{
		cmdInp << argv[2] << '\t' << argv[3] << '\t' << argv[4] << '\n';
		cmdInp >> molIndex >> cutoff >> minTheta >> frameSkip;
	}
	
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
	
	//initialize surface pool using first configuration
	xyzFile.load();
	
	//molecule center of masses and orientations
	position<double> *molCom=new position<double>[molNChains];
	threeVector<double> *orient=new threeVector<double>[molNChains];
	
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
	Cell<double> cellSearch(molCom, molNChains, cutoff, size);
	
	cellSearch.build();
	
	int nSurfaces=0;
	
	//set up surfaces
	for(int i=0;i<molNChains;i++)
	{
		if(molCom[i].type==-1)
		{
			std::vector<int> stack;
			stack.push_back(i);
			molCom[i].type=nSurfaces;
			
			while(stack.size()>0)
			{
				int currentMol=stack.back();
				stack.pop_back();
				threeVector<double> minImg;
				
				//This section will be changed once I verify Cell->query()
				for(int j=cellSearch.query(currentMol,minImg);j!=-1;j=cellSearch.query(currentMol,minImg))
				{
					//check if chain is within cutoff, of the same orientation, and unflagged
					if(molCom[j].type==-1)//flagged?
					{
						threeVector<double> d;
						minImg=0;
						
						d.x=molCom[j].x-molCom[currentMol].x;
						d.y=molCom[j].y-molCom[currentMol].y;
						d.z=molCom[j].z-molCom[currentMol].z;
						
						//set minimum image
						//God this looks stupid when query does this already!
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
						{
							if(dotProduct(orient[currentMol],orient[j])>minCosTheta)//is it oriented?
							{
								//flag it and add it to our stack
								molCom[j].type=nSurfaces;
								stack.push_back(j);
							}
						}
					}
				}
			}
			
			nSurfaces++;
		}
	}
	
	//To correlate and consolidate surfaces
	
	std::vector<int> nSurfaceMolOld(molNChains,0);//number of old surface molecules
	std::vector<int> oldMap(molNChains,0);//map to correlate between surfaces
	
	//count how many exist on a surface
	for(int i=0;i<molNChains;i++)
		nSurfaceMolOld[molCom[i].type]++;
	
	cellSearch.resetState();
	
	//consolidate strange surfaces (i.e. 1 molecule only surfaces)
	for(int i=0;i<molNChains;i++)
	{
		//is it strange?
		if(nSurfaceMolOld[molCom[i].type]<MIN_NUMBER_SURFACE_LIPIDS)
		{
			threeVector<double> minImg;
			double nearestD=cutoffSqr;
			int nearestI=-1;
			//merge it with whatever the nearest surface is, regardless of orientation difference
			for(int j=cellSearch.query(i,minImg);j!=-1;j=cellSearch.query(i,minImg))
			{
				threeVector<double> d;
				
				d.x=molCom[j].x-molCom[i].x;
				d.y=molCom[j].y-molCom[i].y;
				d.z=molCom[j].z-molCom[i].z;
				
				d.x+=minImg.x;
				d.y+=minImg.y;
				d.z+=minImg.z;
				
				int currentD=d.x*d.x+d.y*d.y+d.z*d.z;
				if(currentD<nearestD && nSurfaceMolOld[molCom[j].type]>MIN_NUMBER_SURFACE_LIPIDS)//is it closer than anything so far?
				{
					nearestD=currentD;
					nearestI=j;
				}
			}
			
			//remove the previous designation
			nSurfaceMolOld[molCom[i].type]--;
			
			//Use the nearest index's surface
			molCom[i].type=molCom[nearestI].type;
			
			//add the current designation
			nSurfaceMolOld[molCom[i].type]++;
		}
	}
	
	//remap and consolidate types
	std::vector<int> remapTypes(molNChains,-1);
	int nSurfacesCorrected=0;
	for(int i=0;i<nSurfaces;i++)
	{
		if(nSurfaceMolOld[i]>0)
		{
			remapTypes[i]=nSurfacesCorrected;
			nSurfaceMolOld[nSurfacesCorrected]=nSurfaceMolOld[i];
			nSurfacesCorrected++;
		}
	}
	
	//set the new map
	for(int i=0;i<molNChains;i++)
	{
		molCom[i].type=remapTypes[molCom[i].type];
		oldMap[i]=molCom[i].type;
	}
	
	std::cout << 0;
	for(int i=0;i<nSurfacesCorrected;i++)
		if(nSurfaceMolOld[i]>0)
			std::cout << '\t' << nSurfaceMolOld[i];
		
	std::cout << '\n';
	
	int frame=0;
	double time=0;
	while(xyzFile.load())
	{
		frame++;
		time+=System.readStoreInterval();
		if(frame%frameSkip==0)
		{
			//make sure the size hasn't changed
			if(sizeFile.is_open())
			{
				double sTime=0;
				threeVector<double> sCurrent;
				//for the first iteration, I expect this to exit the loop after one read,
				// otherwise, this is a template for subsequent reads in the loop that is
				// coming up.
				while(sTime<time && !sizeFile.eof())
					sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
				size=sCurrent;
			}
			
			//output surface configura*tion using old surface values
			std::fstream nextFileFrame;
			std::stringstream nextFileName;
			nextFileName << "flipFlop_" << name << "_" << time << ".xyz";
			nextFileFrame.open(nextFileName.str().c_str(),std::ios::out);
			
			int nT=System.readNTypes();
			nextFileFrame << molNChains*molLength << "\nasdf\n";
			for(int i=0;i<molNChains;i++)
				for(int j=molStart+i*molLength;j<molStart+(i+1)*molLength;j++)
					nextFileFrame << (molCom[i].type*nT+p[j].type) << '\t' << p[j].x << '\t' << p[j].y << '\t' << p[j].z << '\n';
				nextFileFrame << '\n';
			nextFileFrame.close();
			
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
			
			nSurfaces=0;
			
			//set up surfaces
			for(int i=0;i<molNChains;i++)
			{
				if(molCom[i].type==-1)
				{
					std::vector<int> stack;
					stack.push_back(i);
					molCom[i].type=nSurfaces;
					
					while(stack.size()>0)
					{
						int currentMol=stack.back();
						stack.pop_back();
						threeVector<double> minImg;
						
						//This section will be changed once I verify Cell->query()
						for(int j=newSearch.query(currentMol,minImg);j!=-1;j=newSearch.query(currentMol,minImg))
						{
							//check if chain is within cutoff, of the same orientation, and unflagged
							if(molCom[j].type==-1)//flagged?
							{
								threeVector<double> d;
								minImg=0;
								
								d.x=molCom[j].x-molCom[currentMol].x;
								d.y=molCom[j].y-molCom[currentMol].y;
								d.z=molCom[j].z-molCom[currentMol].z;
								
								//set minimum image
								//God this looks stupid when query does this already!
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
								{
									if(dotProduct(orient[currentMol],orient[j])>minCosTheta)//is it oriented?
									{
										//flag it and add it to our stack
										molCom[j].type=nSurfaces;
										stack.push_back(j);
									}
								}
							}
						}
					}
					
					nSurfaces++;
				}
			}
			
			std::vector<int> nSurfaceMol(molNChains,0);//current number of surface molecules
			
			//count how many exist on a surface
			for(int i=0;i<molNChains;i++)
				nSurfaceMol[molCom[i].type]++;
			
			newSearch.resetState();
			
			//consolidate strange surfaces (i.e. 1 molecule only surfaces)
			for(int i=0;i<molNChains;i++)
			{
				//is it strange?
				if(nSurfaceMol[molCom[i].type]<MIN_NUMBER_SURFACE_LIPIDS)
				{
					threeVector<double> minImg;
					double nearestD=cutoffSqr;
					int nearestI=-1;
					//merge it with whatever the nearest surface is, regardless of orientation difference
					for(int j=newSearch.query(i,minImg);j!=-1;j=newSearch.query(i,minImg))
					{
						threeVector<double> d;
						
						d.x=molCom[j].x-molCom[i].x;
						d.y=molCom[j].y-molCom[i].y;
						d.z=molCom[j].z-molCom[i].z;
						
						d.x+=minImg.x;
						d.y+=minImg.y;
						d.z+=minImg.z;
						
						int currentD=d.x*d.x+d.y*d.y+d.z*d.z;
						if(currentD<nearestD && nSurfaceMol[molCom[j].type]>MIN_NUMBER_SURFACE_LIPIDS)//is it closer than anything so far?
						{
							nearestD=currentD;
							nearestI=j;
						}
					}
					
					//remove the previous designation
					nSurfaceMol[molCom[i].type]--;
					
					if(nearestI==-1)
					{
						std::cout << "Warning(main): consolidation search couldn't resolve a nearby molecule.\n";
						std::cout << "\tTry increasing cutoff. Defaulting to index 0.\n";
						nearestI=0;
					}
					
					//Use the nearest index's surface
					molCom[i].type=molCom[nearestI].type;
					
					//add the current designation
					nSurfaceMol[molCom[i].type]++;
				}
			}
			
			nSurfacesCorrected=0;
			for(int i=0;i<nSurfaces;i++)
			{
				if(nSurfaceMol[i]>0)
				{
					remapTypes[i]=nSurfacesCorrected;
					nSurfaceMol[nSurfacesCorrected]=nSurfaceMol[i];
					nSurfacesCorrected++;
				}
			}
			
			for(int i=0;i<molNChains;i++)
			{
				molCom[i].type=remapTypes[molCom[i].type];
			}
			
			//std::cout << time << '\t' << nSurfacesCorrected << '\t' << nSurfaces << '\t';
			/*
			std::vector<int> correlationIndex(molNChains,-1);
			 
			//correlate old and new surfaces
			for(int surfaceIndex=0;surfaceIndex<nSurfacesCorrected;surfaceIndex++)
			{
				if(nSurfaceMol[surfaceIndex]>0)
				{
					//correlation count between old and new surfaces
					std::vector<int> nSurfaceMolCorr(molNChains,0);
					 
					//count the correlations to the old map
					for(int i=0;i<molNChains;i++)
					{
						if(molCom[i].type==surfaceIndex)
							nSurfaceMolCorr[oldMap[i]]++;
					}
					int maxCorrelation=0;
					for(int i=0;i<molNChains;i++)
						if(nSurfaceMolCorr[i]>nSurfaceMolCorr[maxCorrelation])
							maxCorrelation=i;
				
					correlationIndex[maxCorrelation]=surfaceIndex;
				}
			}
			*/
			//if(nSurfacesCorrected==2)
			{
				std::cout << time;
				for(int i=0;i<nSurfacesCorrected;i++)
					if(nSurfaceMol[i]>0)
						std::cout << '\t' << nSurfaceMol[i];
					//if(nSurfaceMol[correlationIndex[i]]>0)
					//					std::cout << '\t' << nSurfaceMol[correlationIndex[i]];
					std::cout << '\n';
			}
			for(int i=0;i<molNChains;i++)
				oldMap[i]=molCom[i].type;
		}
	}
	
	if(sizeFile.is_open())
		sizeFile.close();
	delete molCom, orient;
	
	return 0;
}

