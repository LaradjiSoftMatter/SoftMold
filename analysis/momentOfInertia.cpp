/** \brief Calculates eigenvectors for moment of inertia tensor.
 */

#include "../include/MD.h"
#include "../include/system.h"
#include "../include/algorithms/eigen.h"
#include <ctime>

// | m[0][0] m[0][1] m[0][2] |
// | m[1][0] m[1][1] m[1][2] | = det(m) and det(m') [aka conDet(m)]
// | m[2][0] m[2][1] m[2][2] |
double det(std::vector<std::vector<double>> m)
{
	return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-
	       m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+ 
	       m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]); 
}

// sqrt(1.0-b^2/a^2)
double eccentricity(double b, double a)
{
	return sqrt(1.0-(b*b)/(a*a));
}

int main(int argc, char* argv[])
{
	if(argc!=5)
	{
		//so simple, this can't possibly mess it up
		std::cerr << "Calculates eigenvectors for moment of inertia tensor\n";
		std::cerr << "\tusage: " << argv[0] << " name type nanoType cutoff\n";
		return 0;
	}
	
	//read in options
	char *name=argv[1];
	int type,nanoType;
	std::stringstream cmdInp;
	double cutoff;
	cmdInp << argv[2] << ' ' << argv[3] << ' ' << argv[4];
	cmdInp >> type >> nanoType >> cutoff;
	
	//Configuration variables
	Blob<double> System;
	
	//load old variables, then initialize them
	Script<double, Blob <double> > fileIO(name,std::ios::in,&System);
	fileIO.read();
	fileIO.close();
	
	//Get our frames file
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	//This one is now handling allocation of the 'p' pointer.
	std::string newName("frames_");
	newName+=name;
	newName+=".xyz";
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(newName.c_str(), std::ios::in);
	
	std::vector<int> typeI, nanoI;
	for(int i=0;i<nParticles;i++)
	{
		if(p[i].type==type)
			typeI.push_back(i);
		if(p[i].type==nanoType)
			nanoI.push_back(i);
	}
	
	time_t start;
	std::time(&start);
	int frame=0;
	double time=0;
	
	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string sizeName("size_");
	sizeName+=argv[1];
	sizeName+=".dat";
		
	std::fstream sizeFile;
	sizeFile.open(sizeName.c_str(), std::ios::in);
	if(sizeFile.is_open())
	{
		double sTime=-1;
		threeVector<double> sCurrent;
		//for the first iteration, I expect this to exit the loop after one read,
		// otherwise, this is a template for subsequent reads in the loop that is
		// coming up.
		while(sTime<0 && !sizeFile.eof())
			sizeFile >> sTime >> sCurrent.x >> sCurrent.y >> sCurrent.z;
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
	
	std::vector<std::vector<double>> moiELast;
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
				//std::cout << time << '\t' << sTime << '\t' << size.x << '\t' << size.y << '\t' << size.z << std::endl;
				
			}
			size=sCurrent;
		}
		
		double cutoffSqr=cutoff*cutoff;
		
		//center of mass
		threeVector<double> com=0;
		int iCount=0;
		for(auto& i:typeI)
		{
			bool exclude=false;
			for(auto& j:nanoI)
			{
				threeVector<double> d=0;
				d.x=p[i].x-p[j].x;
				d.y=p[i].y-p[j].y;
				d.z=p[i].z-p[j].z;
				d.x+=(d.x>size.x/2.0)?-size.x:0;
				d.y+=(d.y>size.y/2.0)?-size.y:0;
				d.z+=(d.z>size.z/2.0)?-size.z:0;
				d.x+=(d.x<-size.x/2.0)?size.x:0;
				d.y+=(d.y<-size.y/2.0)?size.y:0;
				d.z+=(d.z<-size.z/2.0)?size.z:0;
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					exclude=true;
			}
			if(!exclude)
			{
				threeVector<double> minImg=0;
				minImg.x+=(p[i].x-p[0].x>size.x/2.0)?size.x:0;
				minImg.y+=(p[i].y-p[0].y>size.y/2.0)?size.y:0;
				minImg.z+=(p[i].z-p[0].z>size.z/2.0)?size.z:0;
				minImg.x+=(p[i].x-p[0].x<-size.x/2.0)?-size.x:0;
				minImg.y+=(p[i].y-p[0].y<-size.y/2.0)?-size.y:0;
				minImg.z+=(p[i].z-p[0].z<-size.z/2.0)?-size.z:0;
				com.x+=p[i].x+minImg.x;
				com.y+=p[i].y+minImg.y;
				com.z+=p[i].z+minImg.z;
				iCount++;
			}
		}
		com.x/=iCount;
		com.y/=iCount;
		com.z/=iCount;
		while(com.x>size.x)com.x-=size.x;
		while(com.y>size.y)com.y-=size.y;
		while(com.z>size.z)com.z-=size.z;
		while(com.x<0)com.x+=size.x;
		while(com.y<0)com.y+=size.y;
		while(com.z<0)com.z+=size.z;
		
		//our initial matrix
		std::vector<std::vector<double>> moi(3,std::vector<double>(3,0));
		for(auto& i:typeI)
		{
			bool exclude=false;
			for(auto& j:nanoI)
			{
				threeVector<double> d=0;
				d.x=p[i].x-p[j].x;
				d.y=p[i].y-p[j].y;
				d.z=p[i].z-p[j].z;
				d.x+=(d.x>size.x/2.0)?-size.x:0;
				d.y+=(d.y>size.y/2.0)?-size.y:0;
				d.z+=(d.z>size.z/2.0)?-size.z:0;
				d.x+=(d.x<-size.x/2.0)?size.x:0;
				d.y+=(d.y<-size.y/2.0)?size.y:0;
				d.z+=(d.z<-size.z/2.0)?size.z:0;
				if(d.x*d.x+d.y*d.y+d.z*d.z<cutoffSqr)
					exclude=true;
			}
			//if(!exclude)
			{
				threeVector<double> r;
				r.x=p[i].x-com.x;
				r.y=p[i].y-com.y;
				r.z=p[i].z-com.z;
				r.x+=(r.x>size.x/2.0)?-size.x:0;
				r.y+=(r.y>size.y/2.0)?-size.y:0;
				r.z+=(r.z>size.z/2.0)?-size.z:0;
				r.x+=(r.x<-size.x/2.0)?size.x:0;
				r.y+=(r.y<-size.y/2.0)?size.y:0;
				r.z+=(r.z<-size.z/2.0)?size.z:0;
				moi[0][0]+=r.y*r.y+r.z*r.z;
				moi[1][1]+=r.x*r.x+r.z*r.z;
				moi[2][2]+=r.x*r.x+r.y*r.y;
				moi[0][1]+=r.x*r.y;
				moi[1][2]+=r.z*r.y;
				moi[0][2]+=r.x*r.z;
			}
		}
		moi[1][0]=moi[0][1];
		moi[2][1]=moi[1][2];
		moi[2][0]=moi[0][2];
		//dumpMatrix(moi);
		//for(auto &i:moi)
		//	for(auto &j:i)
		//		j/=typeI.size();
		//eigenvectors and eigenvalues of the matrix
		auto eigenV=eigen(moi);
		std::vector<std::vector<double>> moiE=std::get<0>(eigenV);
		std::vector<double> eigenValues=std::get<1>(eigenV);
		int largest=0,smallest=1;
		for(int i=0;i<eigenValues.size();i++)
			if(eigenValues[i]>eigenValues[largest]) largest=i;
			else if(eigenValues[i]<eigenValues[smallest]) smallest=i;
		int middle;
		if((smallest==0 && largest==1) || (largest==0 && smallest==1)) middle=2;
		if((smallest==0 && largest==2) || (largest==0 && smallest==2)) middle=1;
		if((smallest==1 && largest==2) || (largest==1 && smallest==2)) middle=0;
		for(auto &row:moiE)
		{
			std::vector<double> newRow=row;
			newRow[0]=row[smallest];
			newRow[1]=row[middle];
			newRow[2]=row[largest];
			row=newRow;
		}
		moiELast=moiE;
		//if(detVal>0 && dotProd>0) then it is correct, no transform
		//if(detVal>0 && dotProd<0) then it is correct handedness but axis is flipped
		//if(detVal<0 && dotProd<0) then it is incorrect handedness and axis is correct
		//if(detVal<0 && dotProd>0) then it is incorrect handedness and axis is flipped
			
		double detVal=det(moiE);
		double zDotProd=moiE[2][2];//sqrt(moiE[0][2]*moiE[0][2]+moiE[1][2]*moiE[1][2]+moiE[2][2]*moiE[2][2]);

		//if(detVal>0 && dotProd>0) nothing
		if(detVal>0 && zDotProd<0)//rotate the y z plane around x axis
		{
			moiE[0][2]=-moiE[0][2];
			moiE[1][2]=-moiE[1][2];
			moiE[2][2]=-moiE[2][2];
			moiE[0][1]=-moiE[0][1];
			moiE[1][1]=-moiE[1][1];
			moiE[2][1]=-moiE[2][1];
		}
		if(detVal<0 && zDotProd<0)//flip handedness
		{
			for(auto &row:moiE)
				for(auto &col:row)
					col=-col;
		}
		if(detVal<0 && zDotProd>0)//flip handedness and rotate the y z plane around the x axis
		{
			for(auto &row:moiE)
				for(auto &col:row)
					col=-col;
			moiE[0][2]=-moiE[0][2];
			moiE[1][2]=-moiE[1][2];
			moiE[2][2]=-moiE[2][2];
			moiE[0][1]=-moiE[0][1];
			moiE[1][1]=-moiE[1][1];
			moiE[2][1]=-moiE[2][1];
		}
		
		double yDotProd=moiE[1][1];//sqrt(moiE[0][2]*moiE[0][2]+moiE[1][2]*moiE[1][2]+moiE[2][2]*moiE[2][2]);
		
		if(yDotProd<0) //rotate the x y plane about the z axis
		{
			moiE[0][0]=-moiE[0][0];
			moiE[1][0]=-moiE[1][0];
			moiE[2][0]=-moiE[2][0];
			moiE[0][1]=-moiE[0][1];
			moiE[1][1]=-moiE[1][1];
			moiE[2][1]=-moiE[2][1];
		}
		std::sort(eigenValues.begin(),eigenValues.end());
		std::cerr << time << ' ';
		for(auto e:eigenValues)
			std::cerr << e << ' ';
		std::cerr << det(moiE) << ' ' << moiE[0][0] << ' ' << moiE[1][1] << ' ' << moiE[2][2] << ' ';
		//std::cerr << std::endl;
		
		if(nanoI.size()>1)
		{
			threeVector<double> r;
			r.x=p[nanoI[0]].x-com.x;//p[nanoI[1]].x;
			r.y=p[nanoI[0]].y-com.y;//p[nanoI[1]].y;
			r.z=p[nanoI[0]].z-com.z;//p[nanoI[1]].z;
			r.x+=(r.x>size.x/2.0)?-size.x:0;
			r.y+=(r.y>size.y/2.0)?-size.y:0;
			r.z+=(r.z>size.z/2.0)?-size.z:0;
			r.x+=(r.x<-size.x/2.0)?size.x:0;
			r.y+=(r.y<-size.y/2.0)?size.y:0;
			r.z+=(r.z<-size.z/2.0)?size.z:0;
			r/=magnitude(r);
			
			//std::cout << "7\t" << com.x+r.x*20 << "\t" << com.y+r.y*20 << "\t" << com.z+r.z*20 << std::endl;
			//double dot0=moiE[0][0]
			
			//using the NP as the plane normal:
				double dot1=(r.x*moiE[0][0]+r.y*moiE[1][0]+r.z*moiE[2][0]);
				double dot2=(r.x*moiE[0][1]+r.y*moiE[1][1]+r.z*moiE[2][1]);
				double dot3=(r.x*moiE[0][2]+r.y*moiE[1][2]+r.z*moiE[2][2]);
				double eccent=0;
			//if(std::abs(dot1)>std::abs(dot2) && std::abs(dot2)>std::abs(dot3)) eccent=eccentricity(std::abs(dot3),std::abs(dot2));
			//if(std::abs(dot2)>std::abs(dot1) && std::abs(dot1)>std::abs(dot3)) eccent=eccentricity(std::abs(dot3),std::abs(dot1));
			//if(std::abs(dot3)>std::abs(dot1) && std::abs(dot1)>std::abs(dot2)) eccent=eccentricity(std::abs(dot2),std::abs(dot1));
			//if(std::abs(dot1)>std::abs(dot3) && std::abs(dot3)>std::abs(dot2)) eccent=eccentricity(std::abs(dot2),std::abs(dot3));
			//if(std::abs(dot2)>std::abs(dot3) && std::abs(dot3)>std::abs(dot1)) eccent=eccentricity(std::abs(dot1),std::abs(dot3));
			//if(std::abs(dot3)>std::abs(dot2) && std::abs(dot2)>std::abs(dot1)) eccent=eccentricity(std::abs(dot1),std::abs(dot2));
			//double dot1=r.x*moiE[0][0]+r.y*moiE[0][1]+r.z*moiE[0][2];
			//double dot2=r.x*moiE[1][0]+r.y*moiE[1][1]+r.z*moiE[1][2];
			//double dot3=r.x*moiE[2][0]+r.y*moiE[2][1]+r.z*moiE[2][2];
			
			if(std::abs(eigenValues[0]-eigenValues[1])<std::abs(eigenValues[1]-eigenValues[2]))//prolate
			{
				eccent=eccentricity(sqrt((eigenValues[0]+eigenValues[1])/2.0),sqrt(eigenValues[2]));
				std::cerr << "prolate ";
			}
			else//oblate
			{
				eccent=eccentricity(sqrt(eigenValues[0]),sqrt((eigenValues[1]+eigenValues[2])/2.0));
				std::cerr << "oblate ";
			}
			
			std::cerr << r.x << ' ' << r.y << ' ' << r.z << std::endl;
			std::cout << time << ' ' << eccent << ' ' << dot1 << ' ' << dot2 << ' ' << dot3 << ' ' 
				<< eigenValues[0] << ' ' << eigenValues[1] << ' ' << eigenValues[2] << std::endl;
		}
		else
		{
			std::cout << time << ' ' <<  eigenValues[0] << ' ' << eigenValues[1] << ' ' << eigenValues[2] << std::endl;
		}
		time+=System.readStoreInterval();
		frame++;
	}
	time_t end;
	std::time(&end);
	std::cerr << "Loads at " << static_cast<double>(frame)/static_cast<double>(end-start) << " frames per second ";
	std::cerr << " in " << end-start << " seconds!\n";
	//close our old files
	xyzFile.close();
	
	return 0;
}

