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

#include <fftw3.h>

bool FFT2D(twoVector<double> *c,int nx,int ny,int dir);
bool FFT(int dir,int m,double *x,double *y);
bool Powerof2(int n,int *m,int *twopm);

int main(int argc, char **argv)
{
	if(argc<5)
	{
		std::cerr << argv[0] << " name fromTime nX nY type1 type2 ...\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	double fromTime=0;
	threeVector<int> n;
	cmdArg >> fromTime >> n.x >> n.y;

	std::vector<int> types;
	
	for(int i=4;i<argc;i++)
	{
		int type;
		cmdArg >> type;
		types.push_back(type);
	}
	
	//the variables for the simulation, call with empty constructor
	Blob<double> System;
	
	//file input-output script engine
	//name part (argv[1] in this case) is used with no extension (*.mpd)
	Script<double, Blob <double> > fileIO(argv[1],std::ios::in,&System);
	
	//read the variables in
	fileIO.read();
	
	//close the file that it is reading from
	fileIO.close();
	
	//set time
	double time=0;

	//grab sizes
	threeVector<double> size=System.readSize();
	
	std::string newName("size_");
	newName+=argv[1];
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
		time=sTime;
		size=sCurrent;
		//std::cout << size.x << '\t' << size.y << '\t' << size.z << '\n';
	}
		
	//Get our frames file
	std::string framesName("frames_");
	framesName+=argv[1];
	framesName+=".xyz";
	
	position<double> *p=System.getPositions();
	int nParticles=System.readNParticles();
	
	xyzFormat<double> xyzFile(p,nParticles);
	xyzFile.open(framesName.c_str(), std::ios::in);
	
	std::vector<double> structureFactor(n.x,0.0);
	int frames=0;
	
	twoVector<double> zero;
	zero.x=0;
	zero.y=0;
	std::vector< twoVector<double> > heightMap(n.x*n.y, zero);
	fftw_plan forw3;
	forw3=fftw_plan_dft_2d(n.x, n.y,reinterpret_cast<fftw_complex*>(&heightMap[0]),
		reinterpret_cast<fftw_complex*>(&heightMap[0]), FFTW_FORWARD, FFTW_ESTIMATE);
	
	double avgSize=0;
	while(xyzFile.load())
	{
		std::vector< std::vector< position<double> > > cellList;
		
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
		if(time>=fromTime)
		{
			//make the height map
			threeVector<double> r;
			r.x=size.x/static_cast<double>(n.x);
			r.y=size.y/static_cast<double>(n.y);
			
			
			std::vector<int> countMap(n.x*n.y, 0.0);
			for(int i=0;i<heightMap.size();i++)
				heightMap[i]=zero;
			
			
			double avgZ=0;
			double countZ=0;
			
			for(int i=0;i<System.readNParticles();i++)
			{
				bool match=false;
				for(int tIndex=0;!match && tIndex<types.size();tIndex++)
					match=(p[i].type==types[tIndex]);
				
				if(match)
				{
					threeVector<int> c;
					c.x=floor(p[i].x/r.x);
					c.y=floor(p[i].y/r.y);
					c.x=(c.x>=n.x)?n.x-1:c.x;
					c.y=(c.y>=n.y)?n.y-1:c.y;
					int hash=c.x+c.y*n.x;
					if(hash>=heightMap.size())
					{
						std::cerr << "Out of bounds particle!\n";
						std::cerr << "Size: " << size.x << '\t' << size.y << std::endl;
						std::cerr << "Position: " << p[i].x << '\t' << p[i].y << std::endl;
					}
					heightMap[hash].x+=p[i].z;
					countMap[hash]++;
					avgZ+=p[i].z;
					countZ++;
				}
			}
			
			if(countZ!=0)
				avgZ/=countZ;
			//get the average heightmap
			for(int i=0;i<heightMap.size();i++)
			{
				if(countMap[i]!=0)
				{
					heightMap[i].x/=static_cast<double>(countMap[i]);
					heightMap[i].x-=avgZ;
				}
			}
				
			//transform it
			//FFT2D(&heightMap[0], n.x, n.y, 1);
			fftw_execute(forw3);
			
			std::vector<double> count(n.x,0);
			
			double param=(r.x*r.y)/(size.x*size.y);
			std::vector<double> sF(n.x,0.0);
			
			for(int j=0;j<n.x*n.y;j++)
			{
				twoVector<int> index;
				index.x=j%n.x;
				index.y=j/n.y;
				index.x=(index.x>=n.x/2)?n.x-index.x:index.x;
				index.y=(index.y>=n.y/2)?n.y-index.y:index.y;
				int qq=sqrt(index.x*index.x+index.y*index.y);
				if(qq>sF.size())
				{
					std::cerr << "qq is too large! " << qq << '\t' << index.x 
						<< '\t' << index.y << '\t' << j << std::endl;
					throw 0;
				}
				sF[qq]+=(heightMap[j].x*heightMap[j].x+heightMap[j].y*heightMap[j].y)*param;
				count[qq]++;
			}
			
			std::vector<double> fc;
			for(int i=0;i<n.x;i++)
				fc.push_back(static_cast<double>(i)*2.0*M_PI/(size.x));
			for(int i=1;i<n.x/2;i++)
				if(count[i]!=0)
					std::cout << fc[i] << '\t' << sF[i]/static_cast<double>(count[i]) << std::endl;
				else
					std::cout << fc[i] << '\t' << 0 << std::endl;
			frames++;
			avgSize+=size.x;
			std::cerr << time << std::endl;
		}
		time+=System.readStoreInterval();
	}
	avgSize/=frames;
	fftw_destroy_plan(forw3);
	
	return 0;
}

/*-------------------------------------------------------------------------
   Perform a 2D FFT inplace given a complex 2D array
   The direction dir, 1 for forward, -1 for reverse
   The size of the array (nx,ny)
   Return false if there are memory problems or
      the dimensions are not powers of 2
*/
bool FFT2D(twoVector<double> *c,int nx,int ny,int dir)
{
   int i,j;
   int m,twopm;
   double *x,*y;

   /* Transform the rows */
   x = (double *)malloc(nx * sizeof(double));
   y = (double *)malloc(nx * sizeof(double));
   if (x == NULL || y == NULL)
      return false;
   if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return false;
   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         x[i] = c[i+nx*j].x;
         y[i] = c[i+nx*j].y;
      }
      FFT(dir,m,x,y);
      for (i=0;i<nx;i++) {
         c[i+nx*j].x = x[i];
         c[i+nx*j].y = y[i];
      }
   }
   free(x);
   free(y);

   /* Transform the columns */
   x = (double *)malloc(ny * sizeof(double));
   y = (double *)malloc(ny * sizeof(double));
   if (x == NULL || y == NULL)
      return false;
   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return false;
   for (i=0;i<nx;i++) {
      for (j=0;j<ny;j++) {
         x[j] = c[i+nx*j].x;
         y[j] = c[i+nx*j].y;
      }
      FFT(dir,m,x,y);
      for (j=0;j<ny;j++) {
         c[i+nx*j].x = x[j];
         c[i+nx*j].y = y[j];
      }
   }
   free(x);
   free(y);

   return true;
}

/*-------------------------------------------------------------------------
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform

     Formula: forward
                  N-1
                  ---
              1   \          - j k 2 pi n / N
      X(n) = ---   >   x(k) e                    = forward transform
              N   /                                n=0..N-1
                  ---
                  k=0

      Formula: reverse
                  N-1
                  ---
                  \          j k 2 pi n / N
      X(n) =       >   x(k) e                    = forward transform
                  /                                n=0..N-1
                  ---
                  k=0
*/
bool FFT(int dir,int m,double *x,double *y)
{
   long nn,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   nn = 1;
   for (i=0;i<m;i++)
      nn *= 2;

   /* Do the bit reversal */
   i2 = nn >> 1;
   j = 0;
   for (i=0;i<nn-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<nn;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<nn;i++) {
         x[i] /= (double)nn;
         y[i] /= (double)nn;
      }
   }

   return true;
}

/*-------------------------------------------------------------------------
   Calculate the closest but lower power of two of a number
   twopm = 2**m <= n
   Return TRUE if 2**m == n
*/
bool Powerof2(int n,int *m,int *twopm)
{
   if (n <= 1) {
      *m = 0;
      *twopm = 1;
      return false;
   }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
      return false;
   else
      return true;
}
