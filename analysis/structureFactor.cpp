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

bool FFT2D(twoVector<double> *c,int nx,int ny,int dir);
bool FFT(int dir,int m,double *x,double *y);
bool Powerof2(int n,int *m,int *twopm);

int main(int argc, char **argv)
{
	if(argc!=6)
	{
		std::cerr << argv[0] << " filename sX sY nX nY\n";
		std::cerr << "\tTakes a height map, does a 2D fft, then circularly averages it.\n";
		return 0;
	}
	
	std::stringstream cmdArg;
	for(int i=2;i<argc;i++)
		cmdArg << argv[i] << '\t';
	
	threeVector<int> n;
	threeVector<double> s;
	cmdArg >> s.x >> s.y >> n.x >> n.y;
	
	std::fstream inFile;
	inFile.open(argv[1], std::ios::in);

	std::vector<twoVector<double> > data;
	
	while(inFile.is_open() && !inFile.eof())
	{
		twoVector<double> input;
		input.y=0;
		inFile >> input.x;
		if(!inFile.eof())
			data.push_back(input);
	}
	
	double param=(s.x*s.y)/(n.x*n.y);
	std::vector<double> structureFactor(n.x,0.0);
	int nFrames=data.size()/(n.x*n.y);
	std::vector<double> fc;
	for(int i=0;i<n.x;i++)
		fc.push_back(static_cast<double>(i)*M_PI*2.0/s.x);
	
	std::vector<double> count(n.x,0);
	
	for(int i=0;i<data.size();i+=n.x*n.y)
	{
		FFT2D(&data[i], n.x, n.y, 1);
		
		//std::vector<double> structureFactor(n.x,0.0);
		
		for(int j=0;j<n.x*n.y;j++)
		{
			twoVector<int> index;
			index.x=j%n.x;
			index.y=j/n.y;
			index.x=(index.x>=n.x/2)?n.x-index.x-1:index.x;
			index.y=(index.y>=n.y/2)?n.y-index.y-1:index.y;
			int qq=sqrt(index.x*index.x+index.y*index.y);
			if(qq>structureFactor.size())
			{
				std::cerr << "qq is too large! " << qq << '\t' << index.x << '\t' << index.y << '\t' << j << std::endl;
				throw 0;
			}
			structureFactor[qq]+=(data[j+i].x*data[j+i].x+data[j+i].y*data[j+i].y)*param;
			count[qq]++;
		}
		//for(int j=0;j<n.x/2;j++)
		//	std::cout << fc[j] << '\t' << structureFactor[j] << std::endl;
	}
	
	for(int i=0;i<n.x/2;i++)
		std::cout << fc[i] << '\t' << structureFactor[i]/static_cast<double>(count[i]) << std::endl;
	
	return 0;
}

/*
void stfac(double xlx, double yly, int lmn, position<double> *p, int nP)
{
	int nn=2;
	int l=1024;
	int nsite=l*l;
	int number_q=l/2;
	double data[2*nsite];
	double sf[nsite];
	double cr[nsite];
	double sfunc[2*number_q];
	double num[2*number_q];
	double corel[2*number_q];
	int ndim[2];
	
	double xsize=xlx/double(l);
	double ysize=yly/double(l);
	ndim[0]=l;
	ndim[1]=l;
	double param=xlx*yly/double(nsite)
	for(int i=0;i<l;i++)
	{
		sfunc[i]=0;
		corel[i]=0;
	}
	for(int i=0;i<2*nsite;i++)
	{
		data(i)=0;
	}
	for(int i=0;i<lmn;i++)
	{
		int 
		int ix=(p[i].x/xlx)*l+1;
		int iy=(p[i].y/yly)*l+1;
		int j=ix+(iy-1)*l;
		int ind1=2*j;
		int ind2=ind1-1;
		data[ind1]=0.0;//the imaginary numbers for fourn
		data[ind2]=data[ind2]+1;//my heightmap
	}
	
	//use the transform space
	fourn(data,ndim,2,1);
	
	//now calculate structure factor from the transformed space
	do i=1,nsite
	  ind1=2*i
	  ind2=ind1-1
	  sf(i)=(data(ind1)**2+data(ind2)**2)*param
	enddo
	do i=1,nsite
	  ix=mod(i-1,l)
	  iy=int((i-1)/l)
	  if(ix.gt.l/2)ix=l-ix
	  if(iy.gt.l/2)iy=l-iy
	  qq=dsqrt(dfloat(ix**2+iy**2))
	  kval=qq
	  sfunc(kval)=sfunc(kval)+sf(i)
	  num(kval)=num(kval)+1
	enddo
	open(10,file='stfac.data',access='append')
	do i=1,l/2
	  if(num(i).gt.0)then
	    sfunc(i)=sfunc(i)/num(i)
	  endif
	  qq=pi*i/xlx
	  write(10,*)qq,sfunc(i)
	enddo
	close(10)
	
	//put it back in real space, averaging sf data
	do i=1,nsite
	  ind1=2*i
	  ind2=ind1-1
	  data(ind1)=0.0
	  data(ind2)=sf(i)/nsite
	enddo
	call fourn(data,ndim,2,-1)
	do i=1,nsite
	  ind1=2*i
	  ind2=ind1-1
	  cr(i)=data(ind2)
	enddo
	do i=1,nsite
	  ix=mod(i-1,l)
	  iy=int((i-1)/l)
	  if(ix.gt.l/2)ix=l-ix
	  if(iy.gt.l/2)iy=l-iy
	  rr=dsqrt(dfloat(ix**2+iy**2))
	  kval=rr
	  num(kval)=num(kval)+1
	  corel(kval)=corel(kval)+cr(i)
	enddo
	open(10,file='corel.data',access='append')
	deltax=xlx/l
	do i=1,l/2
	  if(num(i).gt.0)then
	    corel(i)=corel(i)/num(i)
	  endif
	  rr=i*deltax
	  write(10,*)rr,corel(i)
	enddo
	close(10)
}*/

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
