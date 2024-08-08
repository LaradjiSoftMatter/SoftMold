#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "../include/algorithms/dataTypes.h"
/*
template <typename T>
void DFT(T *data, int length)
{
	T *output=new T[length];
	for(int i=0;i<length-1;i++)
	{
		output[i].x=0;
		output[i].y=0;
		for(int j=0;j<length;j++)
		{
			output[i].x+=data[j].x*cos(2.0*M_PI*(double)(i*j)/(double)length);
			output[i].y-=data[j].x*sin(2.0*M_PI*(double)(i*j)/(double)length);
		}
	}
	for(int i=0;i<length;i++)
	{
		data[i].x=abs(output[i].x);
		data[i].y=abs(output[i].y);
	}
	delete output;
}
*/

#define swap(a,b) tempr=(a);(a)=(b);(b)=tempr
/*//NR
void four1(double data[], unsigned long nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    // reverse-binary reindexing
    n = nn<<1;
    j=1;
    for (i=1; i<n; i+=2) {
        if (j>i) {
            swap(data[j], data[i]);
            swap(data[j+1], data[i+1]);
        }
        m = n >> 1;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };

    // here begins the Danielson-Lanczos section
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = isign*(6.28318530717959/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
                tempr = wr*data[j] - wi*data[j+1];
                tempi = wr * data[j+1] + wi*data[j];

                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
	    wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}
*/

//DD
void four1(double* data, unsigned long nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    // reverse-binary reindexing
    n = nn<<1;
    j=1;
    for (i=1; i<n; i+=2) {
        if (j>i) {
            swap(data[j-1], data[i-1]);
            swap(data[j], data[i]);
        }
        m = nn;
        while (m>=2 && j>m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    };

    // here begins the Danielson-Lanczos section
    mmax=2;
    while (n>mmax) {
        istep = mmax<<1;
        theta = -isign*(2*M_PI/mmax);
        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m=1; m < mmax; m += 2) {
            for (i=m; i <= n; i += istep) {
                j=i+mmax;
                tempr = wr*data[j-1] - wi*data[j];
                tempi = wr * data[j] + wi*data[j-1];

                data[j-1] = data[i-1] - tempr;
                data[j] = data[i] - tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wtemp=wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
        }
        mmax=istep;
    }
}

void ZDFT(twoVector<double> *data, int length, int isign)
{
	twoVector<double> *output=new twoVector<double>[length];
	for(int i=0;i<length;i++)
	{
		output[i].x=0;
		output[i].y=0;
		for(int j=0;j<length;j++)
		{
			output[i].x+=data[j].x*cos(isign*2.0*M_PI*(double)(i*j)/(double)length)+data[j].y*sin(isign*2.0*M_PI*(double)(i*j)/(double)length);
			output[i].y+=data[j].y*cos(isign*2.0*M_PI*(double)(i*j)/(double)length)-data[j].x*sin(isign*2.0*M_PI*(double)(i*j)/(double)length);
		}
	}
	if(isign<0)
	{
		for(int i=0;i<length;i++)
			data[i]=output[i]/double(length);
	}
	else
	{
		for(int i=0;i<length;i++)
			data[i]=output[i];
	}
		
	delete output;
}

void ImDFT(double *data, int length)
{
	double *output=new double[length];
	for(int i=0;i<length-1;i++)
	{
		output[i]=0;
		for(int j=0;j<length;j++)
			output[i]-=data[j]*sin(2.0*M_PI*(double)(i*j)/(double)length);
	}
	for(int i=0;i<length;i++)
		data[i]=abs(output[i]);
	delete output;
}

using namespace std;

#define N 8192

int main()
{
	double *input=new double[N*2];
	double *output=new double[N*2];
	
	twoVector<double> *inputDFT=new twoVector<double>[N];
	twoVector<double> *outputDFT=new twoVector<double>[N];
	
	srand(time(NULL));
	
	for(int i=0;i<N;i++)
	{
		input[2*i]=sin((double)i);//2.0*(double)rand()/(double)RAND_MAX-1.0;
		output[2*i]=input[2*i];
		input[2*i+1]=0;
		output[2*i+1]=input[2*i+1];
		inputDFT[i].x=input[2*i];
		inputDFT[i].y=input[2*i+1];
		outputDFT[i]=inputDFT[i];
	}
	
	time_t start=time(NULL);
	{
		four1(output,N,-1);
		four1(output,N,1);
		for(int i=0;i<2*N;i++)
			output[i]/=N;
	}
	time_t end=time(NULL);
	cout << end-start << endl;
	
	start=time(NULL);
	{
		ZDFT(outputDFT,N,-1);
		ZDFT(outputDFT, N,1);
	}
	end=time(NULL);
	cout << end-start << endl;
	
	//DanielsonLanczos<N> fft;
	//fft.apply(output);
	
	//for(int i=0;i<N;i++)
	//{
	//		cout << (double)i << ' ' << input[2*i] << ' ' << output[2*i] << ' ' << output[2*i+1] << ' ' << outputDFT[i].x << ' ' << outputDFT[i].y <<'\n';
	//}

	delete input, output,inputDFT,outputDFT;
	return 0;
}
