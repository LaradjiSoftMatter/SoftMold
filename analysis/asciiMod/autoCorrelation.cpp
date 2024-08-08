//program for averaging column data
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
using namespace std;


#define PARTIAL_SUM_LIMIT 5
double pairSum(double *list, int start, int end)
{
	if(end-start<=PARTIAL_SUM_LIMIT)
	{
		double sum=0;
		for(int i=start;i<end;i++)
			sum+=list[i];
		return sum;
	}
	else
	{
		return pairSum(list, start, start+int((end-start)/2))+pairSum(list, start+int((end-start)/2), end);
	}
}

typedef struct {
	double x,y,std,ms;
} element;

int main(int argc, char **argv)
{
	int step=-1,nCol=2,col=1,nStep=-1,offset=0, discard=0;
	double xStart=0,xEnd=1,xStep=-1;
	bool xStepFlag=false;
	int *fileArgN=new int[argc];
	int nFiles=0;
	bool std=false;
	bool ms=false;
	
	//options
	for(int i=1;i<argc;i++)
	{
		if(argv[i][0]=='-')
		{
			bool found=false;
			if((strcmp(argv[i],"-h")==0) || (strcmp(argv[i],"-help")==0))
			{
				found=true;
				cout << "Usage: " << argv[0] << " File Options\n";
				cout << "Options:\n";
				cout << "-nCol [integer]\n\tNumber of columns present in data.\n";
				cout << "-col [integer]\n\tColumn to average.\n";
			}
			if(strcmp(argv[i],"-nCol")==0)
			{
				found=true;
				i++;
				nCol=atoi(argv[i]);
			}
			if(strcmp(argv[i],"-col")==0)
			{
				found=true;
				i++;
				col=atoi(argv[i]);
			}
			if(!found)
			{
				cout << "Unrecognized option \"" << argv[i] << "\".\n\n";
				
				cout << "Usage: " << argv[0] << " File1 File2... Options\n";
				cout << "Options:\n";
				cout << "-nCol [integer]\n\tNumber of columns present in data.\n";
				cout << "-col [integer]\n\tColumn to average.\n";
				
				return 0;
			}
		}
		else
		{
			fileArgN[nFiles++]=i;
		}
	}
	
	if(nFiles==0)
	{
		cout << "No input files!\n-help for list of options and format.\n";
		return 0;
	}
	
	
	
	
	
	for(int file=0;file<nFiles;file++)
	{
		//derived variables
		int nData;
		fstream input;
		double *line=new double[nCol];
		double *values, *stdVal, *msVal;
		element *data;
	
		//initialize derived variables
		input.open(argv[fileArgN[file]],ios::in);
		
		if(!input.is_open())
		{
			cout << argv[fileArgN[file]] << " not found!\n";
			delete fileArgN;
			return 0;
		}

		
		for(nData=0;!input.eof();nData++)//count number of data points in file
		{
			for(int i=0;i<nCol;i++)
			{
				if(input.eof())
				{
					cout << "End of file encountered while reading!\nPossibly due to wrong number of columns!\n";
					delete line;
					delete fileArgN;
					return 0;
				}
				input >> line[i];
				if(input.eof())
					break;
			}
		}
	
		input.clear();
		input.seekg(ios::beg);//back to beginning of file
	
		if(double((nData-offset)/step)<double(nStep) && nStep!=-1)
		{
			cout << "Not enough data points (" << nData-offset << ") to finish all steps (" << nStep << ") of size " << step << "!\n";
			delete line;
			delete fileArgN;
			return 0;
		}
	
		if(step==-1)
			step=nData-offset-1;
		
		if(nStep==-1)
			nStep=(nData-offset-1)/step;//this becomes one if step initially was -1
		
		data=new element[nStep];
		values=new double[step];
		if(std)
			stdVal=new double[step];
		if(ms)
			msVal=new double[step];
		
		if(!xStepFlag)
			xStep=double(xEnd-xStart)/double(nStep);

	
		
		//move file to offset
		for(int j=0;j<offset;j++)
			for(int i=0;i<nCol;i++)
				input >> line[i];
		
		
		//compute running average
		for(int j=0;j<nStep;j++)
		{
			data[j].x=xStart+(xStep*double(j));
			//average of step
			for(int k=0;k<step;k++)
			{
				for(int i=0;i<nCol;i++)
				{
					input >> line[i];
				}
				values[k]=line[col];
			}
			data[j].y=pairSum(values,discard,step)/double(step-discard);
			//standard deviation of step, if active
			if(std)
			{
				double avg=data[j].y;
				double buf;
				for(int k=discard;k<step;k++)
				{
					buf=values[k]-avg;
					stdVal[k]=buf*buf;
				}
				data[j].std=sqrt(pairSum(stdVal,discard,step)/double(step-discard-1.0));
			}
			if(ms)
			{
				double avg=data[j].y;
				double buf;
				for(int k=discard;k<step;k++)
				{
					buf=values[k]-avg;
					msVal[k]=buf*buf;
				}
				data[j].ms=pairSum(msVal,discard,step)/double(step-discard-1.0);
			}
			cout << data[j].x << '\t' << data[j].y;
			if(std)
				cout << '\t' << data[j].std;
			if(ms)
				cout << '\t' << data[j].ms;
			cout << endl;
		}
	
		delete values;
		delete line;
		delete data;
		if(std)
			delete stdVal;
		if(ms)
			delete msVal;
	}

	delete fileArgN;
	return 0;
}
