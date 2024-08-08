//program for averaging column data
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
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
	double x,y,std,ms,binLen;
} element;

int main(int argc, char **argv)
{
	int step=-1,nCol=2,col=1,nStep=-1,offset=0, discard=0;
	double xStart=0,xEnd=1,xStep=-1;
	bool xStepFlag=false;
	//int *fileArgN=new int[argc];
	std::vector<int> fileArgN;
	bool std=false;
	bool ms=false;
	char *label=NULL;
	//options
	for(int i=1;i<argc;i++)
	{
		if(argv[i][0]=='-' && argv[i][1]=='-')
		{
			bool found=false;
			if((strcmp(argv[i],"--h")==0) || (strcmp(argv[i],"--help")==0))
			{
				found=true;
				cerr << "Usage: " << argv[0] << " File Options\n";
				cerr << "Options:\n";
				cerr << "--step [integer]\n\tNumber of values to be averaged per step.\n";
				cerr << "--nCol [integer]\n\tNumber of columns present in data.\n";
				cerr << "--col [integer]\n\tColumn to average.\n";
				cerr << "--nStep [integer]\n\tNumber of desired steps through data.\n";
				cerr << "--offset [integer]\n\tStart at this offset in data.\n";
				cerr << "--discard [integer]\n\tDiscard this number of data points per step.\n";
				cerr << "--xStart [double]\n\tStart the x column at this value.\n";
				cerr << "--label [string]\n\tLabel the x column at this value. No spaces or quotes!\n";
				cerr << "--xEnd [double]\n\tEnd the x column at this value.\n";
				cerr << "--xStep [double]\n\tHow much to step x column.\n";
				cerr << "--std \n\tOutput standard deviation of data in 3rd column.\n";
				cerr << "--ms \n\tOutput mean square of data in 3rd or 4th column.\n";
			}
			if(strcmp(argv[i],"--step")==0)
			{
				found=true;
				i++;
				step=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--nCol")==0)
			{
				found=true;
				i++;
				nCol=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--col")==0)
			{
				found=true;
				i++;
				col=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--nStep")==0)
			{
				found=true;
				i++;
				nStep=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--offset")==0)
			{
				found=true;
				i++;
				offset=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--discard")==0)
			{
				found=true;
				i++;
				discard=atoi(argv[i]);
			}
			if(strcmp(argv[i],"--xStart")==0)
			{
				found=true;
				i++;
				xStart=atof(argv[i]);
			}
			if(strcmp(argv[i],"--label")==0)
			{
				found=true;
				i++;
				label=argv[i];
			}
			if(strcmp(argv[i],"--xEnd")==0)
			{
				found=true;
				i++;
				xEnd=atof(argv[i]);
			}
			if(strcmp(argv[i],"--xStep")==0)
			{
				xStepFlag=true;
				found=true;
				i++;
				xStep=atof(argv[i]);
			}
			if(strcmp(argv[i],"--std")==0)
			{
				found=true;
				std=true;
			}
			if(strcmp(argv[i],"--ms")==0)
			{
				found=true;
				ms=true;
			}
			
			if(!found)
			{
				cerr << "Unrecognized option \"" << argv[i] << "\".\n\n";
				
				cerr << "Usage: " << argv[0] << " File Options\n";
				cerr << "Options:\n";
				cerr << "--step [integer]\n\tNumber of values to be averaged per step.\n";
				cerr << "--nCol [integer]\n\tNumber of columns present in data.\n";
				cerr << "--col [integer]\n\tColumn to average.\n";
				cerr << "--nStep [integer]\n\tNumber of desired steps through data.\n";
				cerr << "--offset [integer]\n\tStart at this offset in data.\n";
				cerr << "--discard [integer]\n\tDiscard this number of data points per step.\n";
				cerr << "--xStart [double]\n\tStart the x column at this value.\n";
				cerr << "--label [string]\n\tLabel the x column at this value. No spaces or quotes!\n";
				cerr << "--xEnd [double]\n\tEnd the x column at this value.\n";
				cerr << "--xStep [double]\n\tHow much to step x column.\n";
				cerr << "--std \n\tOutput standard deviation of data in 3rd column.\n";
				cerr << "--ms \n\tOutput mean square of data in 3rd or 4th column.\n";
				
				return 0;
			}
		}
		else
		{
			fileArgN.push_back(i);
		}
	}
	
	if(fileArgN.size()==0)
	{
		cerr << "No input files!\n--help for list of options and format.\n";
		return 0;
	}
	
	
	
	
	
	for(int &file:fileArgN)
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
			cerr << argv[fileArgN[file]] << " not found!\n";
			return 0;
		}
		
		
		for(nData=0;!input.eof();nData++)//count number of data points in file
		{
			for(int i=0;i<nCol;i++)
			{
				if(input.eof())
				{
					cerr << "End of file encountered while reading!\nPossibly due to wrong number of columns!\n";
					delete line;
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
			cerr << "Not enough data points (" << nData-offset << ") to finish all steps (" << nStep << ") of size " << step << "!\n";
			delete line;
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
			if(label!=NULL)
				cout << label << '\t';
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

	return 0;
}
