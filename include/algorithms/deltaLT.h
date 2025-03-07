#ifndef MPD_DELTALT
#define MPD_DELTALT

#include <cmath>
#include <iostream>
#include "dataTypes.h"

namespace mpd {

template <typename T>
struct deltaLT {
	//The system timestep and this length step, deltaT and deltaL,
	// together decide how fast this converges to the end length, endL.
	// RelaxStep decides how often it runs. The difference between the 
	// initial size.s[dim] and endL decide how many steps it takes.
	//E.g.
	// size = {25 55 100}
	// deltaT = 0.02
	// deltaL = 0.01
	// endL = 50 
	// dim = 2 (maybe Z or z in mpd file)
	// relaxStep = 10
	//Would result in 10*(100-50)/0.01=50000 steps taken to get
	// from 100 to 50, or 50000*0.02=1000 tau.
	//This would likely increase the system temperature rapidly
	// if the standard for tau is about 19 ns. Volume would be 1/2
	// within 19 us. A more reasonable deltaL would likely be 0.0001, 
	// shrinking the volume within 1.9 ms (100000 tau).
	//Once the Z dimension reaches that size, the system will 
	// no longer change size.
	threeVector<T> deltaL,endL;
	int relaxStep;
	
	//Just ensure they are all zero, we will use deltaL==0 to indicate inactive
	deltaLT():deltaL(0),endL(0),relaxStep(0){}
	
	//Return a new size based on the activity
	constexpr threeVector<T> newSize(threeVector<T> oldSize)
	{
		threeVector<T> diff=oldSize-endL;
		//Check if it is within the last step
		for(int dim=2;dim!=-1;--dim)
		if(std::abs(diff.s[dim])>deltaL.s[dim])
		{
			T direction=diff.s[dim]/std::abs(diff.s[dim]);
			oldSize.s[dim]-=deltaL.s[dim]*direction;
		}
		else
			//this will make a smaller step to the end length
			oldSize.s[dim]=endL.s[dim];
		
		return oldSize;
	}
	
	//this exists somewhere else, or at least it should. 
	// I should probably make a canned version of this
	constexpr threeVector<T> scaleFactor(threeVector<T> oldSize)
	{
		auto nextSize=newSize(oldSize);
		
		threeVector<T> scale=1.0;
		
		for(int dim=2;dim!=-1;--dim)
			scale.s[dim]=1.0+(nextSize.s[dim]-oldSize.s[dim])/oldSize.s[dim];
		
		return scale;
	}
	
	constexpr bool ready(int step)
	{
		return ((deltaL.x!=0 || deltaL.y!=0 || deltaL.z!=0) && step%relaxStep==0);
	}
	
	constexpr bool active()
	{
		return deltaL.x!=0 || deltaL.y!=0 || deltaL.z!=0;
	}
	
	constexpr int nWords() const {return 7;}
	
	std::istream &inStep(std::istream &stream, int wStep)
	{
		switch(wStep)
		{
			case 0: stream >> deltaL.x; break;
			case 1: stream >> deltaL.y; break;
			case 2: stream >> deltaL.z; break;

			case 3: stream >> endL.x; break;
			case 4: stream >> endL.y; break;
			case 5: stream >> endL.z; break;
			case 6: stream >> relaxStep; break;
			default:
				std::cerr << "Forgot to reset input word counter!";
				throw 1;
				break;
		}
		return stream;
	}
	
	std::ostream &outStep(std::ostream &stream, int wStep)
	{
		switch(wStep)
		{
			case 0: stream << deltaL.x << ' '; break;
			case 1: stream << deltaL.y << ' '; break;
			case 2: stream << deltaL.z << ' '; break;
			case 3: stream << endL.x << ' '; break;
			case 4: stream << endL.y << ' '; break;
			case 5: stream << endL.z << ' '; break;
			case 6: stream << relaxStep << '\n'; break;
			default:
				std::cerr << "Forgot to reset output word counter!";
				throw 1;
				break;
		}
		return stream;
	}
};

}
#endif
