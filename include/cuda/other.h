#ifndef MPD_OTHER
#define MPD_OTHER

namespace mpd {
	// compute grid and thread block size for a given number of elements
	void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
	{
		numThreads = min(blockSize, n);
		//numBlocks = iDivUp(n, numThreads);
		numBlocks=(n%numThreads!=0)?(n/numThreads+1):(n/numThreads);
	}
}


#endif
