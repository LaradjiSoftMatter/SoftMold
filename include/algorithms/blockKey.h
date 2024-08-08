#ifndef BLOCKKEY
#define BLOCKKEY

template <typename T, typename I=int>//, int length=3>
class blockKey {
	public:
		blockKey()
		{
			iteration=0;
			length=3;
			nB=blockEncode(0);
		};
		
		
		blockKey(I leng)
		{
			iteration=0;
			length=leng;
			nB=blockEncode(0);
		};
		
		
		blockKey(I iter, I leng)
		{
			iteration=iter;
			length=leng;
			nB=blockEncode(length);
		};
		
		bool operator!= (const blockKey& other) const
		{
			return iteration!=other.iteration;
		};
		
		bool operator== (const blockKey& other) const
		{
			return iteration==other.iteration;
		};
		
		const T& operator* () const
		{
			return nB;
		};
		
		const blockKey& operator++()
		{
			nB=blockEncode(++iteration);
			return *this;
		};
		
		const blockKey& operator++(int)
		{
			nB=blockEncode(++iteration);
			return *this;
		};
		
		const blockKey& operator--()
		{
			nB=blockEncode(--iteration);
			return *this;
		};
		
		const blockKey& operator--(int)
		{
			nB=blockEncode(--iteration);
			return *this;
		};
			
		/** \brief Indexed operator.
		*
		*/
		inline T& operator[] (int index)
		{
			return blockEncode(index);
		};
		
		blockKey begin() const
		{
			return blockKey(0);
		};
		
		blockKey end() const
		{
			return blockKey(length*length*length, length);
		};
		I iteration;
		
	private:
		T nB;
		I length;
		inline T blockEncode(int index)
		{
			T nnB;
			index%=length*length*length;
			nnB.x=(index)%(length)-static_cast<I>(length/2);
			nnB.y=static_cast<I>(index/(length))%(length)-static_cast<I>(length/2);
			nnB.z=static_cast<I>(index/(length*length))-static_cast<I>(length/2);
			return nnB;
		};
};

#endif
