#include <iostream>
#include <vector>
#include <cstdio>
#include "dataTypes.h"

template <typename T>
struct threeVectorNoUnion {
	T x,y,z;
};

int main()
{
	std::cout << "alignment of threeVectorNoUnion<double>: " << alignof (threeVectorNoUnion<double>) << std::endl; 
	std::cout << "alignment of threeVector<double>: " << alignof (threeVector<double>) << std::endl; 
	std::cout << "alignment of fourVector<double>: " << alignof (fourVector<double>) << std::endl; 
	std::cout << "alignment of position<double>: " << alignof (position<double>) << std::endl; 
	std::cout << std::endl;
	
	std::cout << "size of threeVectorNoUnion<double>: " << sizeof(threeVectorNoUnion<double>) << std::endl; 
	std::cout << "size of threeVector<double>: " << sizeof(threeVector<double>) << std::endl; 
	std::cout << "size of fourVector<double>: " << sizeof(fourVector<double>) << std::endl; 
	std::cout << "size of position<double>: " << sizeof(position<double>) << std::endl; 
	std::cout << std::endl;
	
	std::cout << "alignment of threeVectorNoUnion<float>: " << alignof (threeVectorNoUnion<float>) << std::endl; 
	std::cout << "alignment of threeVector<float>: " << alignof (threeVector<float>) << std::endl; 
	std::cout << "alignment of fourVector<float>: " << alignof (fourVector<float>) << std::endl; 
	std::cout << "alignment of position<float>: " << alignof (position<float>) << std::endl; 
	std::cout << std::endl;
	
	std::cout << "size of threeVectorNoUnion<float>: " << sizeof(threeVectorNoUnion<float>) << std::endl; 
	std::cout << "size of threeVector<float>: " << sizeof(threeVector<float>) << std::endl; 
	std::cout << "size of fourVector<float>: " << sizeof(fourVector<float>) << std::endl; 
	std::cout << "size of position<float>: " << sizeof(position<float>) << std::endl;
	std::cout << std::endl;
 
	std::cout << "alignment of threeVectorNoUnion<int>: " << alignof(threeVectorNoUnion<int>) << std::endl; 
	std::cout << "alignment of threeVector<int>: " << alignof(threeVector<int>) << std::endl; 
	std::cout << "alignment of fourVector<int>: " << alignof(fourVector<int>) << std::endl; 
	std::cout << "alignment of position<int>: " << alignof(position<int>) << std::endl; 
	std::cout << std::endl;
	
	std::cout << "size of threeVectorNoUnion<int>: " << sizeof(threeVectorNoUnion<int>) << std::endl; 
	std::cout << "size of threeVector<int>: " << sizeof(threeVector<int>) << std::endl; 
	std::cout << "size of fourVector<int>: " << sizeof(fourVector<int>) << std::endl; 
	std::cout << "size of position<int>: " << sizeof(position<int>) << std::endl; 
	std::cout << std::endl;

	std::cout << "alignment of int: " << alignof(int) << std::endl; 
	std::cout << "alignment of double: " << alignof(double) << std::endl; 
	std::cout << "alignment of float: " << alignof(float) << std::endl; 
	std::cout << "alignment of char: " << alignof(char) << std::endl; 
	std::cout << std::endl;
	
	std::cout << "size of int: " << sizeof(int) << std::endl; 
	std::cout << "size of double: " << sizeof(double) << std::endl; 
	std::cout << "size of float: " << sizeof(float) << std::endl; 
	std::cout << "size of char: " << sizeof(char) << std::endl; 

	std::vector< threeVector<double> > a;
	threeVector<double> b;
	
	b.x=90128734;
	b.y=9128347;
	b.z=932458;
	
	a.push_back(b);
	
	std::cout << sizeof(a[0]) << '\t' << sizeof(b) << '\t' << 
		sizeof(unsigned int) << '\t' << sizeof(double) << std::endl;
	std::cout << alignof(a[0]) << '\t' << alignof(b) << '\t' << 
		alignof(unsigned int) << '\t' << alignof(double) << std::endl;
	
	unsigned int *unalignedA=(unsigned int*)(&a[0]);
	for(int i=0;i<alignof(a[0])/alignof(unsigned int);i++)
		printf("%X\t",unalignedA[i]);
	std::cout << std::endl;
	
	unsigned int *unalignedB=(unsigned int*)(&b);
	for(int i=0;i<alignof(b)/alignof(unsigned int);i++)
		printf("%X\t",unalignedB[i]);
	std::cout << std::endl;
	
	return 0;
}
