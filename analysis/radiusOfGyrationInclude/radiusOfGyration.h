
template <typename IT>
threeVector<double> centerOfMassIndex(IT begin, IT end, position<double> *p, threeVector<double> size)
{
	int n=std::distance(begin,end);
	threeVector<double> com=0;
	for(IT zero=begin;begin!=end;++begin)
	{
		threeVector<double> d;
		d.x=p[*begin].x-p[*zero].x;
		d.y=p[*begin].y-p[*zero].y;
		d.z=p[*begin].z-p[*zero].z;
		threeVector<double> offset=0;
		if(d.x>size.x/2.0) offset.x=-size.x;
		if(d.y>size.y/2.0) offset.y=-size.y;
		if(d.z>size.z/2.0) offset.z=-size.z;
		if(d.x<-size.x/2.0) offset.x=size.x;
		if(d.y<-size.y/2.0) offset.y=size.y;
		if(d.z<-size.z/2.0) offset.z=size.z;
		com.x+=p[*begin].x+offset.x;
		com.y+=p[*begin].y+offset.y;
		com.z+=p[*begin].z+offset.z;
	}
	if(n!=0)
	{
		com.x/=n;
		com.y/=n;
		com.z/=n;
	}
	return com;
}

template <typename IT>
double radiusOfGyrationBasisIndex(IT begin, IT end, position<double> *p, threeVector<double> size, threeVector<double> basis)
{
	auto tmp=begin;
	//rog=(1/N)*sum((r_i-rAvg),i)
	threeVector<double> rAvg=0;
	int n=std::distance(begin,end);
	for(IT zero=begin;begin!=end;++begin)
	{
		threeVector<double> d;
		d.x=p[*begin].x-p[*zero].x;
		d.y=p[*begin].y-p[*zero].y;
		d.z=p[*begin].z-p[*zero].z;
		threeVector<double> offset=0;
		if(d.x>size.x/2.0) offset.x=-size.x;
		if(d.y>size.y/2.0) offset.y=-size.y;
		if(d.z>size.z/2.0) offset.z=-size.z;
		if(d.x<-size.x/2.0) offset.x=size.x;
		if(d.y<-size.y/2.0) offset.y=size.y;
		if(d.z<-size.z/2.0) offset.z=size.z;
		rAvg.x+=p[*begin].x+offset.x;
		rAvg.y+=p[*begin].y+offset.y;
		rAvg.z+=p[*begin].z+offset.z;
	}
	if(n!=0)
	{
		rAvg.x/=n;
		rAvg.y/=n;
		rAvg.z/=n;
	}
	begin=tmp;
	double rog=0;
	for(IT zero=begin;begin!=end;++begin)
	{
		threeVector<double> d;
		d.x=p[*begin].x-p[*zero].x;
		d.y=p[*begin].y-p[*zero].y;
		d.z=p[*begin].z-p[*zero].z;
		threeVector<double> offset=0;
		if(d.x>size.x/2.0) offset.x=-size.x;
		if(d.y>size.y/2.0) offset.y=-size.y;
		if(d.z>size.z/2.0) offset.z=-size.z;
		if(d.x<-size.x/2.0) offset.x=size.x;
		if(d.y<-size.y/2.0) offset.y=size.y;
		if(d.z<-size.z/2.0) offset.z=size.z;
		threeVector<double> r;
		r.x=p[*begin].x-rAvg.x+offset.x;
		r.y=p[*begin].y-rAvg.y+offset.y;
		r.z=p[*begin].z-rAvg.z+offset.z;
		rog+=r.x*r.x*basis.x+r.y*r.y*basis.y+r.z*r.z*basis.z;
	}
	if(n!=0)
		rog/=n;
	return rog;
}

//
template <typename IT>
threeVector<double> radiusOfGyrationAxisIndex(IT begin, IT end, position<double> *p, threeVector<double> size)
{
	threeVector<double> rog=0,basis;
	basis.x=0;basis.y=1;basis.z=1;
	rog.x=radiusOfGyrationBasisIndex(begin,end,p,size,basis);
	basis.x=1;basis.y=0;basis.z=1;
	rog.y=radiusOfGyrationBasisIndex(begin,end,p,size,basis);
	basis.x=1;basis.y=1;basis.z=0;
	rog.z=radiusOfGyrationBasisIndex(begin,end,p,size,basis);
	
	return rog;
}

template <typename IT>
double radiusOfGyrationIndex(IT begin, IT end, position<double> *p, threeVector<double> size)
{
	auto tmp=begin;
	//rog=(1/N)*sum((r_i-rAvg),i)
	threeVector<double> rAvg=0;
	int n=std::distance(begin,end);
	for(IT zero=begin;begin!=end;++begin)
	{
		threeVector<double> d;
		d.x=p[*begin].x-p[*zero].x;
		d.y=p[*begin].y-p[*zero].y;
		d.z=p[*begin].z-p[*zero].z;
		threeVector<double> offset=0;
		if(d.x>size.x/2.0) offset.x=-size.x;
		if(d.y>size.y/2.0) offset.y=-size.y;
		if(d.z>size.z/2.0) offset.z=-size.z;
		if(d.x<-size.x/2.0) offset.x=size.x;
		if(d.y<-size.y/2.0) offset.y=size.y;
		if(d.z<-size.z/2.0) offset.z=size.z;
		rAvg.x+=p[*begin].x+offset.x;
		rAvg.y+=p[*begin].y+offset.y;
		rAvg.z+=p[*begin].z+offset.z;
	}
	if(n!=0)
	{
		rAvg.x/=n;
		rAvg.y/=n;
		rAvg.z/=n;
	}
	begin=tmp;
	double rog=0;
	for(IT zero=begin;begin!=end;++begin)
	{
		threeVector<double> d;
		d.x=p[*begin].x-p[*zero].x;
		d.y=p[*begin].y-p[*zero].y;
		d.z=p[*begin].z-p[*zero].z;
		threeVector<double> offset=0;
		if(d.x>size.x/2.0) offset.x=-size.x;
		if(d.y>size.y/2.0) offset.y=-size.y;
		if(d.z>size.z/2.0) offset.z=-size.z;
		if(d.x<-size.x/2.0) offset.x=size.x;
		if(d.y<-size.y/2.0) offset.y=size.y;
		if(d.z<-size.z/2.0) offset.z=size.z;
		threeVector<double> r;
		r.x=p[*begin].x-rAvg.x+offset.x;
		r.y=p[*begin].y-rAvg.y+offset.y;
		r.z=p[*begin].z-rAvg.z+offset.z;
		rog+=r.x*r.x+r.y*r.y+r.z*r.z;
	}
	if(n!=0)
		rog/=n;
	return rog;
}
