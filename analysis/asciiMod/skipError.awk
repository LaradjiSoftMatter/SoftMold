#!/usr/bin/awk
{
	if((NR-1) % sError==0)
	{
		print $1,$2,$3;
	}
	else
	{
		print $1,$2,0;
	}
}
