#!/usr/bin/awk

FNR==NR{
a[$1]=$2 FS $3;
next;
}{ 
print $0, a[$1];
}
