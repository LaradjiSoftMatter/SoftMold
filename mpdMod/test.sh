#!/bin/bash

g++ -O3 insertModel.cpp -o ~/bin/insertModel
ls #initial directory

bunzip2 in.mpd.bz2
ls #should see in.mpd

g++ circleVert.cpp
./a.out
ls #should see a.out, circleVert.dat, and circleEdge.dat

#Usage: insertModel oldName newName vertices.dat edges.dat kBond lBond pos.x pos.y pos.z
insertModel in out circleVert.dat circleEdge.dat 100 1.25 30 0 30
ls #should see an out.mpd now

#use these two commands to output configuration_name.xyz files
#compile extractXYZ first (this is a useful utility):
#g++ -O3 ../analysis/extractXYZ.cpp -o ~/bin/extractXYZ
#extractXYZ in
#extractXYZ out

MD out
