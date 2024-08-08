Look at test.sh. It shows how to use the insertModel code to insert an arbitrary set of vertices and edges into a simulation.

The general usage for insertModel is:
	insertModel oldName newName vertices.dat edges.dat kBond lBond pos.x pos.y pos.z
Where oldName is an oldName.mpd file for input, and newName is a newName.mpd file for output.
vertices.dat is the positions in the following column format:
	type x y z
where type is a positive integer, and x, y, and z are the coordinates for the vertices.
edges.dat are the bonds in the following column format:
	indexA indexB
where indexA is connected to indexB using the line indices from the vertices.dat file.
kBond and lBond are constants for the harmonic potential:
	U(distance)=(kBond/2.0)*(distance-lBond)^2
The three coordinates, pos.x, pos.y, and pos.z, will place the vertices' center of mass at that coordinate. If any vertex goes over the system edge, it will be placed on the other half of the system using the periodic image convention.
