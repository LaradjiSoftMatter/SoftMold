

# Overview:
In the top  (./), ./generate, ./analysis, ./analysis/asciiMod, ./analysis/xyzMod, and ./analysis/functions directories, you can make the executables in the directory with:
```
	make
```
Install to ~/bin/ with:
```
	make install
```
And clean the directory with:
```
	make clean
```
Running, for example, a vesicle simulation can be performed by:
```
	liposome testLipo $RANDOM 80000 3.45
	MD testLipo
```
Running on a slurm cluster can be performed by creating a submission script, such as submit.sh in this project, and running:
```
	cp submit.sh test/
	cd test
	sbatch submit.sh
```

Alternatively, if you make MD.cpp, or any other executable, on its own:
```
	g++ -O3 -fopenmp MD.cpp -o MD
	mkdir ~/bin/
	cp MD ~/bin/
```

Compiling MDcuda.cu on University of Memphis cluster (this one only implements tension, bonding, bending, and non-bonded interactions):
```
	module load cuda11.2/toolkit
	make cuda
```

The MDcuda executable can be used in place of the MD executable (with no other modifications).
	

# Top (./) directory:
## ./MD.cpp: 
The simulation solver. Runs with openmp or single threaded. Makefile contains build options for the intel (icc) compiler version. Using this is easy if you generate an mpd file in a directory (like test):
```
	make
	make install
	cd generate
	make liposome
	cp liposome ~/bin
	cd ../test
	liposome testLipo $(RANDOM) 5000 3.45
	MD testLipo
```


## ./generate directory:
Generates a variety of configurations (.mpd files). Executables include:
- bendComposition: Creates a flat bilayer with variable lipid length and 2 lipid types.
- bilayer: Creates a flat bilayer.
- build: Creates a flat bilayer.
- continuumSphereAndBilayer: Creates a flat bilayer with a point like nanoparticle.
- continuumSphereAndLiposome: Creates a spherical bilayer with a point like nanoparticle.
- continuumSphereAndLiposomeWrapped: Creates a fully liposome wrapped point like nanoparticle.
- continuumSphereGrid: Creates a flat bilayer with a lattice of point like nanoparticles.
- continuumSphereSurfactant: Creates a flat bilayer with surfactant particles that adhear to a point like nanoparticle and bilayer.
- domain: Creates a flat bilayer with 2 lipid types in a domain.
- gas: Creates a gas of particles.
- lipoCyto: Creates a spherical bilayer with a cytoskeleton.
- liposome: Creates a spherical bilayer.
- nanoSphereAndBilayer: Creates a solid nanoparticle and a flat bilayer.
- nanoSphereGrid: Creates a grid of solid nanoparticles and a flat bilayer.
- nanoSphereModBilayer: 
- rigidNano: Creates a solid nanoparticle.
- scl: simple cubic lattice. Mostly for example.
- setMDConstants: Set molecular dynamics constants. Mostly for example.
- setMDConstants_UminAH: Set molecular dynamics constants. Mostly for example.
- solventBilayer: Creates a flat bilayer with solvent.
- solvent: Creates a gas of solvent.
- tetrahedron: Creates a tetrahedron. Mostly for example.
- tube: Creates a tubular bilayer.

## ./analysis directory:
Do analysis and modifications of systems (with .mpd files). Executables include:
- addBond: Add a bond to a system.
- adjustConstants: Adjust constants of a system.
- adjustMolecule: Adjust a moleucle of a system.
- boundaryPotential: Add a boundary potential to a system.
- chainBend: Calculates bending angle of lipids.
- chainBendLeaflet: Calculates bending angle of lipids per leaflet.
- changeDeltaT: Adjust timestep accuracy.
- changeFrame: Modify current system's frame.
- changeSeed: Modify the system's seed.
- changeTemperature: Modify the system's temperature.
- changeTension: Modify the system's tension.
- changeTime: Modify the system's start and end times.
- cluster: Calculates cluster sizes.
- clusterDomains: Calculates number of cluster domains.
- curvature3D_v7: Calculates curvature energy of a bilayer system.
- cycleXYZ: Cycle through a system's xyz file. Mostly for example.
- cylinderFluctuation: Calculates cylinder fluctuations?
- diffusion: Calculates diffusion. Make sure system store step is short enough to prevent jumps.
- duplicateImage: Duplicates a system across a dimension.
- extractXYZ: Extract the current XYZ position's of a system.
- forceBeadBead: Generates a graph of point like nanoparticles to other point like nanoparticles.
- forceBead: Generates a graph of point like nanoparticles to normal particles.
- frameReadPerSecond: Cycle through XYZ file. Mostly for example. Includes a timer.
- framesBeadParticleCount: Counts how many beads of a particular type there are in a system.
- frames: Cycle frames? Mostly example.
- framesDensity: Calculates density of a frame.
- framesDistanceBetweenIndices: Distances between two particle indices by frame.
- framesFlipFlop: Calculates lipid flip flop through by frame.
- framesPotentialTypePairs: Calculates potential between two types of particles.
- framesPotentialTypePairsList: Calculates potential between two types of particles and dumps a list of those.
- heightMap: Generates a height map from a flat system along z direction.
- heightMapFrames: Generates a height map from a flat system along z direction for every frame.
- lsPotentials: Dumps the non-bonded potential matricies of a system.
- lsProjection: Dumps a graphical projection of a system along 3 primary axii.
- modBond: Modify the a particular bonding interaction.
- modContinuumSphere: Modify a particular point like nanoparticle interaction.
- modContinuumSphereLong: Modify a particular point like nanoparticle interaction along with range (r_c).
- modSize: Modify a system's size.
- momentOfInertia: Calculation the moment of inertia of a system by frame.
- momentOfInertiaXYZ: Generate the axii for the moment of inertia of a system by frame.
- MSD: Mean square displacement. Diffusion.
- nanoPotential: Nanoparticle potential of a solid nanoparticle.
- nanoVsCOM: Nanoparticle center of mass distances by frame.
- neckBleb: Calculates the size of a neck near a bleb.
- netVel: Calculates the net velocity of a system.
- NPAngleFrames: Calculates the angle between point like nanoparticles and center of mass of a liposomes.
- nShapes: Counts the number of blebs in a system. Could be used to count exclusions from some particle.
- particleDistances: Distances between two particles in a system.
- potentialBeadBead: Makes a graph of potential between two point like nanoparticles.
- potentialBead: Makes a graph of potential between a point like nanoparticle and another particle.
- potential: Makes a graph of potential between particles.
- profile: Slices a system along one of the three primary axii.
- projectedCurvature: Projected curvature per layer. Use curvature_v7 instead.
- radial: Radial profile of other particles near a given particle.
- recenterFrames: Recenter the frames around a type, index, or molecule.
- removeSolvent: Remove solvent from a system's XYZ file.
- rollProfile: Obtains the center of mass and rolling angle of a solid nanoparticle.
- rotateXYZ: Rotate an xyz file along an axis.
- rProfile: Radial profile density around a cluster of particles in a system.
- rProfileFrames: Radial profile density around a cluster of particles in a system's XYZ file.
- singlePotential: Non-bonded potential of a single particle.
- singlePotentialType: Non-bonded potential of a single particle type.
- singlePotentialTypeHist: Histogram of a non-bonded potential of a single particle type.
- singlePotentialTypePairs: Non-bonded potential between two particle types.
- spherocity: Calculates the root mean square of a liposome.
- structureFactor: Calculates the structure factor of a flat bilayer.
- structureFactorFrames: Calculates the structure factor of a flat bilayer by frame.
- surface: Counts the number of surfaces in a system of types seperated by some cutoff.
- threeBodyDistribution: Average angular distribution of molecules.
- trajectoryFrame: Dumps trajectories by frame. Uses frames_name.xyz.
- translate: Move a system by some offset.
- twoBodyDistribution: Average distance distribution of molecules.
- typeDistribution: Calculates radial histogram of certain types along surfaces.
- vesicleRadius: Calculates radius of vesicle by frame.
- zPotential: Non-bonded potential histogram along z direction.
- zPotentialFrames: Non-bonded potential histograms along z direction for frames.
- zProfile: Histogram of particle probabilities. Bilayer profile.


## ./analysis/asciiMod directory:
Do calculations on plain text files. Executables include:
- abscissas: Generates the abscissas of a set.
- autoCorrelation: Auto-correlation of a variable.
- averageSets: Average a group of sets with a prefix.
- collect: Collect a group of sets with a prefix.
- gaussianKernelClamped: Gaussian kernal smoothing of a graph with a clamp at the average.
- gaussianKernel: Gaussian kernal smoothing of a graph.
- gaussianKernelProj: Gaussian kernal smoothing of a graph with a clamp at the average and zeroed.
- runningAvg: Running average of sets.
- skipError: Set error to zero every Nth data point.
- splitSets2: Split sets into multiple graphs.
- splitSets: Split sets into multiple graphs.
- sumSeries: Sum a series of sets into a single graph.
- sumSets: Bin values by the x column into averages of the y column.
- wham2: Weigted Histogram Analysis Method for free energies with an offset.
- wham: Weigted Histogram Analysis Method for free energies.


## ./analysis/asciiMod directory:
Do a function style calculation on command line input. Executables include:
- isoApex: Calculates the angle between two legs within a given stride length. Useful for continuumSphereAndVesicle.

## ./analysis/xyzMod directory:
Do operations on xyz files. Executables include:
- frameSkip: Only output every Nth frame.
- framesSeperate: Seperate frames into seperate XYZ files.
- move: Move all the particles by some offset.
- noType: Remove the type and header information. Generates x y and z coordinates only.
- oneTypeCount: Count all of one type.
- oneType: Output only one type.
- oneTypeFrames: Output only one type for frames.
- opticalReduction: Removes hidden particles?
- radiusTube: Radial density of a tube.
- rmFrames: Removes certain frames.
- selectFrames: Select a frame for output.
- selectRange: Select a range of frames for output.
- slice: Slices an XYZ file along an axis.
- sortFrames: Sort particles in frame by type. Useful if particle types change frequently.
- sphereMask: Masks a spherical range for output.
- target: Rescale an XYZ file from another format (ddscat?).
- xyzToBinary: Turns the XYZ file into a binary file. Can shrink files, but need binaryToXyz to make this work properly.


## ./include/system.h
This file includes just one class: Blob. This is the molecular dynamics system blob of parameters and data. It works with Revalee's, et al., model. Usually, you will create a blob with a real/floating type:
```
	Blob<double> system;
```
After which, you would populate it by simply setting values, with the setNAME() or addNAME() member functions, or by a script. Once they are set, you can read the values back with the readNAME() or getNAME() functions. Set and read prefixes indicate single valued parameters. Add and get prefixes indicate multi-valued parameters, and get usually returns a pointer. Getting the size of a system is as simple as:
```
	threeVector<double> size=system.readSize();
```
The important member functions are as follows:
### Parameters
- Blob(); //Constructor called everytime you create a new system.
- molecule<T, fourVector<int> > * addMolecule(molecule<T, fourVector<int> > &value); //Add a molecule residue to the system.
- molecule<T, fourVector<int> > * setMolecule(molecule<T,fourVector<int> > &value, int index); //Set or change a molecule residue in the system.
- molecule<T, fourVector<int> > * delMolecule(int index); //Delete a molecule residue in the system.
- molecule<T, fourVector<int> > * getMolecule(); //Return the pointer to the molecule residues in the system.
- void addParticle(position<T> pos, threeVector<T> vel, threeVector<T> acc); //Add a particle to the system. Includes the whole initial state for equations of motion.
- void setParticle(int index, position<T> pos, threeVector<T> vel, threeVector<T> acc); //Set or change a particle in the system.
- void delParticle(int index); //Delete a particle in the system.
- void allocParticle(int n); //Allocate some particles in the system. This is here to reduce excessive allocations.
- position<T> * getPositions(); //Get the pointer to the positions of particles in the system. 
- threeVector<T> * getVelocities(); //Get the pointer to the velocities of the particles in the system.
- threeVector<T> * getAccelerations(); //Get the pointer to the accelerations of the particles in the system.
- T * addTwoBodyFconst(T value); //Add a non-bonded force constant to the system.
- void setTwoBodyFconst(int index, T value); //Set or change a non-bonded force constant in the system.
- T * delTwoBodyFconst(int index); //Delete a non-bonded force constant in the system.
- T * getTwoBodyFconst(); //Get a pointer to the non-bonded force constants in the system.
- T * addTwoBodyUconst(T value); //Add a non-bonded potential constant to the system.
- void setTwoBodyUconst(int index, T value); //Set or change a non-bonded potential constant in the system.
- T * delTwoBodyUconst(int index); //Delete a non-bonded potential constant in the system.
- T * getTwoBodyUconst(); //Get a pointer to the non-bonded potential constants in the system.
- int readNParticles();//Read the number of particles in the system.
- int readNMolecules();//Read the number of molecular residues in the system.
- T readGamma();//Read gamma for the Langevin thermostat.
- T readInitialTemp();//Read the initial temperature of the system, in reduced units.
- T readFinalTemp();//Read the final temperature of the system, in reduced units.
- int readSeed();//Read the random number seed of the system.
- int readNTypes();//Read the number of particle types in the system.
- threeVector<bool> readPeriodic();//Read the periodic boundaries flag in the system.
- T readCutoff();//Read the force and potential cutoff of the non-bonded interactions.
- threeVector<T> readSize();//Read the size of the system, in reduced units.
- T readInitialTime();//Read the initial time of the system, in tau.
- T readFinalTime();//Read the ending time of the system, in tau.
- T readDeltaT();//Read the timestep of the system, in tau.
- T readStoreInterval();//Read the store interval of the system, in tau.
- T readMeasureInterval();//Read the measure interval of the system, in tau.
- T readDeltaLXY();//Read the maximum step size of the constant pressure length change, in reduced units.
- T readRemoveSolvent();//Read the amount of solvent to remove, as a fraction of total solvent.
- T readTempStepInterval();//Read temperature step interval.
- twoVector<T> readSolventGamma();//Read solvent gamma.
- T readGammaType(int type);//Read gamma type.
- T readTension();//Read constant tension parameter.
- void setGamma(T value);//Set or change gamma.
- void setInitialTemp(T value);//Set or change the initial temperature, in reduced units.
- void setFinalTemp(T value);//Set or change the final temperature, in reduced units.
- void setSeed(int value);//Set or change the random number seed.
- void setNTypes(int value);//Set or change the number of particle types.
- void setPeriodic(threeVector<bool> value);//Set or change the periodic boundaries flag.
- void setCutoff(T value);//Set or change the non-bonded cutoff.
- void setSize(threeVector<T> value);//Set or change the system size, in reduced units.
- void setInitialTime(T value);//Set or change the initial time, in tau.
- void setFinalTime(T value);//Set or change the final time, in tau.
- void setDeltaT(T value);//Set or change the timestep, in tau.
- void setStoreInterval(T value);//Set or change the store interval, in tau.
- void setMeasureInterval(T value);//Set or change the measure interval, in tau.
- void setDeltaLXY(T value);//Set or change the maximum step size of the constant pressure length change, in reduced units.
- void setRemoveSolvent(T value);//Set or change the fraction of solvent removal, as a fraction of total solvent.
- void setTempStepInterval(T value);//Set temperature step interval.
- void setSolventGamma(twoVector<T> value);//Set solvent gamma.
- void setTension(T value);//Set constant tension parameter
### Functions
- void doForce(int i, int j, threeVector<T> &minImg);//Non-bonded force function between particles i and j. Utilizing the minimum image for boundaries.
- T doPotential(int i, int j, threeVector<T> &minImg);//Non-bonded potential function between particles i and j. Utilizing the minimum image for boundaries.
- void doChainForce(int i);//Do CHAIN force bonded interactions on residue i. CHAIN is a combination of two and three body bonded forces.
- void doTorsionForce(int i);//Do TORSION force bonded interactions on residue i.
- void doBondForce(int i);//Do BOND force bonded interactions on residue i.
- void doBendForce(int i);//Do BEND force bonded interactions on residue i.
- void doBeadForce(int i);//Do BEAD force bonded interactions on residue i.
- void doBoundaryForce(int i);//Do BOUNDARY force on some particles in residue i
- fourVector<T> doBeadNeighbors(int i);//Extract neighbor info
- std::vector<threeVector<int> > doBeadChain(int i); //look for bead chains
- T doChainPotential(int i);//Return the CHAIN potential for residue i. 
- T doDihedralPotential(int i);//Return the DIHEADRAL potential for residue i.
- T doTorsionPotential(int i);//Return the Torsion Potential for residue i.
- T doBondPotential(int i);//Return the Bond potential for residue i.
- T doBendPotential(int i);//Return the Bend potential for residue i.
- T doBeadPotential(int i);//Return the Bead potential for residue i.
- T doBoundaryPotential(int i);//Return the BOUNDARY potential on some particles in residue i.
- T doChainDPotential(int i, threeVector<T> aSize);//Return the change in CHAIN potential for residue i when the size changes to a new state.
- T doDihedralDPotential(int i, threeVector<T> aSize);//Return the change in DIHEADRAL potential for residue i when the size changes to a new state.
- T doTorsionDPotential(int i, threeVector<T> aSize);//Return the change in Torsion potential for residue i when the size changes to a new state.
- T doBondDPotential(int i, threeVector<T> aSize);//Return the change in Bond potential for residue i when the size changes to a new state.
- T doBendDPotential(int i, threeVector<T> aSize);//Return the change in Bend potential for residue i when the size changes to a new state.
- T doBeadDPotential(int i, threeVector<T> aSize);//Return the change in Bead potential for residue i when the size changes to a new state.
