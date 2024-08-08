CC=g++
CPPFLAGS=-O3 -std=c++11
EXEC= addBond adjustConstants adjustMolecule boundaryPotential \
chainBend chainBendLeaflet changeDeltaT changeFrame changeSeed changeTemperature \
changeTension changeTime cluster clusterDomains curvature3D_v7 cycleXYZ \
cylinderFluctuation diffusion duplicateImage extractXYZ forceBeadBead forceBead \
frameReadPerSecond framesBeadParticleCount frames framesDensity \
framesDistanceBetweenIndices framesFlipFlop framesPotentialTypePairs \
framesPotentialTypePairsList heightMap heightMapFrames lsPotentials \
lsProjection modBond modContinuumSphere modContinuumSphereLong modSize \
momentOfInertia momentOfInertiaXYZ MSD nanoPotential shellPotential nanoVsCOM neckBleb netVel \
NPAngleFrames nShapes particleDistances potentialBeadBead potentialBead potential \
profile projectedCurvature radial recenterFrames removeSolvent rollProfile \
rotateXYZ rProfile rProfileFrames singlePotential singlePotentialType \
singlePotentialTypeHist singlePotentialTypePairs spherocity structureFactor \
structureFactorFrames surface threeBodyDistribution trajectoryFrame translate \
twoBodyDistribution typeDistribution vesicleRadius zPotential \
zPotentialFrames zProfile widthFrames avgWidthFrames structureFactorEveryFrame

TARGETS: $(EXEC)

ifeq ($(PREFIX),)
PREFIX=$(HOME)
endif

.PHONY: all

all: $(TARGETS)

cluster:
	$(CC) -o cluster cluster.cpp $(CPPFLAGS) -fopenmp

nanoPotential:
	$(CC) -o nanoPotential nanoPotential.cpp $(CPPFLAGS) -fopenmp

shellPotential:
	$(CC) -o shellPotential shellPotential.cpp $(CPPFLAGS) -fopenmp

structureFactorFrames:
	$(CC) -lfftw3 -o structureFactorFrames structureFactorFrames.cpp $(CPPFLAGS)

structureFactorEveryFrame:
	$(CC) -lfftw3 -o structureFactorEveryFrame structureFactorEveryFrame.cpp $(CPPFLAGS)

$(TARGET):
	$(CC) -o $(TARGET) $(TARGET).cpp $(CPPFLAGS)

#buildOpenMP: $(TARGET).cpp
#	$(CC) -o $(TARGET) $(TARGET).cpp $(CPPFLAGS)

.PHONY: clean

clean:
	rm -f $(EXEC)

install: $(EXEC)
	install -d $(PREFIX)/bin/
	install -m 755 $(EXEC) $(PREFIX)/bin/

#test:
#	./test.sh
