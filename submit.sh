#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=400M

# Generate a liposome system
# Usage: liposome name seed nLipids arialDensity
liposome testLipo $RANDOM 80000 3.45

# Run the liposome system
# Usage: MD name
MD testLipo
