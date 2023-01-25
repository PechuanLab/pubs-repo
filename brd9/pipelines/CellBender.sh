#!/bin/bash

#SBATCH -q  long          # submit to long queue
#SBATCH -J  CellBender    # job name
#SBATCH -p  himem #highmemory nodes

cellbender remove-background \
     --input /path/to/cellranger/LIB/outs/raw_feature_bc_matrix \
     --output ./LIB.h5 \
     --expected-cells 5491 \
     --total-droplets-included 15000
