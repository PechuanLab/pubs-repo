#!/bin/bash

#SBATCH -q  long          # submit to long queue
#SBATCH -J  CellBender    # job name
#SBATCH -p  himem #highmemory nodes

cellbender remove-background \
     --input /gstore/scratch/u/pechuanj/NGS4936_KateDCs/SAM24427126/SAM24427126_CellRanger/outs/raw_feature_bc_matrix \
     --output ./SAM24427126_CellBender.h5 \
     --expected-cells 4295 \
     --total-droplets-included 15000
