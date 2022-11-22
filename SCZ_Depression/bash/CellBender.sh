#!/bin/bash

#SBATCH -q  long          # submit to long queue
#SBATCH -J  CellBender    # job name
#SBATCH -p  himem #highmemory nodes

cellbender remove-background \
     --input /gstore/scratch/u/pechuanj/NGS2180_Brain/BrainJasmine/LIB3598482_HITS3419835/outs/raw_feature_bc_matrix \
     --output ./HITS3419835.h5 \
     --expected-cells 3422 \
     --total-droplets-included 12000
