#!/bin/bash

#SBATCH -q  long          # submit to long queue
#SBATCH -J  CellBender    # job name
#SBATCH -p  himem #highmemory nodes

velocyto run10x /gstore/scratch/u/pechuanj/NGS3422/CellRanger/LIB5435059_SAM24392038 /gne/data/dnaseq/analysis/aplle/genomes/refdata-cellranger-igis-4.0/GRCm38/genes/genes.gtf
