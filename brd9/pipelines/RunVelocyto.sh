#!/bin/bash

#SBATCH -q  long          # submit to long queue
#SBATCH -J  Velocyto    # job name
#SBATCH -p  himem #highmemory nodes

velocyto run10x /path/to/CellRanger/LIB5435059_SAM24392038 GRCm38/genes/genes.gtf
