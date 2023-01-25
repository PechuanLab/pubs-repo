#!/bin/bash
#SBATCH -p defq
#SBATCH --qos=medium
#SBATCH -c 16
#SBATCH --mem=32G

echo "sample: $sample, fastqs: $fastqs"
cellranger count \
--id=$sample \
--transcriptome=#Transcriptome \
--fastqs=$fastqs \
--sample=$sample \
--localcores 16 \
--localmem 32
