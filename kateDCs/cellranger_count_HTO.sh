#!/bin/bash
#SBATCH -q long          # submit to long queue
#SBATCH -J cellRanger            # job name
#SBATCH -p  himem #highmemory nodes

DIR_OUT='SAM24427129/'
FASTQ='../Libraries_SAM24427129.csv'
AB='../NGS4936_BarcodeFeature.csv'
DIR_DB='/gstore/project/moussion_lab/kate/GRCm38.transgenes'
DIR_CR='/gne/data/dnaseq/analysis/aplle/software/10x/cellranger-7.0.1/cellranger'
OUT_NAME='SAM24427129_CellRanger'

mkdir $DIR_OUT
cd $DIR_OUT

$DIR_CR count --id=$OUT_NAME --libraries=$FASTQ --feature-ref=$AB --transcriptome=$DIR_DB --localcores=12 --chemistry=threeprime
