#!/bin/bash

#SBATCH -q long          # submit to long queue
#SBATCH -J cellRanger            # job name
#SBATCH -p  himem #highmemory nodes

for line in `cat Libraries.csv`
do
        sample=`echo $line | cut -d, -f1`
        fastqs=`echo $line | cut -d, -f2`
        if [ -d $sample ] ; then
                echo "Directory for $sample already exists, not submitting"
        else
                echo "Will submit for $sample..."
                echo "Fastqs are $fastqs..."
                sbatch -J $sample -e ${sample}_err.txt -o ${sample}_out.txt --export=sample=$sample,fastqs=$fastqs CellRanger.sh
        fi
done
