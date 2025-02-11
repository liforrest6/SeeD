#!/bin/bash -l

#SBATCH -D /group/runciegrp2/Projects/SeeD/Scripts
#SBATCH -o /home/fli21/slurm-log/clump-%j.txt
#SBATCH -e /home/fli21/slurm-log/clump-%j.txt
#SBATCH -J clumping
#SBATCH -t 24:00:00
#SBATCH --mem 40GB
#SBATCH -n 5
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

dir='/group/runciegrp2/Projects/SeeD/'

python3 ${dir}/Scripts/clumpAlternateAlleles-GenomeWide.py


# python3 ${dir}/Scripts/LD_analysis.py