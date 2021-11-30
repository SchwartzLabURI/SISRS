#!/bin/bash
#SBATCH --job-name="chloroplast1"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE TO user email address
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

#CHANGE if not on URI system
module purge

module load BLAST+/2.10.1-gompi-2020a
module load Bowtie2/2.4.4-GCC-10.2.0
module load Python/3.8.6-GCCcore-10.2.0

git clone https://github.com/Kinggerm/GetOrganelle.git

git clone https://github.com/Kinggerm/GetOrganelleDB

python3 ./GetOrganelle/Utilities/get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./GetOrganelleDB/0.0.0
