module load BLAST+/2.11.0-gompi-2020b
module load Bowtie2/2.4.4-GCC-10.2.0

git clone https://github.com/Kinggerm/GetOrganelle.git

git clone https://github.com/Kinggerm/GetOrganelleDB

python3 ./GetOrganelle/Utilities/get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./GetOrganelleDB/0.0.0
