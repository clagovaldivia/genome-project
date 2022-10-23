#!/bin/bash -l
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J met
#SBATCH --mail-type=ALL
#SBATCH --mail-user claudia.gonzalez-valdivia.1872@student.uu.se

export DIR="/home/claudg/genome-project/analyses/5_DNA_bins_prokka"

cd $DIR
list=(1 2)
for A in ${list[@]}
do	
	cd $A
	bins=(1 2 3 4 5 6 7 8 9 10 11 12)
	for x in SRR4342129_bins.$bins.fa.out
 	do
		date=10112022
		grep "eC_number=" PROKKA_${date}.gff | cut -f9 | cut -f1,2 -d ';'| sed 's/ID=//g'| sed 's/;eC_number=/\t/g' > PROKKA.$Sbins.ec
	done
done

for a in 2
do
	bins=(13 14 15 16 17 18 19 20 21 22)
	for a in S*.$bins.*faa
