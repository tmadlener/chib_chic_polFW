#!/bin/bash
timestamp=`date +"%Y%m%d_%H%M%S"`
code_base_dir="/afs/hephy.at/user/j/jnecker/code/chib_chic_polFW/toy_mc_generation"
#out_base_dir="/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc"
out_base_dir="/scratch/jnecker/toymc"
n_events=10000000

echo "Generating chi1"
for CHI_R in "0" "1" "0.666"; do
	outfile="${out_base_dir}/toymc_chib_R-${CHI_R/./p}_${timestamp}.root"
	echo $outfile
	sbatch ${code_base_dir}/generate_chibpol.sh $outfile $n_events 1 $CHI_R 0
done

echo "Generating chi2"
for CHI2_R in "0,0" "1,0" "0,1" "0.4,0.4"; do
    CHI2_R1=$(echo ${CHI2_R} | awk -F',' '{print $1}')
    CHI2_R2=$(echo ${CHI2_R} | awk -F',' '{print $2}')
	outfile="${out_base_dir}/toymc_chib_R1-${CHI2_R1/./p}_R2-${CHI2_R2/./p}_${timestamp}.root"
	echo $outfile
	sbatch ${code_base_dir}/generate_chibpol.sh $outfile $n_events 2 $CHI2_R1 $CHI2_R2
done