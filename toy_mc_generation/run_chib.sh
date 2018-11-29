#!/bin/bash
# $1: 1=do NOT apply selection
timestamp=`date +"%Y%m%d%H%M"`
code_base_dir="/afs/hephy.at/user/j/jnecker/code/chib_chic_polFW/toy_mc_generation"
tmp_out_base_dir="/afs/hephy.at/work/j/jnecker/data/chib_results/masterthesis/toymc/${timestamp}"
n_events=25000000

no_selection=0
if [[ $# -gt 0 ]];then no_selection=$1; fi
suffix="WITHselection"
if [[ $no_selection -eq 1 ]]; then suffix="NOselection"; fi
echo $suffix

mkdir -p $tmp_out_base_dir

echo "Generating chi1"
for CHI_R in "0" "1" "2/3"; do
	outfile="${tmp_out_base_dir}/toymc_chib-R_${CHI_R////o}_${suffix}.root"
	echo $outfile
	sbatch ${code_base_dir}/generate_chibpol.sh $outfile $n_events 1 $CHI_R 0 $no_selection
done

echo "Generating chi2"
for CHI2_R in "0,0" "1,0" "0,1" "2/5,2/5"; do
    CHI2_R1=$(echo ${CHI2_R} | awk -F',' '{print $1}')
    CHI2_R2=$(echo ${CHI2_R} | awk -F',' '{print $2}')
	outfile="${tmp_out_base_dir}/toymc_chib-R1_${CHI2_R1////o}_R2_${CHI2_R2////o}_${suffix}.root"
	echo $outfile
	sbatch ${code_base_dir}/generate_chibpol.sh $outfile $n_events 2 $CHI2_R1 $CHI2_R2 $no_selection
done
