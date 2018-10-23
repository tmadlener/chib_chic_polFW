#!/bin/bash

NGEN_GEN=10000000 # Number of generated events at gen level
NGEN_RECO=$(echo '10 * 100000000' | bc) # number of generated reco level events

dir=${WORK}/ChiPol/chic2_chic1_ratios/toy_mc_data_gen_new/acceptance_correction_maps_data/
exe=${CHIB_CHIC_POLFW_DIR}/polarization/make_corr_map.py


for gen_file in $(ls ${dir}/gen_level/*/toy_data_????????.root); do
    pol_scen=$(echo ${gen_file} | awk -F'/' '{print $(NF-1)}')
    reco_file=${dir}/reco_level/${pol_scen}/toy_data.root
    out_file=${dir}/acc_maps_${pol_scen}.root

    echo "processing "${pol_scen}
    python ${exe} ${gen_file} ${reco_file} ${out_file}
done
