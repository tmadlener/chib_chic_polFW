#!/bin/bash

NGEN_GEN=$(echo '5 * 10000000' | bc) # Number of generated events at gen level
NGEN_RECO=$(echo '10 * 10000000' | bc) # number of generated reco level events

gen_scale=$(echo ${NGEN_RECO} / ${NGEN_GEN=} | bc -l)

bin_var="JpsiPt"
binning="8:20,5" # four bins in pT

dir=${WORK}/ChiPol/chic2_chic1_ratios/toy_mc_data_gen_new/acceptance_correction_maps_data/gen_HX/
exe=${CHIB_CHIC_POLFW_DIR}/polarization/make_corr_maps.py


for gen_file in $(ls ${dir}/gen_level/*/toy_data.root); do
    pol_scen=$(echo ${gen_file} | awk -F'/' '{print $(NF-1)}')
    reco_file=${dir}/reco_level/${pol_scen}/toy_data.root
    out_file=${dir}/acc_maps_${pol_scen}.root

    echo "processing "${pol_scen}

    python ${exe} --bin-variable ${bin_var} --binning ${binning}\
           --scale-gen ${gen_scale}\
           ${gen_file} ${reco_file} ${out_file}
done
