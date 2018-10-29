#!/bin/bash

set -e

NGEN_GEN=$(echo   '5 * 10000000' | bc) # Number of generated events at gen level
NGEN_RECO=$(echo '10 * 10000000' | bc) # number of generated reco level events

gen_scale=$(echo ${NGEN_RECO} / ${NGEN_GEN=} | bc -l)

BIN_VAR="JpsiPt"
BINNING="8:20,5" # four bins in pT

DO_MAPS=1
DO_PLOTS=1
DO_RATIOS=1
DO_RATIO_PLOTS=1

dir=${WORK}/ChiPol/chic2_chic1_ratios/toy_mc_data_gen_new/acceptance_correction_maps_data/gen_HX/
out_base=acc_maps
out_dir=$(dirname ${dir}/${out_base})
exe=${CHIB_CHIC_POLFW_DIR}/polarization/make_corr_maps.py
ratio_exe=${CHIB_CHIC_POLFW_DIR}/misc_scripts/divide_all_hists.py
plot_exe=${CHIB_CHIC_POLFW_DIR}/misc_scripts/plot_all_hists.py


if [[ ${DO_MAPS} -eq 1 ]]; then
    for gen_file in $(ls ${dir}/gen_level/*/toy_data.root); do
        pol_scen=$(echo ${gen_file} | awk -F'/' '{print $(NF-1)}')
        reco_file=${dir}/reco_level/${pol_scen}/toy_data.root
        out_file=${dir}/${out_base}_${pol_scen}.root

        mkdir -p $(dirname ${out_file})
        echo "processing "${pol_scen}

        python ${exe} --bin-variable ${BIN_VAR} --binning=${BINNING}\
               --scale-gen ${gen_scale}\
               ${gen_file} ${reco_file} ${out_file}

    done
fi

if [[ ${DO_PLOTS} -eq 1 ]]; then
    for map_file in ${dir}/${out_base}_*.root; do
        python ${plot_exe} --filter="^proj_" --outdir=$(echo ${map_file} | sed 's/\.root//')_plots\
               --draw-opt="colz" ${map_file}
    done
fi

if [[ ${DO_RATIOS} -eq 1 ]]; then
    for chi1_pol in R_0 R_1; do
        unpol_file=${dir}/${out_base}_chic1_R_2o3.root
        pol_file=${dir}/${out_base}_chic1_${chi1_pol}.root
        outfile=${out_dir}/acc_map_ratios_chic1_${chi1_pol}_R_2o3.root

        python ${ratio_exe} --filter="^proj_acc_map" ${pol_file} ${unpol_file} ${outfile}
    done

    for chi2_pol in R1_0_R2_1 R1_0_R2_0; do
        unpol_file=${dir}/${out_base}_chic2_R1_2o5_R2_2o5.root
        pol_file=${dir}/${out_base}_chic2_${chi2_pol}.root
        outfile=${out_dir}/acc_map_ratios_chic2_${chi2_pol}_R1_2o5_R2_2o5.root

        python ${ratio_exe} --filter="^proj_acc_map" ${pol_file} ${unpol_file} ${outfile}
    done
fi

if [[ ${DO_RATIO_PLOTS} -eq 1 ]]; then
    for map_ratio_file in ${out_dir}/acc_map_ratios*.root; do
        python ${plot_exe} --outdir=$(echo ${map_ratio_file} | sed 's/\.root//')_plots\
               --draw-opt="colz" --zrange=0,2 ${map_ratio_file}
    done
fi
