#!/bin/bash

# script to produce the distribution histograms
hist_producer=${CHIB_CHIC_POLFW_DIR}/toy_mc_generation/create_hists.py

OUTBASEDIR=${WORK}/ChiPol/chic2_chic1_ratios/toy_mc_data_gen/gen_pt_range
# OUTBASEDIR=${WORK}/ChiPol/chic2_chic1_ratios/toy_mc_data_gen/smearing

for f in ${OUTBASEDIR}/*/toy_data_*.root; do
# for f in ${OUTBASEDIR}/*/*/toy_data_*.root; do
    outfile=$(echo $f | sed 's:toy_data:distribution_hists:;s/_[0-9a-zA-z]\{8\}\.root/\.root/') ## removes the random string from generation

    echo "Producing histograms from "${f}
    ${hist_producer} ${f} ${outfile}
done
