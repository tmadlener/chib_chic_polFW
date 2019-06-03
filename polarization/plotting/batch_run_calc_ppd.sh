#!/bin/bash
#SBATCH -J CalcPPD
#SBATCH -D /afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/data_fits/simultaneous_final_fits/ppd_combination
#SBATCH -o /afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/data_fits/simultaneous_final_fits/ppd_combination/batch_logfiles/run_calc_ppd_%A_%a.out

source ${CHIB_CHIC_POLFW_DIR}/scripts/bash_functions.sh

set -x
set -e

ratiofile=${1}
outfile=${2}
nevents=${3}

calc_ppd=${CHIB_CHIC_POLFW_DIR}/polarization/plotting/calc_ppd.py

outdir=$(dirname ${outfile})

run_sandboxed ${outdir} ${calc_ppd} --number-extractions ${nevents} --hists-only \
              ${ratiofile} $(basename ${outfile})
