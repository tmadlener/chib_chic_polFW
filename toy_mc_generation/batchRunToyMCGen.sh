#!/bin/bash
#SBATCH -J ToyMCGeneration
#SBATCH -D /afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/toy_mc_data_gen/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/toy_mc_data_gen/logfiles/runToyMCGen_%A.out

source ${CHIB_CHIC_POLFW_DIR}/scripts/bash_functions.sh

# set -x

## input arguments
genfile=${1}
state=${2}
hel1=${3}
hel2=${4}
nevents=${5}
shift 5

exe=${CHIB_CHIC_POLFW_DIR}/toy_mc_generation/run_chicpolgen

outdir=$(dirname ${genfile})

run_sandboxed ${outdir} ${exe} --genfile $(basename ${genfile}) --nevents ${nevents} \
              --helicity1 ${hel1} --helicity2 ${hel2} --state ${state} ${@}
