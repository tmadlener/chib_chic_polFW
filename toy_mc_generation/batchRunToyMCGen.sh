#!/bin/bash
#SBATCH -J ToyMCGeneration
#SBATCH -D /afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/toy_mc_data_gen_new/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/ChiPol/chic2_chic1_ratios/toy_mc_data_gen_new/logfiles/runToyMCGen_%A.out

source ${CHIB_CHIC_POLFW_DIR}/scripts/bash_functions.sh

set -x

## input arguments
genfile=${1}
state=${2}
hel1=${3}
hel2=${4}
nevents=${5}
shift 5

set +x # To avoid excessive output from arg handling
CREATE_DISTS=0
if $(check_args_flag "--create-dists" ${@}); then
    CREATE_DISTS=1

    # remove the --create-dists argument, because the c++ will not handle it properly
    for arg; do
        shift
        [ ${arg} = "--create-dists" ] && continue
        set -- "$@" "$arg"
    done
fi

RM_SOURCE=0
if $(check_args_flag "--rm-source" ${@}); then
    RM_SOURCE=1

    for arg; do
        shift
        [ ${arg} = "--rm-source" ] && continue
        set -- "$@" "$arg"
    done
fi
set -x

echo "CREATE_DISTS: "${CREATE_DISTS}", RM_SOURCE: "${RM_SOURCE}


exe=${CHIB_CHIC_POLFW_DIR}/toy_mc_generation/run_chicpolgen

outdir=$(dirname ${genfile})

sbid=$(run_sandboxed ${outdir} ${exe} --genfile $(basename ${genfile}) --nevents ${nevents} \
                     --helicity1 ${hel1} --helicity2 ${hel2} --state ${state} ${@})

# create the histogram right after the generation
if [[ ${CREATE_DISTS} -eq 1 ]]; then
    mkhists=${CHIB_CHIC_POLFW_DIR}/toy_mc_generation/create_hists.py
    real_genfile=$(echo ${genfile} | sed 's:\.root:_'"$sbid"'\.root:')
    dist_file=$(echo ${real_genfile} | sed 's:\.root:_dists\.root:')
    ${mkhists} ${real_genfile} ${dist_file}
fi
if [[ ${RM_SOURCE} -eq 1 ]]; then
    # remove the generated file after the histograms have been created
    rm ${real_genfile}
fi
