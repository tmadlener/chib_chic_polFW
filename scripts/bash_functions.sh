#!/bin/bash

## define functions used commonly throughout bash scripts

## check if a given flag is in the passed arguments
## flag is first argument to function, rest is treated as arguments
function check_args_flag() {
    flag=${1}
    shift

    for arg in "${@}"; do
        if [[ ${flag} = ${arg} ]]; then
            return 0 # reversed world in bash ;)
        fi
    done
    return 1
}

## export function to have it available in other scripts
export -f check_args_flag

## chib preparation needs cpp json (header only) library. This function checks if the file exists and if not
## gets it from github. If it exists nothing is done unless the --update flag is passed
function fetch_json_header() {
    if [[ $# -eq 1 ]]; then
        if $(check_args_flag "--update" ${@}) || $(check_args_flag "-u" ${@}); then
            echo "Will download json library again even if already present"
            FORCE_DOWNLOAD=1
        fi
    fi

    ## NOTE: Assuming here that the base dir is already set (so call this after setting it!)
    json_header=${CHIB_CHIC_POLFW_DIR}/general/interface/json.hpp
    if [[ ! -f ${json_header} ]] || [[ ${FORCE_DOWNLOAD} -eq 1 ]]; then
        echo "JSON header library not present or download forced. Fetching it from github"
        wget -O ${json_header} https://github.com/nlohmann/json/releases/download/v3.0.1/json.hpp 2> /dev/null

        if [[ $? -ne 0 ]]; then
           echo "ERROR while getting JSON header"
        fi
    fi
}

## run all the available tests
function run_tests() {
    if [[ $# -eq 0 ]]; then
        python -m unittest discover ${CHIB_CHIC_POLFW_DIR}/python/test
    else
        python -m unittest ${@}
    fi

}

## get a variable length random string that can be used to name temp directories
function rand_str() {
    if [[ $# -eq 0 ]]; then
        length=16
    else
        length=$1
    fi

    echo $(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w ${length} | head -n 1)
}
export -f rand_str

## get the pt bin from the file-ending
## WARNING: NO CHECK OF ARGUMENTS
function get_pt() {
    file=$1
    pt_bin=${file: -6:1} # assuming that pT bin is only 1 digit here. TODO: make this more stable
    echo $pt_bin
}
export -f get_pt

## view the contents of passed pickle files that only use builtin objects and dump them to stdout
function pklview() {
    if [[ $# -eq 0 ]]; then
        echo "usage: pkl_view FILE [FILES]"
        return 64
    fi
    files=${@}
    for file in ${files[@]}; do
        if [[ ! -f ${file} ]]; then
            echo "${file} not found"
            return 1
        fi
        echo "${file}:"
        python -c "import pickle; import pprint; pprint.pprint(pickle.load(open('${file}')))"
    done
}
export -f pklview

## Convert old-style pkl files from fits to json files, so that they can be used by newer versions as well
function pkl2json() {
    if [[ $# -eq 0 ]]; then
        echo "usage: pkl2json PKLFILE [JSONFILE]"
        return 64
    fi
    if [[ $# -ge 1 ]]; then
        pklfile=$1
        shift
        if [[ $# -ge 1 ]]; then
            jsonfile=$2
        else
            jsonfile=$(echo ${pklfile} | sed 's/.pkl/.json/')
        fi
        if [[ ${jsonfile} = ${pklfile} ]]; then
            echo "input and output file are the same"
            return 1
        fi
    fi

    if [[ ! -f ${pklfile} ]]; then
        echo "${pklfile} not found"
        return 1
    fi

    python -c "import pickle; import json; json.dump(pickle.load(open('${pklfile}')), open('${jsonfile}', 'w'), indent=2)"
}
export -f pkl2json

## print the current date and the name of a time point
function print_date() {
    echo "-------------------- "${@}": " $(date) " --------------------"
}

## check if the last action (first argument) has been successful (exited with 0).
## If so, remove the executable (second argument).
## Else leave it in place and exit with the passed exit code.
## If a third argument is passed, the function will cd to the directory specified by that
function cleanup_or_exit() {
    if [ ${1} -eq 0 ]; then
        echo ${2} "exited with 0. Cleaning up executable"
        rm ${2}
    else
        echo ${2} "exited with "${1}". Leaving executable in place"
        print_date "end"
        if [ ${#} -eq 3 ]; then
            cd ${3}
        fi
        exit ${1} # for testing disable the exit
    fi
}

## Run a command in a "pseudo sandbox" by copying the executable (second) argument to the passed directory (first argument)
## and running it there. The rest of the arguments is passed on the the executable
## NOTE: The executable is only guaranteed to not leak out of the sandbox if it produces all of its output in the directory
## it is run in
## NOTE: This uses cleanup_or_exit so this should only be used in scripts as otherwise you will be thrown out
## of your shell if the command you want to sandbox fails
function run_sandboxed() {
    outdir=${1}
    orig_exe=${2}
    shift 2
    args=${@}

    # create a sandbox directory in the output directory
    sandboxid=$(rand_str 8)
    sandboxdir=${outdir}/${sandboxid}

    curr_dir=$(pwd)
    mkdir -p ${sandboxdir}

    # copy the executable into the sandbox directory and run there
    exe=$(basename ${orig_exe})
    cp ${orig_exe} ${sandboxdir}/${exe}
    cd ${sandboxdir}

    logfile=${sandboxdir}/sandbox_run.log

    print_date "start of "${exe}" in "${sandboxdir} >> ${logfile}
    ./${exe} ${args} >> ${logfile}
    cleanup_or_exit $? ${exe} ${curr_dir} >> ${logfile}

    print_date "end" >> ${logfile}

    # if we are still alive here, we have to move the output to the "real" output directory
    echo "moving outputs to "${outdir} >> ${logfile}
    for file in $(ls ${sandboxdir}); do
        mv ${sandboxdir}/${file} ${outdir}/$(echo ${file} | sed 's/\(\..*\)$/_'"$sandboxid"'\1/')
    done

    # echo "removing temp sandboxdir"
    cd ${outdir}
    if ! [[ $(ls -A ${sandboxdir}) ]]; then
        rmdir ${sandboxdir}
    fi

    cd ${curr_dir}

    # return the sandbox id in order to be able to identify which files have been created inside this sandbox
    echo "${sandboxid}"
}

## Function the builds the shared object (assuming a folder is present for it in the general directory)
## First and only argument is the name of the folder
function build_shared_object() {
    local object=${1}
    local cmpfile=$(mktemp ${CHIB_CHIC_POLFW_DIR}/general/${object}/compile.XXXXXXXX)
    make -C ${CHIB_CHIC_POLFW_DIR}/general/${object}/ > ${cmpfile} 2>&1
    local compile_status=$?

    if [[ ${compile_status} -ne 0 ]]; then
        echo "Problem building ""${object}"":"
        echo "================================================================================"
        cat ${cmpfile}
        echo "================================================================================"
        echo "Check output of compilation in: "${cmpfile}
    else
        rm ${cmpfile}
    fi

    return ${compile_status}
}

## Function that calls make in the RooDoubleCB directory so that the shared object is present
function build_shapes() {
    for shape in "RooDoubleCB" "RooErfExponential" "RooPowerlawExponential"; do
        build_shared_object $shape
    done
}

## Function that moves an already existing file to a new place
function mv_existing() {
    file=${1}
    if [ -f ${file} ]; then
        mvfile=$(echo ${file}).$(date '+%d%m%y_%H%M%S')
        echo "${file} already exists. Moving it to ${mvfile}"
        mv ${file} ${mvfile}
    fi
}
export -f mv_existing

## print the passed name and the value of the variable (if it is set)
function print_var() {
    var=${1}
    if [ -n "${!var+set}" ]; then # only print variables that are set (even if possibly empty)
        echo ${var}=${!var}
    fi
}
export -f print_var

## run the costh binned fits and also plot the results
## Arguments:
## 1 - input data file containing the data that are fitted
## 2 - output directory into which the results are stored
## 3 - config json file that contains the fit configuration
## 4 - the binvariable
## 5 - the binning
function run_plot_fits() {
    if [[ $# -lt 5 ]]; then
        echo "usage: run_plot_fits DATAFILE OUTPUTDIR CONFIGFILE BINVAR BINNING [MASSRANGE]"
        return 64
    fi

    local infile=$1
    local outdir=$2
    local config=$3
    local binvar=$4
    local binning=$5

    shift 5
    if [[ $# -ne 0 ]]; then
        local massrange="--massrange "${1}
    fi

    local FITTER=${CHIB_CHIC_POLFW_DIR}/polarization/chic_fitting/costh_binned_massfit.py
    local PLOTTER=${CHIB_CHIC_POLFW_DIR}/polarization/chic_fitting/make_fit_plots.py
    local REPORTGEN=${CHIB_CHIC_POLFW_DIR}/misc_scripts/make_fit_res_summary_binned.py
    local GRAPHPLOTTER=${CHIB_CHIC_POLFW_DIR}/misc_scripts/params_v_costh_plots.py

    local fitresfile=${outdir}/costh_binned_fit_results.root

    mkdir -p ${outdir}
    local logfile=${outdir}/run_fits.log

    echo $(date) > ${logfile}

    python $FITTER config --binning=${binning} --binvariable=${binvar} \
           --fix-shape ${massrange} ${infile} ${fitresfile} ${config} >> ${logfile} 2>&1

    echo $(date) >> ${logfile}

    # rename logfile for plotting
    logfile=${logfile/fits.log/plots.log}

    python $PLOTTER --config --configfile ${outdir}/fit_model.json --outdir ${outdir}/plots/ \
           --fix-shape --graphs --verbose ${fitresfile} > ${logfile} 2>&1


    python $GRAPHPLOTTER --outdir ${outdir}/plots --no-ratio ${outdir}/plots/free_fit_param_graphs.root


    mkdir -p ${outdir}/latex
    cd ${outdir}/latex

    local bin_info_file=../$(basename ${fitresfile}| sed 's/.root/_bin_sel_info.json/')

    python ${REPORTGEN} -o fit_report.tex ../plots/ ${bin_info_file} && \
        ${LATEX_EXE} -interaction=nonstopmode fit_report.tex > /dev/null 2>&1 && \
        rm fit_report.{log,aux}
    cd -

}
export -f run_plot_fits

function run_sim_fits() {
    if [[ $# -lt 3 ]]; then
        echo "usage: run_sim_fits DATAFILE OUTPUTDIR CONFIGFILE"
        return 64
    fi

    local infile=$1
    local outdir=$2
    local config=$3

    local FITTER=${CHIB_CHIC_POLFW_DIR}/polarization/chic_fitting/simultaneous_binned_fit.py
    local PLOTTER=${CHIB_CHIC_POLFW_DIR}/polarization/chic_fitting/plot_simultaneous_binned_fit.py
    local PARAMPLOTTER=${CHIB_CHIC_POLFW_DIR}/misc_scripts/make_proto_param_plots.py
    local REPORTGEN=${CHIB_CHIC_POLFW_DIR}/misc_scripts/make_fit_res_summary.py

    local fitresfile=${outdir}/simultaneous_fit_results.root

    mkdir -p ${outdir}
    local logfile=${outdir}/run_fits.log

    set -x
    echo $(date) > ${logfile}
    python $FITTER --verbosity 1 --outfile ${fitresfile} ${infile} ${config} >> ${logfile} 2>&1 || return 1
    echo $(date) >> ${logfile}

    logfile=${outdir}/run_plots.log
    python $PLOTTER --verbose --graphs --outdir ${outdir}/plots ${fitresfile} ${outdir}/fit_model.json > ${logfile} 2>&1 || return 1

    python $PARAMPLOTTER --outdir ${outdir}/plots --no-ratio ${outdir}/plots/proto_param_graphs.root

    set +x

    mkdir -p ${outdir}/latex
    cd ${outdir}/latex

    python $REPORTGEN --outfile fit_report.tex ../fit_model.json ../plots/ ../simultaneous_fit_results.root && \
        ${LATEX_EXE} -interaction=nonstopmode fit_report.tex > /dev/null 2>&1 && \
        rm fit_report.{aux,log}
    cd -
}
export -f run_sim_fits

function chi2prob() {
    local chi2=${1}
    local ndf=${2}

    python -c "from scipy.stats import chi2; print(chi2.sf("${chi2}","${ndf}"))"
}
export -f chi2prob

# Get a random number possibly using a seed
function get_rand() {
    if [[ $# -gt 0 ]]; then
        local seed=${1}
    else
        local seed="None"
    fi

    echo $(python -c "import random; random.seed("${seed}"); print(random.uniform(7, 100))")
}

# Get the path of a file or directory
# Takes as single argument a relative path to a file and returns the full path
# See: https://stackoverflow.com/a/3915420
function realpath() {
    echo "$(cd "$(dirname "$1")"; pwd -P)/$(basename "$1")"
}
export -f realpath

# Function to merge all root files from different batch jobs into one file and removes the input files
function merge_batch() {
    local basename=${1}
    hadd -f ${basename}.root ${basename}_????????.root && rm ${basename}_????????.root
}
export -f merge_batch

# Function to do apply corrections
function apply_corrections() {
    if [[ $# -lt 4 ]]; then
        echo "usage: apply_corrections FITDIR DATAFILE CORRMAPDIR OUTFIlE"
        return 64
    fi
    local fitdir=${1}
    local datafile=${2}
    local cmapdir=${3}
    local outfile=${4}

    local fitfile=${fitdir}/simultaneous_fit_results.root
    local model=${fitdir}/fit_model.json
    local chi1corrf=${cmapdir}/chic1_R_2o3/toy_data.root
    local chi2corrf=${cmapdir}/chic2_R1_2o5_R2_2o5/toy_data.root

    mkdir -p $(dirname ${outfile})
    python ${CHIB_CHIC_POLFW_DIR}/polarization/plotting/correct_ratio.py --outfile ${outfile} \
           ${fitfile} ${datafile} ${model} ${chi1corrf} ${chi2corrf}
}
export -f apply_corrections
