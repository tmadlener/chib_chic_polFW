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
