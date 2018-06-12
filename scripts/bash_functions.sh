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
## NOTE: This uses cleanup_or_exit so this should only be used in scripts as otherwise you will be thrown out
## of your shell if the command you want to sandbox fails
function run_sandboxed() {
    sandboxdir=${1}
    orig_exe=${2}
    shift 2
    args=${@}

    curr_dir=$(pwd)
    mkdir -p ${sandboxdir}

    # copy the executable into the sandbox directory and make it unique
    exe=$(basename ${orig_exe})_$(rand_str 16)
    cp ${orig_exe} ${sandboxdir}/${exe}
    cd ${sandboxdir}

    print_date "start of "${exe}
    ./${exe} ${args}
    cleanup_or_exit $? ${exe} ${curr_dir}
    print_date "end"
    cd ${curr_dir}
}
