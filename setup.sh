#!/bin/bash
## setup script for setting up some environment and path variables.
## TODO: Find common setup, that works with all the python packages and ROOT etc. and put the initialization here.


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


# first check if setup already happened, and only do it if not or if --force flag is passed.
# In this way, this script can be sourced multiple times or from within bash scripts
if [[ -z "${CHIB_CHIC_POLFW_DIR+x}" ]] || $(check_args_flag "--force" ${@}); then

    ## get the absolute directory of this script, regardless of where it's called from
    export CHIB_CHIC_POLFW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

    ## add the python directory to the python search path
    ## NOTE: 'cmsenv'-ing a CMSSW release completely resets the PYTHON_PATH, so that it might be necessary to call this script again after setting it up
    PYTHONPATH=$PYTHONPATH:${CHIB_CHIC_POLFW_DIR}/python

    ## setup up a common root-version and a common gcc compiler so that the variables can then be used in makefiles, etc.
    ## NOTE: currently uses the default setup of the environment
    export ROOT_CONFIG_BIN=$(which root-config)
    export GCC_BASE_DIR=""


    fetch_json_header ${@}
fi
