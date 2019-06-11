#!/bin/bash
## setup script for setting up some environment and path variables.
## TODO: Find common setup, that works with all the python packages and ROOT etc. and put the initialization here.

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


    ## get the bash helper function definitions
    source ${CHIB_CHIC_POLFW_DIR}/scripts/bash_functions.sh

    fetch_json_header ${@}

    build_shapes

    if $(check_args_flag "--run-tests" ${@}); then
        run_tests
    fi

    ## Set the latex command
    ## Check if xelatex exists and use it, otherwise default to pdflatex
    if type xelatex > /dev/null 2>&1; then
        export LATEX_EXE=xelatex
    else
        export LATEX_EXE=pdflatex
    fi

    # create a random number that can be used to shift the fitted value of delta lambda to avoid unblinding
    export RANDOM_DELTA_LAMBDA_SHIFT=$(get_rand 271828182845904)
fi
