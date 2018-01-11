#!/bin/bash

## setup script for setting up some environment and path variables.
## TODO: Find common setup, that works with all the python packages and ROOT etc. and put the initialization here.

## get the absolute directory of this script, regardless of where it's called from
export CHIB_CHIC_POLFW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

## add the python directory to the python search path
## NOTE: 'cmsenv'-ing a CMSSW release completely resets the PYTHON_PATH, so that it might be necessary to call this script again after setting it up
PYTHONPATH=$PYTHONPATH:${CHIB_CHIC_POLFW_DIR}/python
