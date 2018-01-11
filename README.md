# chi\_b and chi\_c polarization framework

## Intent
This repository is intended to hold all the necessary code for doing a chi\_b and chi\_c polarization analysis. The preparation (i.e. mass- and/or lifetime fits) steps are separated while the later analysis is shared between the two states.

## Usage
Some of the plotting code currently requires **root_pandas**, which is most easily setup by using the environment from a **CMSSW** release later than **CMSSW_9_4_X**.

To setup anything additional do `source .setup.sh` in the base-directory of this repository.

## TODO list
* [ ] Move existing code from other repositories into this one
  - [x] chi\_c tupling (starting from finished mass- and lifetime fit)
  - [ ] existing code for plotting
  - [ ] chi\_b mass fit framework
* [ ] Generalize plotting code to handle both cases
* [ ] Make setup script setting all the necessary environment and `PATH` variables so that no additional setup from CMSSW is necessary.
