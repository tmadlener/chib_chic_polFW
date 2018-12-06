# Toy MC generation

This folder contains (almost) everything that is necessary for generating Toy MC
chi_cJ events (J=0,1,2) with different polarization hypotheses. The only parts
that are not inside this folder are some general utilities like the `ArgParser`
that live in the [general/interface](/general/interface/) directory.

## Introduction and theory background

The Toy MC generation is implemented according to [PRD **83** 096001
(2011)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.83.096001) with
some additions to emulate reconstruction effects (at least on a statistical
level). The core of the original implementation (added in this
[commit](https://github.com/tmadlener/chib_chic_polFW/commit/fdd1a450fef8516904aaf60f064b8fab3efdeed8))
is still in place all of the changes that were later done are mainly concerned
with making the generation easier to configure without having to change `.C` or
`.cc` files and recompiling them first.

In general the generation of an event starts by randomly generating a chi\_c and
a J\psi mass from a Gaussian distribution with the central value at the PDG mass
and the natural width (also from PDG) of the corresponding particle. After the
mass the emission angle of the J/psi in the chi\_c restframe cosTheta and Phi
are generated randomly. Using these values the angular distribution W(cosTheta,
Phi | lambda) is evaluated, where the lambda parameters correspond to the input
settings. Depending on the value of W for a given event it is either accepted or
rejected in such a way that it is ensured that on a statistical basis the
accepted events follow the desired input angular distribution. It is also
possible to do an importance sampling on top of this to have the generation
focus on a given part of the phase space, although this is currently only
implemented in a rather rudimentary way (see [Importance
smapling](#importance-sampling))

After generating the chi\_c the decay to the J/psi + photon and subsequently the
decay of the J/psi into two oppositely charged muons is done by observing all
conservation laws (momentum, energy and angular momentum).

After all the particles have been generated they are smeared and depending on
the used settings the reconstruction efficiencies are evaluated using the
smeared muons and photons. Additionally it is possible to only select certain
events for storage on disk (see below for more details).

### Smearing

To emulate measurement uncertainties due to detector resolution or final state
radiation (FSR) the generated particles are smeared to obtain _reconstruction_
level particles. In smearing the particle momentum is smeared by a random
factor, whereas the mass remains unchanged. To emulate reconstruction effects
the generated muons and photons are smeared and then the _reconstruction_ level
J/psi and chi\_cJ particles are built _bottom-up_ from the smeared muons and
photon.

The muons are smeared using a Gaussian distribution whereas the photons are
smeared using a Crystal Ball distribution. The parameters of these distributions
are currently tuned to give a (rough) agreement with the J/psi (for the muons)
respectively chi\_c1 mass (mainly the photons) distributions as fitted from data
(8 TeV, 2012). Studies show that the exact details of these parameters do not
actually matter in the end.

### Reconstruction efficiencies

To emulate the reconstruction efficiencies each event gets assigned three
efficiency weights; The photon conversion and reconstruction efficiency
(AN-2015/11) as well as the single muon efficiencies (AN-2014/132) are used for
this. For each event the smeared particles are used to evaluate these
efficiencies depending on p\_T and |eta| of the corresponding particle. Each
weight is stored in a separate branch and the total reconstruction efficiency is
the product of the three efficiencies.

If an event is outside of the range where the efficiencies are defined than `-1`
is stored in the branch in order to easily distinguish these events from the
others. If the single muon and photon selections are applied than all the events
that are stored will have efficiencies > 0.

The single muon efficiencies have values between 0 and 1, whereas the photon
efficiencies return their value in %, so that is necessary to multiply them by a
factor of `0.01` to have a total efficiency that can actually be interpreted as
probability for an event to be reconstructed. However, for polarization
measurements the absolute value is not really important, so that this can also
be omitted.

To get a _reconstruction_ level distribution of any given variable the following
(exemplary) `TTree::Draw` command can be used (here the `JpsiPt` is plotted,
`tr` is the name of the `TTree` in the generated `.root` file). See below for a
list of [available branches](#available-branches). **NOTE:** This draw command
does not yet apply any selection and just illustrates how the total weight of
each event is determined.

```c++
tr->Draw("JpsiPt", "gamma_eff_sm * 0.01 * lepP_eff_sm * lepN_eff_sm");
```

### Importance sampling

To increase the efficiency of the generation in phase space regions that are
plagued by low acceptance a rudimentary importance sampling has been
implemented. At the moment only the costheta variable (depending on the J/psi
p\_T) is considered in this importance sampling.

If importance sampling is enabled the sampling kernel is evaluated using the
costheta value of the current event to obtain a value `w_sampling`. This value
is multiplied with the value obtained from evaluating the angular distribution W
for the current event. This product is then compared to a random number `r`
drawn from a uniform distribution between 0 and `r_max` and if it larger than
this number the event is kept, otherwise a new event is generated.

This sampling has to be considered when looking at physical distributions by
weighting each event by `1 / w_sampling`.

This means that to get a _recontruction_ level distribution of any given
variable that has been generated with importance sampling the following
`TTree::Draw` command has to be used (see the
[efficiency](#reconstruction-efficiencies) section for some more details):

```c++
tr->Draw("JpsiPt", "gamma_eff_sm * 0.01 * lepP_eff_sm * lepN_eff_sm / w_sampling");
```

## Code organization

The code is organized as follows: The main file is the
[chicpolgen.C](chicpolgen.C) file that includes some additional files from the
[interface](interface) directory. These included files are mainly helper files
that concern reading in efficiencies from `.root` files and later on providing
an easy interface to evaluating these efficiencies as a function of p\_T and eta
(see the [efficiencies.h](interface/efficiencies.h) file for this). A second
aspect they serve is to have an easily (at runtime) configurable way of
selecting certain events depending on the J/psi p\_T and rapidity as well as on
the single muon and photon p\_T and rapidity (see the
[select.h](interface/select.h) file for how this is implemented). A third header
file deals with providing an easy-to-use interface for smearing particles. This
is all done in [smearing.h](interface/smearing.h). This header file currently
still contains the now unused `SmearingProvider` class. Currently only the
`smearParticleGaus` and `smearParticleTF1` functions are used for smearing.

The [run_chicpolgen.cc](run_chicpolgen.cc) file is only a wrapper around
`chicpolgen.C` that does nothing else then parsing the command line arguments,
packing them into several configuration `struct`s and then calling `chicpolgen`
with them.

## Building

In order to run `run_chicpolgen` it has to be built first. This should be as
simple as calling `make` in this directory.

## Running

To run the Toy MC generation it is enough to run `./run_chicpolgen` (without any
arguments) after building it. By default this will produce an output file
`chicpolgen.root` containing 3 million unpolarized chi_c1 events with a p\_T
between 7 and 30 GeV and |y| < 1.3 assuming the HX frame is the natural
polarization frame. For a full list of the available command line arguments see
the [runtime configuration](#runtime-configuration-of-the-generation).

### Examples of commonly used generation settings

Since the full list of available configuration flags is quite large here are
some useful combinations of commonly used arguments

#### Specifying a polarization scenario

The desired J=(0,1,2) value is specified with the `--state` argument. The
desired polarization can be defined using the `--helicity1` and `--helicity2`
arguments. The value passed with `--helicity1` corresponds to `gen_config::R`
and the value passed with `--helicity2` corresponds to `gen_config::R2` in
`chicpolgen.C`.

The physical meaning of these parameters correspond to the fractions of the |Jm|
with m=(1,2) eigenstates of the total chi state in PRD **83** 096001 (2011).
Specifically `R` corresponds to to the sum of |b\_+1|^2 + |b\_-1|^2 in Eq. 21
and Eq. 22 and `R2` corresponds to the sum of |b\_+2|^2 + |b\_-2|^2 in Eq. 22 in
PRD **83** 096001 (2011). These equations can also be used to get the
corresponding values of lambda_theta for a given value of J.

This brings the following implications:
- for J=1, only `R` will be used (and has a physical meaning).
- for J=2, both `R` and `R2` will be used
- In both cases the sum of all amplitudes has to be 1. The |b_0|^2 is
  automatically determined by using this condition. Thus, it has to be ensured
  that the values that are used as inputs make sense physically as there are no
  (or only very limited) checks for this in the Toy MC generation. Specifically
  it has to be made sure that the sum of `R` and `R2` does not exceed 1.

The following is a list of settings to obtain polarization scenarios for which
there is only a polar anisotropy, since they are pure Eigenstates of Jz.

- unpolarized chi\_c1 (lambda_theta = 0): `--helicity1 0.6666666666667`
- maximally transversely polarized chi\_c1 (lambda\_theta = 1): `--helicity1 0`
- maximally longitudinally polarized chi\_c1 (lambda\_theta = -1/3): `--helicity1 1`
- unpolarized chi\_c2 (lambda_theta = 0): `--helicity1 0.4 --helicity2 0.4`
- maximally transversely polarized chi\_c2 (lambda\_theta = 1): `--helicity1 0
  --helicity2 1`
- maximally longitudinally polarized chi\_c2 (lambda\_theta = -3/5):
  `--helicity1 0 --helicity2 0`


#### Correction map histograms without storing any branches

To store only the histograms that are necessary for determining correction maps
the following settings can be used. These are (apart from specifying an output
file, the number of events to generate and enabling the importance sampling) the
same settings that have been used to obtain the correction maps for the narrow
pT range between 8 and 9 GeV.

```bash
./run_chicpolgen --storeHists true --storeBranches none \
     --jpsiSel true --psiPtMin 8 --psiPtMax 9 --psiRapMin 0.5 --psiRapMax 1.0 \
     --muonSel true --muonEffs EffFiles/single_muon_effs_noTracking_L3ptg2_final_fit.root \
     --photonSel true --photonEffs EffFiles/photon_effs_param.root
```

#### Generating a test sample with only the necessary branches

For testing different correction maps it is not necessary to store all branches.
A subset of branches suffices to be able to do a rather large set of tests. To
generate such a sample the following settings can be used. The settings below
already apply all the selections that are also applied on data and thus all
events in this sample can be used for testing. This makes it possible to omit
storing the single muon and photon kinematics and, thus, save disk space.

```bash
./run_chicpolgen --storeBranches JpsiPt JpsiRap costh_HX phi_HX gamma_eff_sm lepP_eff_sm lepP_eff_sm \
     --jpsiSel true --muonSel true --photonSel true \
     --muonEffs EffFiles/single_muon_effs_noTracking_L3ptg2_final_fit.root \
     --photonEffs EffFiles/photon_effs_param.root
```

## Runtime configuration of the generation

The settings of the generation can be changed via command line arguments. The
full list of arguments and a short explanation is provided here (each flag is
followed by the expected type and the default value for the flag):

- **`--genfile`** (`string`, `"chicpolgen.root"`): The name of the output file that will
  be generated.
- **`--helicity1`** (`float`, `2./3.`): The fraction of the |Jz|=1 eigenstate.
  Has to be between 0 and 1. (Only applicable for the J=(1,2) states).
- **`--helicity2`** (`float`, `0`): The fraction of the |Jz|=2 eigenstate. Has
  to be between 0 and 1. (Only applicable for the J=2 state).
- **`--state`** (`int`, `1`): The value of J in chi\_cJ. Can be either 0, 1 or 2 (other
  values are accepted, but the generation will go into an infinite loop)
- **`--nevents`** (`int`, `3000000`): The number of events to generate
- **`--ptmin`** (`float`, `7.0`): The minimum p\_T of the generated chi\_cJ 
- **`--ptmax`** (`float`, `30`): The maximum p\_T of the generated chi\_cJ 
- **`--ymin`** (`float`, `0`): The minimum **absolute** rapidity of the generated chi\_cJ
- **`--ymax`** (`float`, `0`): The maximum **absolute** rapidity of the generated chi\_cJ
- **`--CSframe`** (`bool`, `false`): Use the CS frame as natural frame instead of the HX
  frame.
- **`--naccept`** (`int`, `0`): The number of events that have to be accepted before the
  generation stops early. This can be useful if a certain amount of accepted
  events is necessary, so that is possible to set the `--nevents` flag to a
  large enough value and run the generation until enough events are accepted.
  Accepted events are counted after all possible selections (see below for some
  more details)
- **`--muonEffs`** (`string`, `""`): The file name where the muon efficiencies can be
  found. If an argument is passed here than the muon efficiencies will be
  evaluated for each event.
- **`--photonEffs`** (`string`, `""`): The file name where the photon efficiencies can
  be found. If an argument is passed here than the photon efficiencies will be
  evaluated for each event.
- **`--storeBranches`** (`list of strings`, `"all"`): List of branches that
  should be stored to the output file. Can be used to store only a subset of
  branches to the output file to save disk space by not storing unnecessary
  branches. The list of branches to store has to be passed as list of strings
  that are separated by whitespace. If `"all"` is passed all branches are
  stored, if `"none"` is passed, no branches are stored.
- **`--storeHists`** (`bool`, `false`): Enable the storing of costheta-phi histograms in
  the HX, PX and CS frames at three different stages during the generation:
  1. At _generation_ (`gen`) level, directly after the J/psi selection. These are named 
  2. At _acceptance_ (`acc`) level, directly after the single muon and photon
     selection (additionally to the J/psi selection)
  3. At _reconstruction_ (`reco`) level, where each event that passes the single
     muon and photon selection (additionally to the J/psi selection) is weighted
     by the product of the single muon and photon efficiencies when filling the
     histogram.
- **`--jpsiSel`** (`bool`, `false`): Apply the J/psi selection. All events that do
  not pass the J/psi selection will not be written to the output `TTree` (if it
  is filled).
- **`--psiPtMin`** (`float`, `8.0`): Minimum p\_T of the J/psi
- **`--psiPtMax`** (`float`, `20.0`): Maximum p\_T of the J/psi
- **`--psiRapMin`** (`float`, `0`): Minimum **absolute** rapidity of the J/psi
- **`--psiRapMax`** (`float`, `1.2`): Maximum **absolute** rapidity of of the J/psi
- **`--muonSel`** (`bool`, `false`): Apply the muon selection and store only those
  events for which both muons pass this selection. The muon selection is
  implemented in the `LooseMuonSelector` class [here](interface/select.h) and
  selects muons with the following kinematics:
  - p\_T > 3.5 GeV, if |eta| < 1.2
  - p\_T > 3.5 - (|eta| - 1.2) * 2.5 GeV, if 1.2 < |eta| < 1.6
- **`--photonSel`** (`bool`, `false`): Apply the photon selection and store only
  those events that pass it. The photon selection is implemented in terms of the
  `MinPtMaxEtaSelector` class [here](interface/select.h) and selects photons
  with the following kinematics:
  - p\_T > 0.41 GeV and |eta| < 1.5
- **`--sampling`** (`bool`, `false`): Use an importance sampling kernel for
  sampling the generation. In the current implementation this uses the
  `StepSamplingKernel` which is defined in [chicpolgen.C](chicpolgen.C).

### Available branches
The following branches are available currently and are stored if
`--storeBranches all` is used (some others need to be uncommented in
`chicpolgen.C` in order for them to be stored):

- **`gen_chicPt`**: The generated chi\_c p\_T
- **`gen_JpsiPt`**: The generated J/psi p\_T
- **`gen_chicRap`**: The generated chi\_c rapidity
- **`gen_JpsiRap`**: The generated J/psi rapidity
- **`gen_chicMass`**: The generated chi\_c mass
- **`gen_JpsiMass`**: The generated J/psi mass
- **`gen_photonPt`**: The generated photon p\_T
- **`gen_photonPl`**: The generated z-component of the photon momentum 
- **`gen_photonEta`**: The generated photon (psuedo-)rapidity
- **`gen_muPPt`**: The generated positively charged muon p\_T
- **`gen_muPEta`**: The generated positively charged muon pseudo-rapidity
- **`gen_muNPt`**: The generated negatively charged muon p\_T
- **`gen_muNEta`**: The generated negatively charged muon pseudo-rapidity
- **`cosTH_psi`**: The generated cosTheta of the J/psi in the chi\_c restframe
- **`costh_chihe`**: The costheta of the dimuon decay w.r.t. the J/psi direction
  seen in the chi\_c rest frame
- **`phi_chihe`**: The phi of the dimuon decay w.r.t. the J/psi direction seen
  in the chi\_c rest frame
- **`gen_costh_HX`**: The generated costheta (of the dimuon) in the HX frame
- **`gen_phi_HX`**: The generated phi (of the dimuon) in the HX frame
- **`gen_costh_CS`**: The generated costheta (of the dimuon) in the CS frame 
- **`gen_phi_CS`**: The generated phi (of the dimuon) in the CS frame
- **`chicPt`**: The chi\_c pT after smearing (_reconstructed_)
- **`chicRap`**: The chi\_c rapidity after smearing (_reconstructed_)
- **`mumugammaMass`**: The mumu+gamma invariant mass after smearing (_reconstructed_)
- **`chicMass`**: The Q-value (mumu+gamma mass - J/psi mass J/psi PDG mass)
  after smearing (_reconstructed_). This is basically the same as is used in
  real data.
- **`photonPt`**: The photon p\_T after smearing (_reconstructed_)
- **`photonEta`**: The photon (pseudo-)rapidity after smearing (_reconstructed_)
- **`JpsiPt`**: The J/psi p\_T after smearing (_reconstructed_) 
- **`JpsiRap`**: The J/psi rapidity after smearing (_reconstructed_) 
- **`JpsiMass`**: The J/psi mass after smearing (_reconstructed_) 
- **`muPPt`**: The positively charged muon p\_T after smearing (_reconstructed_)
- **`muPEta`**: The positively charged muon pseudo-rapidity after smearing
  (_reconstructed_)
- **`muNPt`**: The negatively charged muon p\_T after smearing (_reconstructed_) 
- **`muNEta`**: The negatively charged muon pseudo-rapidity after smearing
  (_reconstructed_)
- **`Q_value_gen`**: The Q-value at generation level 
- **`costh_HX`**: costheta after smearing in the HX frame (_reconstructed_) 
- **`phi_HX`**: phi after smearing in the HX frame (_reconstructed_) 
- **`costh_CS`**: costheta after smearing in the CS frame (_reconstructed_)   
- **`phi_CS`**: phi after smearing in the CS frame (_reconstructed_)  
- **`costh_PX`**: costheta after smearing in the PX frame (_reconstructed_)  
- **`phi_PX`**: phi after smearing in the PX frame (_reconstructed_) 
- **`cosTH_HX_sm`**: cosTheta after smearing in the HX frame (_reconstructed_)  
- **`cosTH_PX_sm`**: cosTheta after smearing in the HX frame (_reconstructed_)   
- **`cosTH_CS_sm`**: cosTheta after smearing in the HX frame (_reconstructed_)  
- **`w_sampling`**: The sampling weight (only present if importance sampling is
  enabled)
- **`lepP_eff_sm`**: The muon efficiency evaluated using the (smeared)
  positively charged muon as input
- **`lepN_eff_sm`**: The muon efficiency evaluated using the (smeared)
  positively charged muon as input
- **`gamma_eff_sm`**: The photon efficiency evaluated using the (smeared) photon
  as input


## Compile time configuration of the generation

At the beginning (after the includes) of `chicpolgen.C` there are three
preprocessor variables that allow to change the behavior of the generation at
compile time (i.e. changing these requires to re-compile `run_chicpolgen` via
`make` in order for the changes to be effective). Two of these slightly change
the physics of the generation while the third one is solely for instrumenting
the generation to obtain some timing measurements of the different steps of
generating one events. The variables are the following (with the default
values):
- **`GENRAPIDITY`** (`1`): If `1` then the chi\_cJ will be generated flat in
  rapidity. If `0` they will be generated flat in the pseudo-rapidity.
- **`GENPTM`** (`1`): If `1` then the p\_T/M will be drawn from a realistic
  p\_T/M distribution using the per event value of M and the p\_T value is
  obtained from there. If `0` then the PDG value of the specified state will be
  used while drawing the p\_T/M value. In practice, since the width (as reported
  by the PDG) for the chi states are all rather small the difference between the
  two settings is rather small.
- **`TIMING_INSTRUMENTATION`** (`0`): Allows to get some additional information
  about the generation process itself when set to `1`. The following branches
  become available if enabled: 
  - **`t_gen`**: The time it takes for an event to be generated (and accepted by
    the sampling). This includes the time it takes to generate the initial
    chi\_cJ and the subsequent decay down to the single muons and photons.
  - **`n_gen`**: The number of times an event as to be _re-generated_ before
    being accepted by the sampling. To get a better estimate of how much time it
    actually takes to generate and decay an event `t_gen` should be divided by
    this number.
  - **`t_smear`**: The time it takes for each event to go through the smearing
    (and _reconstruction_) machinery.
  - **`t_eff`**: The time it takes for each event to evaluate the efficiencies
    and (if applicable) also to fill the output histograms.

## Open points

- [ ] Currently building emits a warning from an included header (that is not
even strictly necessary at the moment). Although this can be safely ignored,
this should be fixed at some point.

## Current limitations

- [ ] The Toy MC currently is only able to generate an polar anisotropy in the
      natural frame. The phi direction will always be generated flat (in the
      natural frame).
