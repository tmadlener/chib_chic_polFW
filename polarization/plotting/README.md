# Plotting and PPD scanning
This folder contains scripts for obtaining corrected ratios as well as plotting
them and also doing a PPD scan of the ratio and then plotting the result of
that.

*NOTE*: Currently only some of the scripts are described and also for the
described scripts only the major aspects are covered.

## `calc_ppd.py`

A short (and not yet complete) description of how to use `calc_ppd.py`. The
description currently serves only to document what the different histograms that
are created by the script (unless the `--no-hists` flag is passed).

Depending on which flags are used the histograms might change their meaning, so
in the following the different possible use cases are described.

### `--costh`
The scanning assumes that the passed `ratiofile` contains the fitted (and
corrected) ratio `r_chic2_chic1_v_costh_HX_fold_bin_0` and uses this to scan the
PPD assuming the costh ratio as the analytic function.

The costh ratio analytic function is:

```python
N * (1 + lth_2 * costh**2) / (1 + lth_1 * costh**2)
```

The flags for the `--max-costh` and `--costh-ratio` are ignored.

For the variables `lth` (lambda\_theta(chic1)), `lth2` (lambda\_theta(chic2)),
`dlth` (Delta\_lambdatheta), `norm_costh` (normalization) the following
histograms are produced:
- `prior_1d_X`: The prior distribution without taking into account the
  importance sampling of the normalization. These priors only make sense for the
  angular variables
- `prior_1d_X_norm`: The prior distribution taking into account the importance
  sampling of the normalization (See [Normalization importance
  sampling](#Normalization-importance-sampling)).
- `prior_1d_X_2d_const`: The prior distribution imposing the 2d constraints as
  defined by Eq. 27 and Eq. 28, and shown in Fig. 3 (PRD **83**, 096001 (2011))
  on the angular variables. The importance sampling of the normalization is not
  taken into account here.
- `prior_1d_X_norm_2d_const`: The prior distributions imposing the 2d
  constraints on the angular variables and also taking into account the
  importance sampling of the normalization.
- `ppd_1d_X`: The posterior probability density (PPD) of the scan given the
  input ratio assuming flat priors of `lth` and `lth2` in their allowed ranges.
- `ppd_1d_X_2d_const`: The PPD of the scan given the input ratio assuming priors
  that are flat in the 2d space defined by the 2d constraints on the angular
  distributions (see above).

### `--phi`
The scanning assumes that the passed `ratiofile` contains the fitted (and
corrected) ratio `r_chic2_chic1_v_phi_HX_fold_bin_0` and uses this to scan the
PPD assuming the phi ratio as analytic function.

The phi ratio analytic function is:
```python
N * (1 + kappa1 * cos(2 * phi)) / (1 + kappa2 * cos(2 * phi))
```

where kappa is defined as:

```python
kappa = (3 - max_costh**2) / (3 + lth * max_costh**2) * lph
```

This definition should make the usage of the `--max-costh` argument immediately
obvious in case of obtaining the PPD for a phi ratio. This value has to match
the one that has been used in the preselection of the events that go into the
phi fit.

For the phi ratio the list of produced histograms depends on whether the
`--costh-ratio` argument is used. If it is **not** used then for the variables
`lph` (lambda\_phi(chic1)), `lph2` (lambda\_phi(chic2)), `dlph`
(Delta\_lambdaphi) and `norm_phi` (normalization) the same histograms with the
same interpretations as for the [`--costh`](#--costh) flag are produced
(replacing `lth` and `lth2` with `lph` and `lph2` respectively where necessary).

#### with `--costh-ratio`
If a file containing a costh ratio is passed with the `--costh-ratio` argument
then instead this costh ratio is also considered in the fit and effectively
instead of a 1d fit a *1d + 1d* fit is done and information from both ratios is
leveraged. The 1d results can still be obtained as can be seen by the list of
histograms that are being produced. 

For the variables `lth` (lambda\_theta(chic1)), `lth2` (lambda\_theta(chic2)),
`dlth` (Delta\_lambdatheta), `norm_costh` (normalization along costh), `lph`
(lambda\_phi(chic1)), `lph2` (lambda\_phi(chic2)), `dlph` (Delta\_lambdaphi) and
`norm_phi` (normalization along phi), `ltilde` (lambda\_tilde(chic1)), `ltilde2`
(lambda\_tilde(chic2)) and `dltilde` (Delta\_lambdatilde) the following
histograms are produced:
- `prior_1d_X`: The prior distribution without taking into account the
  importance sampling of either normalization (costh or phi). These priors only
  make sense for the angular variables.
- `prior_1d_X_2d_const`: The prior distribution without taking into account the
  importance sampling, but imposing the 2d constraints.
- `prior_1d_X_norm`: The priors taking into account the importance sampling of
  both costh and phi.
- `prior_1d_X_norm_costh`: The priors only taking into account the importance
  sampling of the costh normalization. These are equivalent to the
  `prior_1d_X_norm` histograms for the [`--costh`](#--costh) only scan.
- `prior_1d_X_norm_phi`: The priors taking into account the importance sampling
  of the phi normalization. These are equivalent to the `prior_1d_X_norm`
  histograms for the `--phi` only scan.
- `ppd_1d_X`: The PPD of the scan given the costh and phi ratio as inputs
  assuming flat priors for `lth`, `lth2`, `lph` and `lph2` in their respective
  allowed ranges.
- `ppd_1d_X_2d_const`: The PPD of the scan given the costh and phi ratio as
  inputs assuming that the priors are flat in the 2d space allowed by the 2d
  constraints on the angular distributions.
- `ppd_1d_X_costh`: The PPD of the scan only using information from the costh
  ratio, assuming flat priors for `lth` and `lth2`. These are the same as
  `ppd_1d_X` in the `--costh` only scan.
- `ppd_1d_X_phi`: The PPD of the scan only using information from the costh
  ratio, assuming flat priors for `lph` and `lph2` as well as for `lth` and
  `lth2`. These are the same as the `ppd_1d_X` in the `--phi` only scan.
- `ppd_1d_X_(costh|phi)_2d_const`: Same as the two above, but imposing the 2d
  constraints on the priors.

### Normalization importance sampling
In order to increase the "efficiency" of the PPD scanning an importance sampling
is used for the normalization of the scanned ratio. The normalization values are
drawn from a normal distribution and the central value of this normal
distribution is determined from the input ratio. 

In case of a costh ratio the central value is simply put to the value of the
first bin in the ratio, since the analytic costh ratio assumes the value of the
normalization at costh = 0.

In case of a phi ratio a constant function is fitted to the ratio and the fitted
constant value is used as central value for the Gaussian importance sampling
kernel. The reason for this approach is that the phi ratio assumes the value of
the normalization for phi = 45 and this should give a pretty good approximation
of that even for non-flat phi ratios.

The two importance sampling kernels are completely independent of each other, so
if a costh and a phi ratio is passed they will not interfere with each other.
