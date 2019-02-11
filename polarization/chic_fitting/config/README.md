# JSON config files for the ConfigFitModel
This README describes the (valid) structure of a JSON config file that can be
used as an input for the `ConfigFitModel` as it is defined in
[config_fitting.py](/python/utils/config_fitting.py).

## Necessary fields

The following fields are necessary to fully define a fit model.

### `fit_variable`
The name of the `RooRealVar` that is fitted

### `full_model`
An object describing the full model with the following fields:
- `name`: The name of the `RooAbsPdf` that describes the full model
- `type`: How the `sub_models` (see below) should be combined. **NOTE:** This will
be used directly in the `RooWorkspace::factory` expression.

### `sub_models`
A list of models that will be combined into the full model. Each model will be
multiplied by its `event_yield` variable before it is used in the combination to
the full model. Each of the defined models in this list will be plotted in the
`FitModel.plot` function.

Each of the models in the list is an object which has the following fields:
- `name`: The name of the `RooAbsPdf` that is used in the workspace
- `event_yield`: The number of events that will be obtained in the fit for this
model. This has to be a valid expression that defines this variable in the
`RooWorkspace` when it is used in the `RooWorkspace::factory` unless the
variable is already present (i.e. it is defined with some other expression; see
[`expression strings`](#expression_strings)). If the variable is not defined
directly here, than the name of the variable is enough.
- `expression`: The expression string describing this model, which is directly
passed to the `RooWorkspace::factory` after substituting the fit variable (see
above) and the name of the model (using `.format(fit_variable, name)` on the
expression string).
- `plot_config`: A list containing three values:
  - The line style that will be used in the fit plots for this model
  - The line color that will be used in the fit plots for this model
  - The legend entry that will be used in the fit plots for this model

### `plot_config`
The configuration steering the plotting. Currently has only one field:
- `legpos`: A list with 4 floats that fixes the position of the legend in the
plot. The numbers are directly passed to the `TLegend` constructor.

## Optional fields

The following fields can be optionally defined.

### `expression_strings`
Expression strings that are passed to the `RooWorkspace::factory` **before** the
`sub_models` and the `full_model` strings. Each string has to be a valid
expression and can be used to define variables and other conditions that can be
done via the `RooWorkspace::factory` interface. This can also be used to define
sub models that should not show up in the fit plots or that are used as
ingredients to other composite models that only show up in the fit plots after
they have been composed into one model.

### `fix_vars`
A list of dictionaries/objects where the keys/field names are the variables and
the values are the values to which the variable should be fixed before fitting.
Each of the variables that is fixed here has to be already present in the
`RooWorkspace`.

### `floating_costh`
A list of variable names that define the variables that are left to float in a
costh binned fit even after other parameters have been fixed to a costh
integrated fit. Each of the variables used here has to be already present in the
`RooWorkspace`. **NOTE:** The event yields of the `sub_models` will be left to
float automatically and do not have to go here. The only exception to this is
when the event yields depend on some other variable. Then the variable(s) that
are the ones that are floating in the fit have to be put here, otherwise the
yield variable will not be floating.
