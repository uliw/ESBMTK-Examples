# ESBMTK Examples

This repo holds example code for the [Earth Science Box Modeling Toolkit.](https://github.com/uliw/esbmtk) Please see the <https://esbmtk.readthedocs.io/en/latest/> for ESBMTK installation and usage instruction.

To run these examples, you will need at least ESBMTK 0.14.2.x

Download the examples via:

```sh
git clone git@github.com:uliw/ESBMTK-Examples.git
```

or as zipfile via <https://github.com/uliw/ESBMTK-Examples/archive/refs/heads/main.zip>


## Examples from the Manual

This directory contains the example code that is shown in the manual. For details and explanations see <https://esbmtk.readthedocs.io/en/latest/>


## Boudreau 2010 et al.

This directory contains the ESBMTK implementation of the Boudreau et al. 2010 model see <https://doi.org/10.1029/2009gb003654>

-   `model_walkthrough.ipynb` is Jupyter Notebook explaining the model setup.
-   `boudreau2010.py` is the model definition
-   `steady_state_plots.py` runs the model as a steady state problem
-   `is92a_comparison_plots.py` runs the model with the is92a scenario.


## Hypsography

This directory contains a script to calculate the hypsometric data used by ESBMTK. The data is being retrieved [PyGMT](https://www.pygmt.org/), which needs to be installed manually via conda or pip.

-   `hypsometry.py`
