# import classes from the esbmtk library
from esbmtk import (
    Model,  # the model class
    Reservoir,  # the reservoir class
)
# define the basic model parameters
M = Model(
    stop="3 Myr",  # end time of model
    max_timestep="1 kyr",  # upper limit of time step
    element=["Carbon"],  # list of element definitions
)

Reservoir(
    name="GO",  # Name of reservoir group
    geometry=[0, -10000, 1],  # upper, lower
    concentration={M.DIC: "0 mmol/kg"},
    register=M,
)

M.hyp.volume(0,-200)
M.hyp.area_dz(0,-250)
