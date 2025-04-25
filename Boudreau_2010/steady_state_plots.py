"""Calculate the steady state concentrations for the Boudreau Model.

Plot the results against published values.
"""

import importlib.metadata

import boudreau_2010 as bd
import numpy as np
from esbmtk import carbonate_system_2_pp, data_summaries
from packaging import version

current_version = importlib.metadata.version("esbmtk")
min_required_version = "0.14.2"
if version.parse(current_version) < version.parse(min_required_version):
    raise ValueError("Please update esbmtk to version 0.14.2 or higher")

run_time = "1000 kyr"
time_step = "100 yr"  # this is max timestep
rain_ratio = 0.3
alpha = 0.6

# import the model definition
M = bd.initialize_model(rain_ratio, alpha, run_time, time_step)
M.run()
M.save_state(directory="steady_state")
# --------- Create figures and print diagnostic data ---------------#
""" Some of the tracer and critical depth interval data are not
part of the ODE model. Rather, we have to calculate them in a
post-processing step. CS2 requires the CaCO3 export flux, so we
query the model for the number used earlier, and then feed this
to the carbonate_system_2_pp which will compute the missing tracers.
"""
CaCO3_export = M.CaCO3_export.to(f"{M.f_unit}").magnitude
carbonate_system_2_pp(M.D_b, CaCO3_export, 200, 10999)
M.save_data()

""" Create plots showing the results of the model run.
ESBMTk provides the function data_summaries than be used to
create a list of ESBMTK objects that is then used as input to
the ESBMTK plot() method. This method will plot the respective
ESBTMK objects into a common figure.
"""

species_names = [M.DIC, M.TA, M.pH, M.CO3, M.zcc, M.zsat, M.zsnow]
box_names = [M.L_b, M.H_b, M.D_b]
pl = data_summaries(M, species_names, box_names, M.L_b.DIC)
pl += [M.CO2_At]
# plot the model results, but do not render the plot
plt, fig, axs = M.plot(
    pl,
    fn="steady_state.pdf",
    title="ESBMTK Preindustrial Steady State",
    no_show=True,
)

# add comparisons with published data
data = [1952e-6, 2153e-6, 2291e-6, 2288e-6, 2345e-6, 2399e-6]

# Add digitized data to model results
# FIXME: This seems to have no effect

v = 0
u = 0

axs = np.array(axs).reshape(4, 2)
for i in range(2):
    for j in range(3):
        d = data[u]
        axs[j, i].scatter(10, d, color=f"C{j}")
        u = u + 1
        axs[j, i].autoscale(enable=None, axis="y")

# CO32- values after Boudreau et al. 2010, Tab 3
axs[1, 1].scatter(1, 234e-6, color="C0")
axs[1, 1].scatter(1, 138e-6, color="C1")
axs[1, 1].scatter(1, 86e-6, color="C2")

# # zcc, zsat, zsnow values after Boudreau et al. 2010, Tab 3
axs[2, 0].scatter(1, 4750, color="C0")  # zcc
axs[2, 0].set_ylim([4700, 4900])
axs[2, 1].scatter(1, 3715, color="C0")  # zsat
axs[2, 1].set_ylim([3600, 3800])
axs[2, 1].set_title("zsat")  # not sure why this is needed
axs[3, 0].scatter(1, 4750, color="C0")  # zsnow
axs[3, 0].set_ylim([4700, 4900])

fig.tight_layout()
plt.show(block=False)
fig.savefig("steady_state.pdf")

# Printout the final concentrations
m = [M.L_b.DIC, M.L_b.TA, M.H_b.DIC, M.H_b.TA, M.D_b.DIC, M.D_b.TA]
for i, _n in enumerate(data):
    d = data[i]
    print(f"Delta {m[i].full_name} = {(m[i].c[0] - d) * 1e6:.2f} [umol/kg]")

# Print a table with the equilibrium constants for each box
ks = "K0, K1, K2, KW, KB"
kl = ks.split(", ")
rl = [M.L_b, M.H_b, M.D_b]

print(r"\nBox, {ks}, DIC[\mu{{}}mol/kg], TA[\mu{{}}mol/kg]")
for b in rl:
    print(f"{b.name}, ", end="")
    for k in kl:
        v = getattr(b.swc, k)
        print(f"{v:.4e}, ", end="")
    dic = b.DIC.c[-1] * 1e6
    ta = b.TA.c[-1] * 1e6
    print(f"{dic:.0f}, {ta:.0f}")
