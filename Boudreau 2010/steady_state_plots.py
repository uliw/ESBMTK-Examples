from math import log10
from esbmtk import data_summaries, ExternalData
from esbmtk import carbonate_system_2_pp
import boudreau_2010 as bd

run_time = "1000 kyr"
time_step = "100 yr"  # this is max timestep
rain_ratio = 0.3
alpha = 0.6

M = bd.initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step)
M.read_state(directory="init_data_3")  # get steady state
M.run()
M.save_state(directory="init_data_3")

""" Some of the tracers and critical depth interval data are not
part of the ODE model. Rather, we have to calculate them in a
post-processing step. CS2 requires the CaCO3 export flux, so we
query the model for the number used earlier, and then feed this
to the carbonate_system_2_pp which will the missing tracers.
"""
CaCO3_export = M.CaCO3_export.to(f"{M.f_unit}").magnitude
carbonate_system_2_pp(M.D_b, CaCO3_export, 200, 10999)
M.save_data()

""" Create some plots shiwing the results of the model run.
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
# data = "dic_l dic_h dic_d TA_l TA_h TA_d hplus_l hplus_h hplus_d".split(" ")
data = [1952e-6, 2153e-6, 2291e-6, 2288e-6, 2345e-6, 2399e-6]
# Add digitized data to model results
v = 0
u = 0
for i in range(2):
    for j in range(3):
        d = data[u]
        axs[0, i].scatter(10, d, color=f"C{j}")
        u = u + 1
    axs[v, i].autoscale(enable=None, axis="y")
    axs[v, i].autoscale()

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

m = [M.L_b.DIC, M.L_b.TA, M.H_b.DIC, M.H_b.TA, M.D_b.DIC, M.D_b.TA]
data = [1952e-6, 2288e-6, 2153e-6, 2345e-6, 2291e-6, 2399e-6]
for i, n in enumerate(data):
    d = data[i]
    print(f"Delta {m[i].full_name} = {(m[i].c[0] - d)*1e6:.2f} [umol/kg]")

# Print a table with the equilibrium constants for each box
ks = "K0, K1, K2, KW, KB"
kl = ks.split(", ")
rl = [M.L_b, M.H_b, M.D_b]
print(f"\nBox, {ks}, DIC[\mu{{}}mol/kg], TA[\mu{{}}mol/kg]")
for b in rl:
    print(f"{b.name}, ", end="")
    for k in kl:
        v = getattr(b.swc, k)
        print(f"{v:.4e}, ", end="")
    dic = b.DIC.c[-1] * 1e6
    ta = b.TA.c[-1] * 1e6
    print(f"{dic:.0f}, {ta:.0f}")
