"""Compare ESBMTK results with published data.

Sepcifically Boudreau et al. see https://doi.org/10.1029/2009gb003654

Authors: Uli Wortmann & Tina Tsan

Copyright (C), 2024 Ulrich G. Wortmann & Tina Tsan

     This program is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import boudreau_2010 as bd
import numpy as np
from esbmtk import (
    DataField,
    ExternalData,
    Signal,
    Source,
    Species2Species,
    carbonate_system_2_pp,
    gas_exchange_fluxes,
)

run_time = "3800 yr"
time_step = "1 month"
rain_ratio = 0.3
alpha = 0.6

# import the base model
M = bd.initialize_model(rain_ratio, alpha, run_time, time_step)

# Create a Signal Instance for the carbon pulse
M.CP = Signal(
    name="CP",  # name of signal
    species=M.CO2,  # species
    filename="IS92a-scenario.csv",
    scale=0.877,  # required to match the total carbon release
    register=M,
)

# Add a source for the carbon pulse
M.Carbon_Pulse = Source(name="Carbon_Pulse", species=M.CO2)
M.C_CP = Species2Species(
    source=M.Carbon_Pulse,
    sink=M.CO2_At,
    rate="0 mol/yr",
    signal=M.CP,  # list of processes
)

M.read_state(directory="steady_state")  # read pre-industrial steady state
M.run()
M.save_data()

# ------------------ Create the figures ----------------------- #

# get CaCO3_export in mol/year and post process the data
CaCO3_export = M.CaCO3_export.to(f"{M.f_unit}").magnitude
carbonate_system_2_pp(M.D_b, CaCO3_export, 200, 10999)

# get the air sea gas-exchange fluxes
GEX_L = gas_exchange_fluxes(M.L_b.DIC, M.CO2_At, "4.8 m/d")
GEX_H = gas_exchange_fluxes(M.H_b.DIC, M.CO2_At, "4.8 m/d")

# import the digitized data and create data fields to compare the
# results
data = "dic_l dic_h dic_d TA_l TA_h TA_d hplus_l hplus_h hplus_d Fburial"
data = data + " zsat zcc zsnow EL EH pco2 Cpulse"
for n in data.split(" "):
    ExternalData(
        name=f"ef_{n}",
        filename=f"digitized/{n}.csv",
        legend=n,
        register=M,
    )

DataField(
    name="df_pulse",
    x1_data=[M.time, M.ef_Cpulse.x],
    y1_data=[M.CP.signal_data.m, M.ef_Cpulse.y],
    y1_label=["This model", "Fig2A"],
    y1_color=["C0", "C1", ""],
    y1_style=["solid", "dotted"],
    y1_legend="C [mol/year]",
    title="g)",
    register=M,
)

DataField(
    name="df_GEX",
    x1_data=[M.time, M.time, M.ef_EL.x, M.ef_EH.x],
    y1_data=[GEX_L, GEX_H, M.ef_EL.y, M.ef_EH.y],
    y1_label=["Low-Lat gex", "High Lat gex", "d_LL", "d_HL"],
    y1_color=["C0", "C1", "C0", "C1", ""],
    y1_style=["solid", "solid", "dotted", "dotted"],
    y1_legend="Gas Exchange Flux [mol/yr]",
    title="d)",
    register=M,
)

DataField(
    name="df_DIC",
    x1_data=[M.time, M.time, M.time, M.ef_dic_l.x, M.ef_dic_h.x, M.ef_dic_d.x],
    y1_data=[
        M.L_b.DIC.c * 1000,
        M.H_b.DIC.c * 1000,
        M.D_b.DIC.c * 1000,
        M.ef_dic_l.y * 1000,
        M.ef_dic_h.y * 1000,
        M.ef_dic_d.y * 1000,
    ],
    y1_label=["Low latitude", "High latitude", "Deep box", "d_L", "d_H", "d_D"],
    y1_legend="DIC [mmol/L]",
    y1_color=["C0", "C1", "C2", "C0", "C1", "C2"],
    y1_style=["solid", "solid", "solid", "dotted", "dotted", "dotted"],
    title="a)",
    register=M,
)

DataField(
    name="df_TA",
    x1_data=[M.time, M.time, M.time, M.ef_TA_l.x, M.ef_TA_h.x, M.ef_TA_d.x],
    y1_data=[
        M.L_b.TA.c * 1000,
        M.H_b.TA.c * 1000,
        M.D_b.TA.c * 1000,
        M.ef_TA_l.y * 1000,
        M.ef_TA_h.y * 1000,
        M.ef_TA_d.y * 1000,
    ],
    y1_label=["Low latitude", "High latitude", "Deep box", "d_L", "d_H", "d_D"],
    y1_color=["C0", "C1", "C2", "C0", "C1", "C2"],
    y1_style=["solid", "solid", "solid", "dotted", "dotted", "dotted"],
    y1_legend="TA [mmol/L]",
    title="b)",
    register=M,
)

DataField(
    name="df_pH",
    x1_data=[M.time, M.time, M.time, M.ef_hplus_l.x, M.ef_hplus_h.x, M.ef_hplus_d.x],
    y1_data=[
        M.L_b.pH.c,
        M.H_b.pH.c,
        M.D_b.pH.c,
        -np.log10(M.ef_hplus_l.y),
        -np.log10(M.ef_hplus_h.y),
        -np.log10(M.ef_hplus_d.y),
    ],
    y1_label=["Low latitude", "High latitude", "Deep box", "d_L", "d_H", "d_D"],
    y1_color=["C0", "C1", "C2", "C0", "C1", "C2"],
    y1_style=["solid", "solid", "solid", "dotted", "dotted", "dotted"],
    y1_legend="pH",
    register=M,
    title="c)",
)

DataField(
    name="df_Depths",
    x1_data=[M.time, M.time, M.time, M.ef_zsat.x, M.ef_zcc.x, M.ef_zsnow.x],
    y1_data=[
        -M.D_b.zsat.c,
        -M.D_b.zcc.c,
        -M.D_b.zsnow.c,
        -M.ef_zsat.y,
        -M.ef_zcc.y,
        -M.ef_zsnow.y,
    ],
    y1_label=["zsat", "zcc", "zsnow", "d_zsat", "d_zcc", "d_zsnow"],
    y1_color=["C0", "C1", "C2", "C0", "C1", "C2"],
    y1_style=["solid", "solid", "solid", "dotted", "dotted", "dotted"],
    y1_legend="Depth (m)",
    title="e)",
    register=M,
)

DataField(
    name="df_burial",
    x1_data=[M.time, M.time, M.ef_Fburial.x, M.ef_Fburial.x],
    y1_data=[
        M.D_b.Fburial.c,  # burial from model
        M.D_b.Fdiss.c,  # dissolution from model
        M.ef_Fburial.y,  # burial as digitized
        CaCO3_export - M.ef_Fburial.y,  # dissolution from digitized burial
    ],
    y1_label=["Fburial", "Fdiss", "Fburial_d", "Fdiss_d"],
    y1_color=["C0", "C1", "C0", "C1"],
    y1_style=["solid", "solid", "dotted", "dotted"],
    y1_legend="C [mol/yr]",
    title="h)",
    register=M,
)

DataField(
    name="df_atm",
    x1_data=[M.time, M.ef_pco2.x],
    y1_data=[M.CO2_At.c * 1e6, M.ef_pco2.y],
    y1_color=["C0", "C1"],
    y1_style=["solid", "dotted"],
    y1_label=["pCO2", "d_pCO2"],
    y1_legend="ppm",
    title="f)",
    register=M,
)

DataField(
    name="df_Carbon_pulse",
    y1_data=M.CP.signal_data.m,
    y1_label="CO2 input",
    y1_legend="C mol/yr",
    register=M,
)

M.plot(
    [
        M.df_DIC,
        M.df_TA,
        M.df_pH,
        M.df_GEX,
        M.df_Depths,
        M.df_atm,
        M.df_pulse,
        M.df_burial,
    ],
    fn="comparison_w_d.pdf",
    title="ESBMTK vs Boudreau 2010 et al.",
)
