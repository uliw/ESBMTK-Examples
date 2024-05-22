""" The ESBMTK implmentation of the Boudreau et al. 2010
model see https://doi.org/10.1029/2009gb003654

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


def initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step):
    from esbmtk import (
        Model,
        Q_,
        GasReservoir,
        create_bulk_connections,
        create_reservoirs,
        build_ct_dict,
        add_carbonate_system_1,
        add_carbonate_system_2,
        Species2Species,
        ConnectionProperties,
    )

    M = Model(
        stop=run_time,  # end time of model
        timestep=time_step,  # time step
        element=[  # list of elements we consider in the model
            "Carbon",
            "Boron",
            "Hydrogen",
            "Phosphor",
            "Oxygen",
            "misc_variables",
        ],
        mass_unit="mol",
        volume_unit="l",
        concentration_unit="mol/kg",
        opt_k_carbonic=13,  # Use Millero 2006
        opt_pH_scale=3,  # 1:total, 3:free scale
        opt_buffers_mode=2,  # carbonate, borate water alkalinity only
    )

    # -------------------- Set up box parameters ------------------------ #
    """ Boudreau et al defined their reservoirs by area and and volume, and
    explicitly assign temperature, pressure and salinity.
    """
    box_parameter: dict = {  # name: [[geometry], T, P]
        "H_b": {  # High-Lat Box
            "g": {"area": "0.5e14m**2", "volume": "1.76e16 m**3"},  # geometry
            "T": 2,  # temperature in C
            "P": 17.6,  # pressure in bar
            "S": 35,  # salinity in psu
        },
        "L_b": {  # Low-Lat Box
            "g": {"area": "2.85e14m**2", "volume": "2.85e16 m**3"},
            "T": 21.5,
            "P": 5,
            "S": 35,
        },
        "D_b": {  # Deep Box
            "g": {"area": "3.36e14m**2", "volume": "1.29e18 m**3"},
            "T": 2,
            "P": 240,
            "S": 35,
        },  # sources an sinks
        "Fw": {"ty": "Source", "sp": [M.DIC, M.TA]},
        "Fb": {"ty": "Sink", "sp": [M.DIC, M.TA]},
    }

    """ Boudreau at al 2010, provide volume and initial concentrations
    for of DIC and TA for each box. Here we use a shortcut and calculate
    the average concentration of DIC and TA and start the model
    from a state where DIC and TA concentrations are uniform.
    If you need box specific initial conditions use the output of
    build_concentration_dicts as starting point, but it is usually
    faster to just run the model to steady state and use this as a
    starting condition.
    """
    V_t = 1.76e16 + 2.85e16 + 1.29e18
    DIC_0 = (2.153 * 1.76e16 + 1.952 * 2.85e16 + 2.291 * 1.29e18) / V_t
    TA_0 = (2.345 * 1.76e16 + 2.288 * 2.85e16 + 2.399 * 1.29e18) / V_t

    initial_conditions: dict = {
        # species: concentration, Isotopes, delta value
        M.DIC: [Q_(f"{DIC_0} mmol/kg"), False, 0],
        M.TA: [Q_(f"{TA_0} mmol/kg"), False, 0],
    }

    # initialize reservoirs
    create_reservoirs(box_parameter, initial_conditions, M)

    # ------------------------------ Transport Processes ----------------- #
    # get a list of all species
    species_list: list = list(initial_conditions.keys())

    """define the mixing between high latitude box and deep water
    through a dictionary that specifies the respective source and sink
    reservoirs, connection id,  the connection type, the scaling factor
    and the list of species that will be affected.
    """
    connection_dict = {
        "H_b_to_D_b@mix_down": {  # source_to_sink@id
            "ty": "scale_with_concentration",  # type
            "sc": Q_("30 Sverdrup"),  # * M.H_b.swc.density / 1000,  # scale
            "sp": species_list,
        },
        "D_b_to_H_b@mix_up": {
            "ty": "scale_with_concentration",
            "sc": Q_("30 Sverdrup"),  # * M.D_b.swc.density / 1000,
            "sp": species_list,
        },
    }
    create_bulk_connections(connection_dict, M, mt="1:1")

    """ Specify the upwelling connnections. In order to save some typing
    we use a helper function to create the connection_dictionary by
    first creating a dictionary that contains the respective scaling
    factors, and then use the build_ct_dict() function to merge this
    data into a connection_dictionary.
    """
    conveyor_belt_transport = {  # advection
        "L_b_to_H_b@thc": Q_("25 Sverdrup"),  # * M.L_b.swc.density / 1000,
        # Downwelling
        "H_b_to_D_b@thc": Q_("25 Sverdrup"),  # * M.H_b.swc.density / 1000,
        # Upwelling
        "D_b_to_L_b@thc": Q_("25 Sverdrup"),  # * M.D_b.swc.density / 1000,
    }
    connection_dict: dict = build_ct_dict(
        conveyor_belt_transport,  # dict with scaling factors
        {
            "ty": "scale_with_concentration",  # connection type
            "sp": species_list,  # affected species
        },
    )
    create_bulk_connections(connection_dict, M, mt="1:1")

    """ Organic matter flux P - 200 Tm/yr
    POM = Particulate Organic matter. Since this model uses a fixed rate we
    can declare this flux with the rate keyword. Boudreau 2010 only considers
    the effect on DIC, and ignoresthe effect of POM on TA.

    Note the "bp" keyword: This specifies that the connection will remove the respective
    species form the source reservoir, but will bypass the addition to the sink
    reservoir. It is used here for the CaCO3 export, since carbonate remineralization
    his handled by the carbonate_system_2(). CS2 will calculate the amount of CaCO3
    that is dissolved and add this to the deep-box. The amount that buried is thus implicit.
    This is equivalent to export all CaCO3 into the deeb-box (i.e. bp="None"),
    and then creating an explicit connection that describes burial into the sediment.
    Since these a fixed rates, they could also be combined into one flux for DIC and
    one flux for TA.
    """
    M.OM_export = Q_("200 Tmol/a")
    M.CaCO3_export = Q_("60 Tmol/a")

    # Fluxes going into deep box
    connection_dict = {
        "L_b_to_D_b@POM": {  # DIC from organic matter
            "sp": M.DIC,
            "ty": "Regular",
            "ra": M.OM_export,
        },
        "L_b_to_D_b@PIC_DIC": {  # DIC from CaCO3
            "sp": M.DIC,
            "ty": "Regular",
            "ra": M.CaCO3_export,
            "bp": "sink",
        },
        "L_b_to_D_b@PIC_TA": {  # TA from CaCO3
            "sp": M.TA,
            "ty": "Regular",
            "ra": M.CaCO3_export * 2,
            "bp": "sink",
        },
    }
    create_bulk_connections(connection_dict, M, mt="1:1")

    # -------------------- Carbonate System- virtual reservoir definitions -- #
    """ To setup carbonate chemistry, one needs to know the soource adn sink of the
    export production fluxes. We thus keep two lists one for the surface boxes
    and one for the deep boxes.
    
    Carbonate system ` calculates tracer like CO2aq,  CO3, and H+, for the
    surface boxes.

    Carbonate System ` calculates the above tracers, and additionally
    the critical depth intervals (Saturation, CCD, and snowline), as well
    as the amount of carbonate that it dissolved.
    """
    surface_boxes: list = [M.L_b, M.H_b]
    deep_boxes: list = [M.D_b]
    ef = M.flux_summary(filter_by="PIC_DIC", return_list=True)

    add_carbonate_system_1(surface_boxes)

    add_carbonate_system_2(
        r_db=deep_boxes,  # list of reservoir groups
        r_sb=surface_boxes,  # list of reservoir groups
        carbonate_export_fluxes=ef,  # list of export fluxes
        z0=-200,  # depth of shelf
        alpha=alpha,  # dissolution coefficient
    )
    # -------------------- Atmosphere -------------------------
    GasReservoir(
        name="CO2_At",
        species=M.CO2,
        species_ppm="280 ppm",
        register=M,
    )

    """ GasExchange connections currently do not support the setup
    with the ConnectionsProperties class, since they connect CO2 to
    DIC which fools the automatic species matching logic. As such
    we use the Species2Species class to create the connection
    explicitly.
    """
    pv = "4.8 m/d"  # piston velocity
    Species2Species(
        source=M.CO2_At,  # Reservoir Species
        sink=M.H_b.DIC,  # Reservoir Species
        species=M.CO2,
        ref_species=M.H_b.CO2aq,
        solubility=M.H_b.swc.SA_co2,  # float
        area=M.H_b.area,
        id="H_b",
        piston_velocity=pv,
        water_vapor_pressure=M.H_b.swc.p_H2O,
        register=M,
        ctype="gasexchange",
    )

    Species2Species(
        source=M.CO2_At,  # Reservoir Species
        sink=M.L_b.DIC,  # Reservoir Species
        species=M.CO2,
        ref_species=M.L_b.CO2aq,
        solubility=M.L_b.swc.SA_co2,  # float
        area=M.L_b.area,
        id="L_b",
        piston_velocity=pv,
        water_vapor_pressure=M.L_b.swc.p_H2O,
        register=M,
        ctype="gasexchange",
    )

    # create the weathering fluxes
    ConnectionProperties(
        source=M.Fw,
        sink=M.L_b,
        rate={M.DIC: "12 Tmol/a", M.TA: "24 Tmol/a"},
        id="weathering",
        species=[M.DIC, M.TA],
        ctype="regular",
    )

    return M


if __name__ == "__main__":
    from math import log10
    from esbmtk import data_summaries, ExternalData
    from esbmtk import carbonate_system_2_pp

    run_time = "1000 kyr"
    time_step = "100 yr"  # this is max timestep
    rain_ratio = 0.3
    alpha = 0.6

    M = initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step)
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

    # import the digitized data and create data fields to compare the
    # results. Names must match with filenames
    data = "dic_l dic_h dic_d TA_l TA_h TA_d hplus_l hplus_h hplus_d".split(" ")
    # data = data + " zsat zcc zsnow EL EH pco2 Cpulse"
    for n in data:
        ExternalData(
            name=f"ef_{n}",
            filename=f"digitized/{n}.csv",
            legend=n,
            register=M,
        )

    # Add digitized data to model results
    v = 0
    for i in range(2):
        for j in range(3):
            d = getattr(M, f"ef_{data[v]}")
            if "hplus" in data[v]:
                axs[0, i].scatter(1, -log10(d.y[0]), color=f"C{j}")
            else:
                axs[0, i].scatter(1, d.y[0], color=f"C{j}")
            v = v + 1

    # CO32- values after Boudreau et al. 2010, Tab 3
    axs[1, 1].scatter(1, 234e-6, color="C0")
    axs[1, 1].scatter(1, 138e-6, color="C1")
    axs[1, 1].scatter(1, 86e-6, color="C2")
    # # zcc, zsat, zsnow values after Boudreau et al. 2010, Tab 3
    axs[2, 0].scatter(1, 4750, color="C0")  # zcc
    axs[2, 0].set_ylim([4700, 4900])  
    axs[2, 1].scatter(1, 3715, color="C0") # zsat
    axs[2, 1].set_ylim([3600, 3800])
    axs[2, 1].set_title("zsat") # not sure why this is needed
    axs[3, 0].scatter(1, 4750, color="C0") # zsnow
    fig.tight_layout()
    plt.show(block=False)
    fig.savefig("steady_state.pdf")
