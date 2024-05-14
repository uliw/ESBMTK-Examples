import numpy as np


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
        Connect,
        ConnectionProperties,
    )

    M = Model(
        stop=run_time,  # end time of model
        timestep=time_step,  # time step
        element=[
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
        bio_pump_functions=0,
        opt_k_carbonic=13,  # Use Millero 2006
        opt_pH_scale=3,  # 1:total, 3:free scale
        opt_buffers_mode=2,  # carbonate, borate water alkalinity only
        # parse_model=False,
    )

    # -------------------- Set up box parameters ------------------------ #
    """ boxes are defined by area and depth interval here we use a 
    dictionary to define the box geometries. The next column is
    temperature in deg C, followed by pressure in bar the geometry is
    [upper depth datum, lower depth datum, area percentage], or,
    as in this case, we give area and volume explicitly.
    """

    box_parameter: dict = {  # name: [[geometry], T, P]
        "H_b": {  # High-Lat Box
            "g": {"area": "0.5e14m**2", "volume": "1.76e16 m**3"},
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
    starting condition
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
    # species_list: list = [M.DIC, M.TA]

    """ define the mixing between high latitude box and deep water
    through a dictionary that specifies the respective source and sink
    reservoirs, connection id,  the connection type, the scaling factor
    and the list of species that will be affected.
    """
    connection_dict = {
        "H_b_to_D_b@mix_down": {  # source_to_sink@id
            "ty": "scale_with_concentration",  # type
            "sc": Q_("30 Sverdrup") * M.H_b.swc.density / 1000,  # scale
            "sp": species_list,
        },
        "D_b_to_H_b@mix_up": {
            "ty": "scale_with_concentration",
            "sc": Q_("30 Sverdrup") * M.D_b.swc.density / 1000,
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
        "L_b_to_H_b@thc": Q_("25 Sverdrup") * M.L_b.swc.density / 1000,
        # Downwelling
        "H_b_to_D_b@thc": Q_("25 Sverdrup") * M.H_b.swc.density / 1000,
        # Upwelling
        "D_b_to_L_b@thc": Q_("25 Sverdrup") * M.D_b.swc.density / 1000,
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
    POM = Particulate Organic matter. Unlike mix, this is really just a name to keep
    things organized. Since this model uses a fixed rate we can declare this flux
    with the rate keyword. Boudreau 2010 only considers the effect on DIC, and ignores
    the effect of POM on TA.

    Note the "bp" keyword: This specifies that the connection will remove the respective
    species form the source reservoir, but will bypass the addition to the sink
    reservoir. It is used here for the CaCO3 export, since carbonate remineralization
    his handled by the carbonate_system_2(). CS2 will calculate the amount of CaCO3
    that is dissolved and add this to the deep-box. The amount that buried is thus implicit.
    This is equivalent to export all CaCO3 into the deeb-box (i.e. bp="None"),
    and then creating an explicit connection that describes burial into the sediment.
    """
    M.OM_export = Q_("200 Tmol/a")
    # M.CaCO3_export = M.OM_export * rain_ratio
    M.CaCO3_export = Q_("60 Tmol/a")

    # Fluxes going into deep box
    connection_dict = {
        "L_b_to_D_b@POM": {  # organic matter
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
    """ To setup carbonate chemistry, one needs to know the export production,
    and where it ends up. We thus keep two lists one for the surface boxes
    and one for the deep boxes.
    
    Carbonate system one calculates tracer like CO2aq,  CO3, and H+, for the
    surface boxes.

    Carbonate System two calculates the aabove tracers, and additionally
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
        reservoir_mass="1.833E20 mol",
        species_ppm="280 ppm",
        register=M,
    )

    pv = "4.8 m/d"  # piston velocity
    Connect(
        source=M.CO2_At,  # Reservoir
        sink=M.H_b.DIC,  # ReservoirGroup
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

    Connect(
        source=M.CO2_At,  # Reservoir
        sink=M.L_b.DIC,  # ReservoirGroup
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

    # Connect(  # CaCO3 weathering
    #     source=M.Fw.DIC,  # source of flux
    #     sink=M.L_b.DIC,
    #     reservoir_ref=M.CO2_At,
    #     ctype="weathering",
    #     id="wca",
    #     scale=1,
    #     ex=0.2,
    #     pco2_0="280 ppm",
    #     rate="12 Tmol/a",
    #     register=M,
    # )

    # Connect(  #
    #     source=M.Fw.TA,  # source of flux
    #     sink=M.L_b.TA,  # target of flux
    #     ctype="scale_with_flux",
    #     ref_flux=M.flux_summary(filter_by="wca", return_list=True)[0],
    #     scale=2,
    #     id="wca_ta",
    # )

    # Connect(
    #     source=M.Fw.DIC,
    #     sink=M.L_b.DIC,
    #     rate="12 Tmol/a",
    #     register=M,
    # )
    # ConnectionProperties(
    #     source=M.Fw,
    #     sink=M.L_b,
    #     rate="12 Tmol/a",
    #     id="dic_weathering",
    #     species=[M.DIC],
    #     ctype="regular",
    # )

    # ConnectionProperties(
    #     source=M.Fw,
    #     sink=M.L_b,
    #     rate="24 Tmol/a",
    #     register=M,
    #     id="ta_Weathering",
    #     species=[M.TA],
    #     ctype="regular",
    # )
    ConnectionProperties(
        source=M.Fw,
        sink=M.L_b,
        rate={M.DIC: "12 Tmol/a", M.TA: "24 Tmol/a"},
        id="weathering",
        species=[M.DIC, M.TA],
        ctype="regular",
    )

    return M


# this goes at the end of the file
if __name__ == "__main__":
    from esbmtk import data_summaries
    from esbmtk import carbonate_system_2_pp

    run_time = "1000 kyr"
    time_step = "1000 yr"
    rain_ratio = 0.3
    alpha = 0.72
    alpha = 0.6

    M = initialize_esbmtk_model(rain_ratio, alpha, run_time, time_step)

    M.read_state(directory="init_data_2")
    M.run()
    # M.save_state(directory="init_data_2")

    # get CaCO3_export in mol/year
    CaCO3_export = M.CaCO3_export.to(f"{M.f_unit}").magnitude
    carbonate_system_2_pp(M.D_b, CaCO3_export, 200, 6000)

    M.save_data()

    species_names = [M.DIC, M.TA, M.CO3, M.pH, M.zcc, M.zsat, M.zsnow]
    box_names = [M.L_b, M.H_b, M.D_b]
    pl = data_summaries(M, species_names, box_names, M.L_b.DIC)
    pl += [M.CO2_At]
    M.plot(pl, fn="baseline_w_d.pdf")

    # -- test results versus pyCO2sys ---------------- #
    # update swc parameters to the concentrations at
    # the end of the model run
    M.D_b.swc.update_parameters()
    M.L_b.swc.update_parameters()
    M.H_b.swc.update_parameters()
    print()
    print(f"D_b: swc pH_total = {M.D_b.swc.pH_total:.4f}, M =  {M.D_b.pH.c[-1]:.4f}")
    print(
        f"D_B  swc co2aq = {M.D_b.swc.co2aq*1e6:.4f}, M =  {M.D_b.CO2aq.c[-1]*1.e6:.4f}"
    )
    print(f"D_b: swc co3 = {M.D_b.swc.co3*1e6:.4f}, M =  {M.D_b.CO3.c[-1]*1.e6:.4f}\n")

    print(
        f"L_B  swc co2aq = {M.L_b.swc.co2aq*1e6:.4f}, M =  {M.L_b.CO2aq.c[-1]*1.e6:.4f}"
    )
    print(f"L_b: swc co3 = {M.L_b.swc.co3*1e6:.4f}, M =  {M.L_b.CO3.c[-1]*1.e6:.4f}\n")

    print(
        f"H_B  swc co2aq = {M.H_b.swc.co2aq*1e6:.4f}, M =  {M.H_b.CO2aq.c[-1]*1.e6:.4f}"
    )
    print(f"H_b: swc co3 = {M.H_b.swc.co3*1e6:.4f}, M =  {M.H_b.CO3.c[-1]*1.e6:.4f}\n")

    # pco2_w =
    print(f"L_b: pCO2: swc pco2 = {M.L_b.swc.pCO2:.0f}, M =  {M.CO2_At.c[-1]*1.e6:.0f}")
    print(
        f"H_b: pCO2: swc pco2 = {M.H_b.swc.pCO2:.0f}, M =  {M.CO2_At.c[-1]*1.e6:.0f}\n"
    )

    print(f"zsat should = 3715, actual = {M.D_b.zsat.c[-1]:.2f}")
    print(f"zcc should = 4750, actual = {M.D_b.zcc.c[-1]:.2f}")
    print(f"zsnow should = 4750, actual = {M.D_b.zsnow.c[-1]:.2f}")

    # print(f"H_b area % = {M.H_b.area/M.hyp.area_dz(-200, -4152)")
    # M.hyp.show_data(1000,-6000)
    M.hyp.read_data()
