from __future__ import annotations
from numpy import array as npa
from numba import njit
from esbmtk import carbonate_system_1 ,gas_exchange ,carbonate_system_2

# @njit(fastmath=True)
def eqs(t, R, M, gpt, toc, area_table, area_dz_table, Csat_table) -> list:
        '''Auto generated esbmtk equations do not edit
        '''

# ---------------- write computed reservoir equations -------- #
# that do not depend on fluxes
        dCdt_M_L_b_Hplus, dCdt_M_L_b_CO2aq = carbonate_system_1(
            (R[2]),
            (R[3]),
            R[6],
            R[7],
            gpt[0],
        )

        dCdt_M_H_b_Hplus, dCdt_M_H_b_CO2aq = carbonate_system_1(
            (R[0]),
            (R[1]),
            R[8],
            R[9],
            gpt[1],
        )

        M_C_CO2_At_to_H_b_CO2_H_b__F = gas_exchange(
            (R[12]),
            (R[0]),
            (R[9]),
            gpt[3],
        )

        M_C_CO2_At_to_L_b_CO2_L_b__F = gas_exchange(
            (R[12]),
            (R[2]),
            (R[7]),
            gpt[4],
        )


# ---------------- write all flux equations ------------------- #
        M_CG_H_b_to_D_b_DIC_mix_down__F = toc[6] * R[0]
        M_CG_H_b_to_D_b_TA_mix_down__F = toc[8] * R[1]
        M_CG_D_b_to_H_b_DIC_mix_up__F = toc[10] * R[4]
        M_CG_D_b_to_H_b_TA_mix_up__F = toc[12] * R[5]
        M_CG_L_b_to_H_b_DIC_thc__F = toc[14] * R[2]
        M_CG_L_b_to_H_b_TA_thc__F = toc[16] * R[3]
        M_CG_H_b_to_D_b_DIC_thc__F = toc[18] * R[0]
        M_CG_H_b_to_D_b_TA_thc__F = toc[20] * R[1]
        M_CG_D_b_to_L_b_DIC_thc__F = toc[22] * R[4]
        M_CG_D_b_to_L_b_TA_thc__F = toc[24] * R[5]
        M_CG_L_b_to_D_b_DIC_POM__F = toc[27]
        M_CG_L_b_to_D_b_PIC_DIC__F = toc[29]
        M_CG_L_b_to_D_b_PIC_TA__F = toc[31]
        M_CG_Fw_to_L_b_DIC_weathering__F = toc[42]
        M_CG_Fw_to_L_b_TA_weathering__F = toc[44]

# ---------------- write computed reservoir equations -------- #
# that do depend on fluxes
        M_D_b_DIC_db_cs2, M_D_b_TA_db_cs2, dCdt_M_D_b_Hplus, dCdt_M_D_b_zsnow = carbonate_system_2(
            M_CG_L_b_to_D_b_PIC_DIC__F,
            (R[4]),
            (R[5]),
            (R[2]),
            R[10],
            R[11],
            gpt[2],
        )


# ---------------- write input only reservoir equations -------- #

# ---------------- write regular reservoir equations ------------ #
        dCdt_M_H_b_DIC = (
            - M_CG_H_b_to_D_b_DIC_mix_down__F
            + M_CG_D_b_to_H_b_DIC_mix_up__F
            + M_CG_L_b_to_H_b_DIC_thc__F
            - M_CG_H_b_to_D_b_DIC_thc__F
            + M_C_CO2_At_to_H_b_CO2_H_b__F
        )/toc[0]

        dCdt_M_H_b_TA = (
            - M_CG_H_b_to_D_b_TA_mix_down__F
            + M_CG_D_b_to_H_b_TA_mix_up__F
            + M_CG_L_b_to_H_b_TA_thc__F
            - M_CG_H_b_to_D_b_TA_thc__F
        )/toc[1]

        dCdt_M_L_b_DIC = (
            - M_CG_L_b_to_H_b_DIC_thc__F
            + M_CG_D_b_to_L_b_DIC_thc__F
            - M_CG_L_b_to_D_b_DIC_POM__F
            - M_CG_L_b_to_D_b_PIC_DIC__F
            + M_C_CO2_At_to_L_b_CO2_L_b__F
            + M_CG_Fw_to_L_b_DIC_weathering__F
        )/toc[2]

        dCdt_M_L_b_TA = (
            - M_CG_L_b_to_H_b_TA_thc__F
            + M_CG_D_b_to_L_b_TA_thc__F
            - M_CG_L_b_to_D_b_PIC_TA__F
            + M_CG_Fw_to_L_b_TA_weathering__F
        )/toc[3]

        dCdt_M_D_b_DIC = (
            + M_CG_H_b_to_D_b_DIC_mix_down__F
            - M_CG_D_b_to_H_b_DIC_mix_up__F
            + M_CG_H_b_to_D_b_DIC_thc__F
            - M_CG_D_b_to_L_b_DIC_thc__F
            + M_CG_L_b_to_D_b_DIC_POM__F
            + M_D_b_DIC_db_cs2
        )/toc[4]

        dCdt_M_D_b_TA = (
            + M_CG_H_b_to_D_b_TA_mix_down__F
            - M_CG_D_b_to_H_b_TA_mix_up__F
            + M_CG_H_b_to_D_b_TA_thc__F
            - M_CG_D_b_to_L_b_TA_thc__F
            + M_D_b_TA_db_cs2
        )/toc[5]

        dCdt_M_CO2_At = (
            - M_C_CO2_At_to_H_b_CO2_H_b__F
            - M_C_CO2_At_to_L_b_CO2_L_b__F
        )/toc[38]


# ---------------- write isotope reservoir equations ------------ #

# ---------------- bits and pieces --------------------------- #
        return [
            dCdt_M_H_b_DIC,  # 0
            dCdt_M_H_b_TA,  # 1
            dCdt_M_L_b_DIC,  # 2
            dCdt_M_L_b_TA,  # 3
            dCdt_M_D_b_DIC,  # 4
            dCdt_M_D_b_TA,  # 5
            dCdt_M_L_b_Hplus,  # 6
            dCdt_M_L_b_CO2aq,  # 7
            dCdt_M_H_b_Hplus,  # 8
            dCdt_M_H_b_CO2aq,  # 9
            dCdt_M_D_b_Hplus,  # 10
            dCdt_M_D_b_zsnow,  # 11
            dCdt_M_CO2_At,  # 12
        ]
