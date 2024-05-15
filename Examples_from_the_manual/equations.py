from __future__ import annotations
from numpy import array as npa
from numba import njit

# @njit(fastmath=True)
def eqs(t, R, M, gpt, toc, area_table, area_dz_table, Csat_table) -> list:
        '''Auto generated esbmtk equations do not edit
        '''

# ---------------- write computed reservoir equations -------- #
# that do not depend on fluxes

# ---------------- write all flux equations ------------------- #
        M_CG_S_b_to_D_b_thc_up_DIC_thc_up__F = toc[4] * R[0]
        M_CG_S_b_to_D_b_thc_up_DIC_thc_up__F_l =  toc[4] * R[0] * R[1] / R[0]
        M_CG_S_b_to_D_b_thc_up_PO4_thc_up__F = toc[5] * R[2]
        M_CG_D_b_to_S_b_thc_down_DIC_thc_down__F = toc[6] * R[3]
        M_CG_D_b_to_S_b_thc_down_DIC_thc_down__F_l =  toc[6] * R[3] * R[4] / R[3]
        M_CG_D_b_to_S_b_thc_down_PO4_thc_down__F = toc[7] * R[5]
        M_CG_weathering_to_S_b_weathering_PO4_weathering__F = toc[9]
        M_CG_weathering_to_S_b_weathering_DIC_weathering__F = toc[11]
        M_CG_weathering_to_S_b_weathering_DIC_weathering__F_l =  toc[11] * 1000 / (0.0112372 * (0 + 1000) + 1000)
        M_CG_S_b_to_D_b_primary_production_PO4_primary_production__F = toc[13] * R[2]
        M_CG_S_b_to_D_b_OM_production_DIC_OM_production__F = toc[14] * R[2]
        M_CG_S_b_to_D_b_OM_production_DIC_OM_production__F_l =  R[1] * toc[14] * R[2] / (0.972 * R[0] + R[1] - 0.972 * R[1])
        M_CG_D_b_to_burial_burial_DIC_burial__F = toc[16] * M_CG_S_b_to_D_b_primary_production_PO4_primary_production__F
        M_CG_D_b_to_burial_burial_DIC_burial__F_l =  R[4] * toc[16] * M_CG_S_b_to_D_b_primary_production_PO4_primary_production__F / (1.0 * R[3] + R[4] - 1.0 * R[4])
        M_CG_D_b_to_burial_burial_PO4_burial__F = toc[18] * M_CG_S_b_to_D_b_primary_production_PO4_primary_production__F

# ---------------- write computed reservoir equations -------- #
# that do depend on fluxes

# ---------------- write input only reservoir equations -------- #

# ---------------- write regular reservoir equations ------------ #
        dCdt_M_S_b_DIC = (
            - M_CG_S_b_to_D_b_thc_up_DIC_thc_up__F
            + M_CG_D_b_to_S_b_thc_down_DIC_thc_down__F
            + M_CG_weathering_to_S_b_weathering_DIC_weathering__F
            - M_CG_S_b_to_D_b_OM_production_DIC_OM_production__F
        )/toc[0]

        dCdt_M_S_b_PO4 = (
            - M_CG_S_b_to_D_b_thc_up_PO4_thc_up__F
            + M_CG_D_b_to_S_b_thc_down_PO4_thc_down__F
            + M_CG_weathering_to_S_b_weathering_PO4_weathering__F
            - M_CG_S_b_to_D_b_primary_production_PO4_primary_production__F
        )/toc[1]

        dCdt_M_D_b_DIC = (
            + M_CG_S_b_to_D_b_thc_up_DIC_thc_up__F
            - M_CG_D_b_to_S_b_thc_down_DIC_thc_down__F
            + M_CG_S_b_to_D_b_OM_production_DIC_OM_production__F
            - M_CG_D_b_to_burial_burial_DIC_burial__F
        )/toc[2]

        dCdt_M_D_b_PO4 = (
            + M_CG_S_b_to_D_b_thc_up_PO4_thc_up__F
            - M_CG_D_b_to_S_b_thc_down_PO4_thc_down__F
            + M_CG_S_b_to_D_b_primary_production_PO4_primary_production__F
            - M_CG_D_b_to_burial_burial_PO4_burial__F
        )/toc[3]


# ---------------- write isotope reservoir equations ------------ #
        dCdt_M_S_b_DIC_l = (
            - M_CG_S_b_to_D_b_thc_up_DIC_thc_up__F_l
            + M_CG_D_b_to_S_b_thc_down_DIC_thc_down__F_l
            + M_CG_weathering_to_S_b_weathering_DIC_weathering__F_l
            - M_CG_S_b_to_D_b_OM_production_DIC_OM_production__F_l
        )/toc[0]

        dCdt_M_D_b_DIC_l = (
            + M_CG_S_b_to_D_b_thc_up_DIC_thc_up__F_l
            - M_CG_D_b_to_S_b_thc_down_DIC_thc_down__F_l
            + M_CG_S_b_to_D_b_OM_production_DIC_OM_production__F_l
            - M_CG_D_b_to_burial_burial_DIC_burial__F_l
        )/toc[2]


# ---------------- bits and pieces --------------------------- #
        return [
            dCdt_M_S_b_DIC,  # 0
            dCdt_M_S_b_DIC_l,  # 1
            dCdt_M_S_b_PO4,  # 2
            dCdt_M_D_b_DIC,  # 3
            dCdt_M_D_b_DIC_l,  # 4
            dCdt_M_D_b_PO4,  # 5
        ]
