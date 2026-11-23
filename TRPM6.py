from defs import *
from values import *
import numpy as np
import math


def trpm6(cell, ep, memb_id, act, area):
    F_si = 9.6485e4

    delta_vol = (cell.ep[0] - cell.ep[1]) * EPref
    E_nernst = (RT * np.log(abs(cell.conc[16, 1] / cell.conc[16, 0]))) / (2 * F)
    Mg_inhib = 0.51

    # Intracellular Mg2+ has an inhibitory effect on TRPM6 activity
    trpm6_Mgi = 1 / (1 + (cell.conc[16,1] / Mg_inhib) ** 2)

    # Lowering extracellular pH increases TRPM6 inward currents
    # single channel conductance of TRPM6 at pH7.4 = 83.6 pS
    trpm6_pH = 56.6 * (2 - 1 / (1 + ((7.4 - cell.pH[0]) / (7.4 - 5.5))))

    flux_trpm6 = trpm6_Mgi * trpm6_pH  * (delta_vol - E_nernst) / (2 * F_si)

    # Concentration-dependence
    f_con = 1 + (cell.conc[16, 0]) / (0.8 + cell.conc[16, 0])

    fluxMg = area[memb_id[0]][memb_id[1]] * act * flux_trpm6 * f_con

    return [16], [fluxMg]