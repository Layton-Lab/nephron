from defs import *
from values import *
import numpy as np
import math

def calculate_conduct(cell):
    gf_star = 59 #ps
    gs_star = 29 #ps

    ph = np.zeros(NC)
    if cell.segment == 'DCT':
        for j in [0,1,4,5]:
            ph[j] = -np.log(abs(cell.conc[11][j])/1.0e3)/np.log(10.0)
    else:
        for j in range(NC):
            ph[j] = -np.log(abs(cell.conc[11][j])/1.0e3)/np.log(10.0)

    pHi = ph[1]
    pHe = ph[0]

    gf = gf_star + \
        (gf_star/(gf_star+gs_star))*((91-58)/(7.4-5.4)) *\
        (pHe - 7.4)

    gs = gs_star + \
        (gs_star/(gf_star+gs_star))*((91-58)/(7.4-5.4)) *\
        (pHe - 7.4)

    k1_t = 42.7 

    if pHi < 7.4:
        k2_t = 55.9 + ((173.3-55.9)/(7.0-7.4))*(pHi-7.4)
    else:
        k2_t = 55.9 + ((30.4-55.9)/(8.4-7.4))*(pHi-7.4)

    k3_t = 0.1684*np.exp(0.6035*pHi)
    k4_t = 58.7


    Ps = 1/(1+ k2_t/k1_t + k3_t/k4_t)
    Pf = (k3_t/k4_t)*Ps

    g_trpv5 = Pf*gf + Ps*gs
    

    return g_trpv5

def trpv5(cell,i,ep,memb_id,act,area):
    


    N_trpv5 = act # ** value based on experiment 


    #g_trpv5 = 59 # pS
    F_si = 9.6485e4

    g_trpv5 = calculate_conduct(cell)

    delta_vol = (cell.ep[0] - cell.ep[1])*EPref
    E_nernst = (RT*np.log(abs(cell.conc[15,1]/cell.conc[15,0])))/(2*F)
    C_inhib = 74e-6

    f_trpv5 = 1/(1+ (cell.conc[15,1]/C_inhib)) 

    flux_trpv5 = f_trpv5 * ((N_trpv5 * g_trpv5)*1e-9)/(Cref*href) * ( delta_vol - E_nernst)/(2*F_si)

    # Luminal [Ca2+] Effect:
    f_con = 1 + (cell.conc[15, 0] / 1.98) ** 2.37

    fluxCa = area[memb_id[0]][memb_id[1]] * flux_trpv5 * f_con

    # at DCT2 and CNT only
    if (cell.segment == 'DCT' and i < 2.0/3.0*cell.total):
        fluxCa = 0



    return [15], [fluxCa]

