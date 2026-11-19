from defs import *
from values import *
import numpy as np
import math

#%%
def ncx1(cell,i, ep,memb_id,act,area):
    
    Km_NCX1 = 0.125e-3 

    Km_Nao = 87.5 ##
    Km_Cao = 1.3 ##

    Km_Nai = 12.29 ## 
    Km_Cai = 3.59e-3 ##

    k_sat = 0.27
    gamma = 0.35

    Phi_R = np.exp((gamma-1)*cell.ep[1]*EPref*F/RT)
    Phi_F = np.exp(gamma*cell.ep[1]*EPref*F/RT)

    m = memb_id[1]

    G = (cell.conc[0,m]**3)*(cell.conc[15,1]) + \
        (cell.conc[0,1]**3)*(cell.conc[15,m]) + \
        (Km_Nao**3)*cell.conc[15,1] + \
        (Km_Cao)*(cell.conc[0,1]**3) + \
        (Km_Nai**3)*(cell.conc[15,m])*(1+(cell.conc[15,1]/Km_Cai)) + \
        (Km_Cai)*(cell.conc[0,m]**3)*(1+(cell.conc[0,1]**3/Km_Nai**3))

    cell_term = (cell.conc[15,1]/Km_NCX1)

    fluxNCX1_Ca = (cell_term**2/(1+cell_term**2)) * \
                  (Phi_R*(cell.conc[0,m]**3)*cell.conc[15,1] - \
                   Phi_F*(cell.conc[0,1]**3)*cell.conc[15,m])/ \
                   (G*(1+k_sat*Phi_R))

    JNCX1_Ca = area[memb_id[0]][memb_id[1]] * act * fluxNCX1_Ca

    if (cell.segment == 'DCT') and (i < 2.0/3.0*cell.total):
        JNCX1_Ca = 0
        
    JNCX1_Na = -3*JNCX1_Ca


    return [0,15], [JNCX1_Na, JNCX1_Ca]

#%%