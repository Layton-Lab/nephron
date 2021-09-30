#---------------------------------------------------------------------72
#---------------------------------------------------------------------72
#    COMPUTE WATER FLUXES
#---------------------------------------------------------------------72
#---------------------------------------------------------------------72
#     The hydraulic and oncotic pressures are made non-dimensional
#     by dividing by RT*Cref

from values import *
from defs import *
import numpy as np
from set_params import set_torq_params

def compute_water_fluxes (cell):

    # need to take care of these constants
    if cell.segment == 'CNT' or cell.segment == 'CCD':
        compl = 0.3
    else:
        compl = 0.1
    # also, should define range of relevant compartments for each cell

    # local variables
    PRES = np.zeros(NC)
    ONC  = np.zeros(NC)
    jvol = np.zeros(NC*NC).reshape((NC,NC))

    PM=cell.pres[0]
  
    
    if cell.segment == 'PT' or cell.segment == 'S3':
        Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(cell.species,cell.sex,cell.preg)
        PB = PbloodPT        
    else:
        PB=0

    PRES[0] = PM/(RTosm*Cref)
    PRES[1] = PRES[0]
    PRES[2] = PRES[0]
    PRES[3] = PRES[0]
    PRES[5] = PB/(RTosm*Cref)
    PRES[4] = PRES[5]+(cell.vol[4]/cell.volref[4]-1)/compl/(RTosm*Cref) 
    if cell.species == 'rat' or cell.species == 'mou':
        if cell.segment == 'SDL'or cell.segment == 'LDL' or cell.segment == 'LAL':
            PRES[1] = PRES[5]

    ONC[0]=cell.cimpref[0]
    ONC[1]=cell.cimpref[1]*cell.volref[1]/cell.vol[1]
    ONC[4]=cell.cimpref[4]
    ONC[5]=cell.cimpref[5]
    if cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
        ONC[2]=cell.cimpref[2]*cell.volref[1]/cell.vol[2]
        ONC[3]=cell.cimpref[3]*cell.volref[1]/cell.vol[3]    

    for k in range(NC-1):
        for l in range(k+1,NC):
            osm = sum(cell.sig[:,k,l]*(cell.conc.T[k,:]-cell.conc.T[l,:]))
                        
            jvol[k][l] = cell.area[k][l] * cell.dLPV[k][l] * (PRES[k]-PRES[l]-ONC[k]+ONC[l]-osm)

            jvol[l][k] = -jvol[k][l]

            # if (cell.flag==2) and (i==0) and (j==1):
            #     print(cell.area[i][j] * cell.dLPV[i][j],'\t',PRES[i]-PRES[j],'\t',-ONC[i]+ONC[j],'\t',-osm)
            #     print('\n')
                
            # if (cell.flag==2) and (i==0) and (j==4):
            #     print(cell.area[i][j] * cell.dLPV[i][j],'\t',PRES[i]-PRES[j],'\t',-ONC[i]+ONC[j],'\t',-osm)
            #     print('\n')

    return jvol


    
