import numpy as np
from defs import *
from values import *

# compute electro-convective-diffusive fluxes

def compute_ecd_fluxes (cell,jvol):

    # function return values
    jsol = np.zeros([NS,NC,NC])
    delmu = np.zeros([NS,NC,NC])

    
    # local variables
    dmu = np.zeros([NS,NC])
    if cell.segment == 'PT' or cell.segment == 'S3' or cell.segment =='SDL' or cell.segment =='LDL' or cell.segment =='LAL' or cell.segment == 'mTAL' or cell.segment == 'cTAL' or cell.segment == 'MD' or cell.segment == 'DCT' or cell.segment == 'IMCD':
        crange=[0,1,4,5]
    elif cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
        crange=[0,1,2,3,4,5]
    for i in range(NS):
        for k in crange:
            dmu[i,k] = RT*np.log(abs(cell.conc[i,k]))+zval[i]*F*EPref*cell.ep[k]
            
    
    for i in range(NS):
        for k in crange:
            for l in crange[crange.index(k):]:

                XI = zval[i]*F*EPref/RT*(cell.ep[k]-cell.ep[l])                
                dint = np.exp(-XI)
               
                if (abs(dint-1)<1e-6):
                    jsol[i][k][l] = cell.area[k][l]*cell.h[i][k][l]*(cell.conc[i][k]-cell.conc[i][l])

                else:
                    jsol[i][k][l] = cell.area[k][l]*cell.h[i][k][l]*XI*(cell.conc[i][k]-cell.conc[i][l]*dint)/(1-dint)             

                # convective component
                concdiff = cell.conc[i][k]-cell.conc[i][l]
                if (abs(concdiff)>1e-6):
                    concmean=concdiff/np.log(abs(cell.conc[i][k]/cell.conc[i][l]))
                    dimless=(Pfref*Vwbar*Cref)/href
                    convect=(1.0e0-cell.sig[i][k][l])*concmean*jvol[k][l]*dimless
                    jsol[i][k][l]=jsol[i][k][l]+convect

                # define driving force for coupled fluxes
                delmu[i][k][l] = dmu[i][k]-dmu[i][l]
                
    return jsol,delmu

        
