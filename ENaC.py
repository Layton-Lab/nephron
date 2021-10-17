import numpy as np
from values import *
from defs import *

def ENaC(cell,i,memb_id,hNaMP,area,jvol):
# hNaMP is the activity of ENaC from .dat file
    
    facNaMP=(30.0/(30.0+cell.conc[0,0]))*(50.0/(50.0+cell.conc[0,1]))
    
    
    if cell.segment =='DCT':
        # Partial overlap of NCC and ENaC in DCT
        LzD=i-1
        x1=LzD+1
        xn=cell.total
        xdct2=2.0/3.0

        if x1/xn<xdct2:
            ENaCexp=0.0
        else:
            ENaCexp=(x1/xn-xdct2)/(1-xdct2)            
        #Dependence of Lumen-Cell permeability to Na (ENac) on Na concentrations and on flow
        facphMP=1.0
               
        hENaC=ENaCexp*hNaMP*facNaMP*facphMP

    elif cell.segment == 'CNT':
        #Dependence of Lumen-Cell permeability to Na (ENac) on Na concentrations and on flow
        facCaMP=1.0
        facphMP=1.0
        if cell.species == 'rat':
            flow_ref = 2.0e-6

            NaMPq0=cell.vol_init[0]-(flow_ref)/60/Vref
        elif cell.species == 'mou':
            NaMPq0=cell.vol_init[0]-(1.6e-6)/60/Vref 
        elif cell.species == 'hum':
            NaMPq0=cell.vol_init[0]-(2.0e-6)/60/Vref
        facFvMP=max(0.01,1+3*((cell.vol[0]/NaMPq0)-1))
        
        # # tracking
        # fname = 'tracking_all.txt'
        # f0 = open(fname, 'a')
        # f0.write('cell.vol_init[0]: ' + str(cell.vol_init[0]) + ', flow_ref/60/Vref: ' + str((flow_ref/60/Vref)) + ', cell_vol[0]: '+ str(cell.vol[0]) + ', NaMPq0: '+str(NaMPq0)+ ', cell_vol[0]/NaMPq0: '+str(cell.vol[0]/NaMPq0) + '\n')
        # f0.close()


        # fname1 = 'tracking_NaMPq0.txt'
        # f1 = open(fname1, 'a')
        # f1.write(str(NaMPq0)+'\n')
        # f1.close()

        # fname2 = 'tracking_facFvMP.txt'
        # f2 = open(fname2, 'a')
        # f2.write(str(facFvMP) + '\n')
        # f2.close()

        hENaC=hNaMP*facNaMP*facCaMP*facphMP*facFvMP

    elif cell.segment == 'CCD':
        #Dependence of Lumen-Cell permeability to Na (ENac) on Na concentrations and on flow
        facphMP=1.0
        if cell.species == 'rat':
            flow_ref = 0.1e-6

            NaMPq0=cell.vol_init[0]-(flow_ref)/60/Vref
        elif cell.species == 'mou':
            NaMPq0=cell.vol_init[0]-(0.08e-6)/60/Vref
        elif cell.species == 'hum':
            NaMPq0=cell.vol_init[0]-(0.1e-6)/60/Vref
        facFvMP=max(0.01,1+3*((cell.vol[0]/NaMPq0)-1))
                
        hENaC=hNaMP*facNaMP*facphMP*facFvMP

    elif cell.segment == 'OMCD':
        #Dependence of Lumen-Cell permeability to Na (ENac) on Na concentrations and on flow
        facphMP=1.0
        
        hENaC=hNaMP*facNaMP*facphMP

    XI=zval[0]*F*EPref/RT*(cell.ep[0]-cell.ep[1])
    dint=np.exp(-XI)
    if (abs(1-dint)<1e-6):
        dJENaC=cell.area[0,1]*hENaC*(cell.conc[0,0]-cell.conc[0,1])
    else:
        dJENaC=cell.area[0,1]*hENaC*XI*(cell.conc[0,0]-cell.conc[0,1]*dint)/(1-dint)

    if cell.segment=='CNT' or cell.segment == 'CCD':
        concdiff=cell.conc[0,0]-cell.conc[0,1]
        if abs(concdiff)>1e-6:
            concmean=(cell.conc[0,0]-cell.conc[0,1])/np.log(abs(cell.conc[0,0]/cell.conc[0,1]))
            dimless=(Pfref*Vwbar*Cref)/href
            convect=(1-cell.sig[0,0,1])*concmean*jvol[0,1]*dimless
            dJENaC=dJENaC+convect

    return [0],[dJENaC]
