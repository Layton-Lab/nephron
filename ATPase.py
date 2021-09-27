#---------------------------------------------------------------------72
#---------------------------------------------------------------------72
#    ACTIVE TRANSPORT VIA ATPases
#---------------------------------------------------------------------72
#---------------------------------------------------------------------72

from defs import *
from values import *
import numpy as np
import math

#---------------------------------------------------------------------72
#    Na-K-ATPase
#---------------------------------------------------------------------
#   conc[NS][NC]: all concentrations in all compartments
#   ep[NC]:       all potential in all compartments
#   memb_id:      specifies interface, (lumen,cell), (cell,bath), or (cell,LIS)
#   CT:           activity level
#   area:         surface area of all membranes
#
#   returns:      solute IDs and corresponding fluxes
#

def nakatpase(cell,ep,memb_id,act,area):
    if cell.diabete != 'Non':
        if cell.segment == 'PT' or cell.segment == 'S3' or cell.segment == 'cTAL' or cell.segment == 'DCT' or cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
            act = act*1.1
        elif cell.segment == 'mTAL':
            act = act*1.5
        elif cell.segment == 'IMCD':
            act = act*2.5
    
    if cell.segment == 'MD':
        AffNa = 1+cell.conc[1][memb_id[0]]/10
    else:
        AffNa = 0.2*(1.0e0+cell.conc[1][memb_id[0]]/8.33e0)
    actNa = cell.conc[0][memb_id[0]]/(cell.conc[0][memb_id[0]]+AffNa)

    if cell.segment == 'PT' or cell.segment == 'S3':
        AffK = 0.1*(1.0+cell.conc[0][memb_id[1]]/18.5e0)
    else:
        AffK = 0.1*(1.0+cell.conc[0][5]/18.5e0)
    
    if cell.segment == 'CNT' or cell.segment == 'CCD' or cell.segment == 'OMCD':
        AffNH4 = AffK/0.2
    elif cell.segment == 'IMCD':
        AffNH4 = 5*AffK
    else:
        AffNH4 =  AffK
    
    if cell.segment == 'PT' or cell.segment == 'S3':
        actK5 = (cell.conc[1][memb_id[1]]+cell.conc[10][memb_id[1]])/(cell.conc[1][memb_id[1]]+cell.conc[10][memb_id[1]]+AffK)
    else:
        actK5 = (cell.conc[1][memb_id[1]])/(cell.conc[1][memb_id[1]]+AffK)

    ro5 = (cell.conc[10][memb_id[1]]/AffNH4)/(cell.conc[1][memb_id[1]]/AffK)
    
       
    #print('activity of NaKATPase = %f' %act)
    dJactNa5 = area[memb_id[0]][memb_id[1]]*act*(actNa**3.0e0)*(actK5**2.0e0)
    dJactK5 = -2.0e0/3.0e0*dJactNa5/(1.e0+ro5)

    fluxNaKsod = dJactNa5
    fluxNaKpot = dJactK5
    fluxNaKamm = dJactK5*ro5
    
    return [0,1,10],[fluxNaKsod,fluxNaKpot,fluxNaKamm]

#---------------------------------------------------------------------72
#     H-ATPase
#     See Strieter & Weinstein paper for signs (AJP 263, 1992)
#---------------------------------------------------------------------72

def hatpase(cell,ep,memb_id,act,area):

    # parameters
    if cell.segment == 'PT' or cell.segment == 'S3':
        dmuATPH = 1.45   # PT value, maybe not for other segments
        steepATPH = 0.40   # PT value, maybe not for other segments
    elif cell.segment == 'CNT' or cell.segment == 'CCD':
        dmuATPH = 2.1
        if memb_id[0]==0:
            steepATPH=0.4
        elif memb_id[0]==3:
            steepATPH=-0.4
    elif cell.segment == 'OMCD':
        dmuATPH = 2.1
        steepATPH = 0.4

    dmu0 = RT*np.log(abs(cell.conc[11,memb_id[0]]))+zval[11]*F*EPref*ep[memb_id[0]]
    dmu1 = RT*np.log(abs(cell.conc[11,memb_id[1]]))+zval[11]*F*EPref*ep[memb_id[1]]
    delmu = dmu0-dmu1
    DactH=1.0+math.exp(steepATPH*(delmu-dmuATPH))
    
    if memb_id[0]==0:
        fluxHATPase=-area[memb_id[0]][memb_id[1]]*act/DactH
    elif memb_id[0]==3:
        fluxHATPase=area[memb_id[0]][memb_id[1]]*act/DactH
    else:
        fluxHATPase=-area[memb_id[0]][memb_id[1]]*act/DactH
    return [11],[fluxHATPase]
#------------------------------------------------------------------------
#    H-K-ATPase
#------------------------------------------------------------------------

def hkatpase(cell,memb_id,act,area):

    hkconc=np.zeros(4)
    hkconc[0]=cell.conc[1,memb_id[1]]
    hkconc[1]=cell.conc[1,memb_id[0]]
    hkconc[2]=cell.conc[11,memb_id[1]]
    hkconc[3]=cell.conc[11,memb_id[0]]

    Amat = np.matrix(fatpase(Natp,hkconc))

    IAmat = Amat.I

    c7=IAmat[6,0]
    c8=IAmat[7,0]

    dkf5=40.0
    dkb5=200.0

    hefflux=cell.area[memb_id[0],memb_id[1]]*act*(dkf5*c7-dkb5*c8)
    if cell.segment=='IMCD':
        coeff=1.0
    else:
        coeff=2.0
    return [1,11],[coeff*hefflux,-coeff*hefflux]

def fatpase(n,hkconc):
    Amat=np.zeros([n,n])
    Amat[0,]=1

    dkf1 = 1.3e7
    dkb1 = 6.5
    dkf2a = 8.9e3
    dkb2a = 7.3e4
    dkf2b = 8.9e3
    dkb2b = 7.3e4
    dkf3a = 5.3e9
    dkb3a = 6.6e2
    dkf3b = 5.3e9
    dkb3b = 6.6e2
    dkf4 = 5.0e1
    dkb4 = 2.5e6
    dkf5 = 4.0e1
    dkb5 = 2.0e2
    dkf6a = 5.0e7
    dkb6a = 8.0e12
    dkf6b = 5.0e7
    dkb6b = 8.0e12
    dkf7a = 2.6e10
    dkb7a = 1.5e8
    dkf7b = 2.6e10
    dkb7b = 1.5e8
    dkf8 = 5.4e1
    dkb8 = 3.2e1
    dkf9 = 1.75
    dkb9 = 3.5e1
    dkf10 = 5.0e4
    dkb10 = 5.0e1
    dkf11 = 5.0e2
    dkb11 = 5.0

    catp = 2.0e-3
    cadp = 0.04e-3
    cpi = 5.0e-3

    kin = hkconc[0]*1e-3
    kout = hkconc[1]*1e-3
    hin = hkconc[2]*1e-3
    hout = hkconc[3]*1e-3

    Amat[1,1] = dkf2a
    Amat[1,2] = -(dkb2a*kin+dkf2b)
    Amat[1,3] = dkb2b*kin
    Amat[2,2] = dkf2b
    Amat[2,3] = -(dkb2b*kin+dkf3a*hin)
    Amat[2,4] = +dkb3a
    Amat[3,3] = +dkf3a*hin
    Amat[3,4] = -(dkb3a+dkf3b*hin)
    Amat[3,5] = +dkb3b
    Amat[4,4] = +dkf3b*hin
    Amat[4,5] = -(dkb3b+dkf4)
    Amat[4,6] = +dkb4*cadp
    Amat[5,5] = +dkf4
    Amat[5,6] = -(dkb4*cadp+dkf5)
    Amat[5,7] = +dkb5
    Amat[6,6] = +dkf5
    Amat[6,7] = -(dkb5+dkf6a)
    Amat[6,8] = dkb6a*hout
    Amat[7,7] = +dkf6a 
    Amat[7,8] = -(dkb6a*hout+dkf6b) 
    Amat[7,9] = +dkb6b*hout
    Amat[8,8] = +dkf6b 
    Amat[8,9] = -(dkb6b*hout+dkf7a*kout) 
    Amat[8,10] = +dkb7a
    Amat[9,9] = +dkf7a*kout 
    Amat[9,10] = -(dkb7a+dkf7b*kout) 
    Amat[9,11] = +dkb7b
    Amat[10,10] = +dkf7b*kout 
    Amat[10,11] = -(dkb7b+dkf8) 
    Amat[10,12] = +dkb8*cpi
    Amat[11,0] =  +dkb9 
    Amat[11,11] = +dkf8 
    Amat[11,12] = -(dkb8*cpi+dkf9+dkf10*catp)
    Amat[11,13] = +dkb10
    Amat[12,1] = +dkb11 
    Amat[12,12] = +dkf10*catp 
    Amat[12,13] = -(dkb10+dkf11) 
    Amat[13,0] = +dkf1*catp
    Amat[13,1] = -(dkb1+dkf2a+dkb11)
    Amat[13,2] = +dkb2a*kin
    Amat[13,13] = dkf11
    
    return Amat