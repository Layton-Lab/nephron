#
# Glucose transporters, include SGLT1/SGLT2, GLUT1/GLUT2

from values import *
import math
import numpy as np

#
# SGLT cotransporters modules
#---------------------------------------------------------------------
#This subroutine to computes the fluxes across SGLT (Na-glucose)
#cotransporters, based on the model of Parent et al. (J Memb Biol 1992)
#with the modifications described by Eskandari et al. (J Memb Biol 2005)
#---------------------------------------------------------------------
#   conc[NS][NC]: all concentrations in all compartments
#   ep[NC]:       all potential in all compartments
#   memb_id:      specifies interface, (lumen,cell), (cell,bath), or (cell,LIS)
#   CT:           activity level
#   area:         surface area of all membranes
#
#   returns:      solute IDs and corresponding fluxes
#

def sglt1(cell,ep,memb_id,CT,area):

    if cell.diabete != 'Non':
        CT = CT*0.67

    # if cell.segment == 'PT':
    #     sglt1 = 0.0
    # elif cell.segment == 'S3':
    #     sglt1 = 1.0
    #assign concentrations (in M, not mM)
    nao = cell.conc[0][memb_id[0]]*1.0e-3 
    nai = cell.conc[0][memb_id[1]]*1.0e-3
    gluo = cell.conc[14][memb_id[0]]*1.0e-3
    glui = cell.conc[14][memb_id[1]]*1.0e-3
    pot = F*(ep[memb_id[1]]-ep[memb_id[0]])*EPref/RT

    #The reaction rates for SGLT1 come from Wright et al. Physiol Rev 2011

    delta = 0.70
    alphap = 0.30
    alphapp = 0.00
    dk12 = 140000.0 * math.exp(-pot*alphap)*(nao**2)
    dk21 = 300.0 * math.exp(+pot*alphap)
    dk23 = 45000.0*gluo
    dk32 = 20.0
    dk34 = 50.0
    dk43 = 50.0
    dk45 = 800.0
    dk54 = 190000.0*glui 
    dk56 = 5.0 * math.exp(-pot*alphapp)
    dk65 = 2250.0 * math.exp(+pot*alphapp)*(nai**2) 
    dk61 = 25.0 * math.exp(-pot*delta) 
    dk16 = 600.0 * math.exp(+pot*delta) 
    dk25 = 0.01
    dk52 = 0.0005
    nstoich = 2

    #Compute the KAT terms and the concentrations of each enzyme state

    A = np.mat([[dk12+dk16+dk61,dk61-dk21,dk61,dk61,dk61],
                [-dk12,dk21+dk23+dk25,-dk32,0,-dk52],
                [0,-dk23,dk32+dk34,-dk43,0],
                [0,0,-dk34,dk43+dk45,-dk54],
                [dk65,dk65-dk25,dk65,dk65-dk45,dk52+dk54+dk56+dk65]])
    b = np.array([dk61,0,0,0,dk65])
    x = np.linalg.solve(A,b)

    C1=x[0]
    C2=x[1]
    C3=x[2]
    C4=x[3]
    C5=x[4]
    C6=1.0-x[0]-x[1]-x[2]-x[3]-x[4]

    #Compute the sodium and glucose fluxes

    fluxsglt=[0,0]
    fluxsglt[0] = area[memb_id[0]][memb_id[1]]*CT*nstoich*(dk34*C3-dk43*C4+dk25*C2-dk52*C5)#*sglt1
    fluxsglt[1] = area[memb_id[0]][memb_id[1]]*CT*1.0*(dk34*C3-dk43*C4)#*sglt1
    #print(fluxsglt[0],fluxsglt[1])
    return [0,14],fluxsglt

#---------------------------------------------------------------------72
#   Compute SGLT2 fluxes
#---------------------------------------------------------------------72
#   NET formulation from Alan's 2007 PT model
#---------------------------------------------------------------------72
#   conc[NS][NC]: all concentrations in all compartments
#   ep[NC]:       all potential in all compartments
#   memb_id:      specifies interface, (lumen,cell), (cell,bath), or (cell,LIS)
#   CT:           activity level
#   area:         surface area of all membranes
#
#   returns:      solute IDs and corresponding fluxes
#
def sglt2(cell,ep,memb_id,CT,area):

    if cell.diabete != 'Non':
        CT = CT*1.38

    zeta = F*EPref*(ep[memb_id[0]]-ep[memb_id[1]])/RT
    affglu=4.90e0
    affna=25.0e0
    nal = cell.conc[0][memb_id[0]]/affna
    nac = cell.conc[0][memb_id[1]]/affna
    glul = cell.conc[14][memb_id[0]]/affglu
    gluc = cell.conc[14][memb_id[1]]/affglu
    denom = (1.0 + nal + glul + nal*glul)*(1.0 + nac*gluc) + (1.0 + nac + gluc + nac*gluc)*(1.0 + nal*glul*math.exp(zeta))
    

    fluxsglt=[0,0]
    fluxsglt[0] = (area[memb_id[0]][memb_id[1]]*CT*(cell.conc[0][0]*cell.conc[14][0]*math.exp(zeta) - cell.conc[0][1]*cell.conc[14][1])/denom)
    fluxsglt[1] = fluxsglt[0]

    return [0,14],fluxsglt


#---------------------------------------------------------------------72
#   Compute GLUT1 fluxes
#---------------------------------------------------------------------72

def glut1(cell,ep,memb_id,CT,area):

    affglut1 = 2.0
    Ro5 = affglut1*(cell.conc[14][memb_id[0]]-cell.conc[14][memb_id[1]])/(affglut1+cell.conc[14][memb_id[1]])/(affglut1+cell.conc[14][memb_id[0]])
    fluxglut1 = area[memb_id[0]][memb_id[1]]*CT*Ro5
    return 14,fluxglut1

#---------------------------------------------------------------------72
#   Compute GLUT2 fluxes
#---------------------------------------------------------------------72

def glut2(cell,ep,memb_id,CT,area):

    affglut2 = 17.0
    Ro5 = affglut2*(cell.conc[14][memb_id[0]]-cell.conc[14][memb_id[1]])/(affglut2+cell.conc[14][memb_id[1]])/(affglut2+cell.conc[14][memb_id[0]])
    fluxglut2 = area[memb_id[0]][memb_id[1]]*CT*Ro5

    return 14,fluxglut2
