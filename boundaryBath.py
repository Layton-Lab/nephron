from values import *
import numpy as np
from defs import *

def female_conc(cell, i):
    # Concentrations of K, Cl are lower in female rat.
    cell.conc[1,5]=cell.conc[1,5]-1
    cell.conc[2,5]=cell.conc[2,5]-1
    if cell.segment == 'PT':
        if cell.preg == 'mid':
            # conc of Na, K, Cl, HCO3 different in MP rat plasma
            cell.conc[0,0] = cell.conc[0,0]*0.94
            cell.conc[1,0] = cell.conc[1,0]*1.0 #no change in K in MP
            cell.conc[2,0] = cell.conc[2,0]*0.94
            cell.conc[3,0] = cell.conc[3,0]*0.94

            # also change in starting for bath
            cell.conc[0,5] = cell.conc[0,0]
            cell.conc[1,5] = cell.conc[1,0] 
            cell.conc[2,5] = cell.conc[2,0]
            cell.conc[3,5] = cell.conc[3,0]
        elif cell.preg == 'late':
            # conc of Na, K, Cl, HCO3 different in LP rat plasma
            cell.conc[0,0] = cell.conc[0,0]*0.95
            cell.conc[1,0] = cell.conc[1,0]*1.15 # upregulated K in LP
            cell.conc[2,0] = cell.conc[2,0]*0.95
            cell.conc[3,0] = cell.conc[3,0]*0.95

            # also change in starting for bath
            cell.conc[0,5] = cell.conc[0,0]
            cell.conc[1,5] = cell.conc[1,0]
            cell.conc[2,5] = cell.conc[2,0]
            cell.conc[3,5] = cell.conc[3,0]
    

def boundaryBath(cell,i):
    if cell.species == 'rat':
        if cell.inhib == 'NHE3-50':
            cell.oi[0] = TotSodOI_noinhib
            cell.oi[1] = TotPotOI_noinhib
            cell.oi[2] = TotCloOI_noinhib
            cell.oi[8] = TotureaOI_50inhib
            cell.pap[0] = TotSodPap_50inhib
            cell.pap[1] = TotPotPap_noinhib
            cell.pap[2] = TotCloPap_50inhib
            cell.pap[8] = TotureaPap_50inhib
        elif cell.inhib == 'NHE3-80':
            cell.oi[0] = TotSodOI_noinhib
            cell.oi[1] = TotPotOI_noinhib
            cell.oi[2] = TotCloOI_noinhib		
            cell.oi[8] = TotureaOI_80inhib
            cell.pap[0] = TotSodPap_80inhib
            cell.pap[1] = TotPotPap_80inhib
            cell.pap[2] = TotCloPap_80inhib
            cell.pap[8] = TotureaPap_80inhib
        elif cell.inhib == 'NKCC2-70':
            cell.oi[0] = TotSodOI_70NKCCinhib
            cell.oi[1] = TotPotOI_70NKCCinhib
            cell.oi[2] = TotCloOI_70NKCCinhib		
            cell.oi[8] = TotureaOI_70NKCCinhib
            cell.pap[0] = TotSodPap_70NKCCinhib
            cell.pap[1] = TotPotPap_70NKCCinhib
            cell.pap[2] = TotCloPap_70NKCCinhib
            cell.pap[8] = TotureaPap_70NKCCinhib
        elif cell.inhib == 'NKCC2-100':
            cell.oi[0] = TotSodOI_100NKCCinhib
            cell.oi[1] = TotPotOI_100NKCCinhib
            cell.oi[2] = TotCloOI_100NKCCinhib		
            cell.oi[8] = TotureaOI_100NKCCinhib
            cell.pap[0] = TotSodPap_100NKCCinhib
            cell.pap[1] = TotPotPap_100NKCCinhib
            cell.pap[2] = TotCloPap_100NKCCinhib
            cell.pap[8] = TotureaPap_100NKCCinhib
        elif cell.inhib == 'SNB-70':
            cell.oi[0] = TotSodOI_70NKCCinhib
            cell.oi[1] = TotPotOI_70NKCCinhib
            cell.oi[2] = TotCloOI_70NKCCinhib		
            cell.oi[8] = TotureaOI_70NKCCinhib
            cell.pap[0] = TotSodPap_70NKCCinhib
            cell.pap[1] = TotPotPap_70NKCCinhib
            cell.pap[2] = TotCloPap_70NKCCinhib
            cell.pap[8] = TotureaPap_70NKCCinhib
        elif cell.inhib == 'SNB-100':
            cell.oi[0] = TotSodOI_100NKCCinhib
            cell.oi[1] = TotPotOI_100NKCCinhib
            cell.oi[2] = TotCloOI_100NKCCinhib		
            cell.oi[8] = TotureaOI_100NKCCinhib
            cell.pap[0] = TotSodPap_100NKCCinhib
            cell.pap[1] = TotPotPap_100NKCCinhib
            cell.pap[2] = TotCloPap_100NKCCinhib
            cell.pap[8] = TotureaPap_100NKCCinhib
            
    if cell.segment=='cTAL' or cell.segment == 'MD' or cell.segment=='DCT' or cell.segment=='PT' or cell.segment == 'CNT' or cell.segment == 'CCD':
        if cell.species == 'rat':
            pHplasma = 7.323
            phpap = 7.0
        elif cell.species == 'mou':
            pHplasma = 7.323
            phpap = 7.0
        elif cell.species == 'hum':
            pHplasma = 7.4
            phpap = 7.3
        else:
            print('cell.species: ' + str(cell.species))
            raise Exception('what is species?')

        if cell.unx == 'Y':
            pHplasma = 7.4-0.1
            phpap = 7.3-0.1
        
        cell.pH[5] = pHplasma

        facbic = np.exp(np.log(10)*(cell.pH[5]-pKHCO3))
        facpho = np.exp(np.log(10)*(cell.pH[5]-pKHPO4))
        facamm = np.exp(np.log(10)*(cell.pH[5]-pKNH3))
        fachco2 = np.exp(np.log(10)*(cell.pH[5]-pKHCO2))

        cell.conc[0,5] = cell.cm[0]
        cell.conc[1,5] = cell.cm[1]
        cell.conc[2,5] = cell.cm[2]
        cell.conc[3,5] = cell.cm[3]
        cell.conc[4,5] = cell.cm[4]
        cell.conc[5,5] = cell.cm[5]
             
        cell.conc[4,5] = cell.conc[5,5]/342.07
        cell.conc[3,5] = cell.conc[4,5]*10**(cell.pH[5]-pKHCO3)
        cell.conc[6,5] = cell.cm[6]*facpho/(1+facpho)
        cell.conc[7,5] = cell.cm[6]/(1+facpho)

        cell.conc[8,5] = cell.cm[8]
        
        if cell.segment=='cTAL':
            TotAmmCT = 1.0
            TotAmmCM = 1.5
            pos = (1*(cell.total-i))/(1*cell.total)
        elif cell.segment=='MD':
            TotAmmCT = 1.0
            TotAmmCM = 1.5
            pos = 1/200
        elif cell.segment=='DCT':
            TotAmmCT = 0.1
            if cell.type !='sup':
                TotAmmCM = 2.9
            else:
                TotAmmCM = 1.0
            pos = (1.0*(cell.total-i))/(1.0*cell.total)
        elif cell.segment =='PT':
            TotAmmCT = 0.1
            TotAmmCM = 1.5
            if cell.species == 'hum':
                total = cell.total*0.9
            elif cell.species == 'rat':
                total = cell.total*0.88
            elif cell.species == 'mou':
                total = cell.total*0.88
            else:
                print('cell.species: ' + str(cell.species))
                raise Exception('what is species?')
            pos = i/total
        elif cell.segment == 'CNT':
            TotAmmCT = 0.1
            TotAmmCM = 1.5
            pos = (1*(cell.total-i))/(1*cell.total)
        elif cell.segment == 'CCD':
            TotAmmCT = 0.1
            TotAmmCM = 1.5
            pos = (1*i)/(1*cell.total)            

        if cell.segment == 'CNT':
            Ammtotz = TotAmmCT
        else:
            Ammtotz = TotAmmCT+(TotAmmCM-TotAmmCT)*pos
                
        cell.conc[9,5] = Ammtotz*facamm/(1+facamm)
        cell.conc[10,5] = Ammtotz/(1+facamm)
        cell.conc[2,5] = cell.cm[2]+(Ammtotz/(1+facamm))

        cell.conc[11,5] = np.exp(-np.log(10)*cell.pH[5])*1000.0

        cell.conc[12,5] = cell.cm[12]*fachco2/(1+fachco2)
        cell.conc[13,5] = cell.cm[12]/(1+fachco2)
        cell.conc[14,5] = cell.cm[14]      

        elecS = 0
        for j in range(NS):
            elecS = elecS+zval[j]*cell.conc[j,5]

        cell.conc[2,5] = cell.conc[2,5]+elecS
        if cell.sex=='female' and cell.species == 'rat':
            female_conc(cell, i)

    elif cell.segment=='mTAL' or cell.segment=='SDL' or cell.segment=='OMCD':

        if cell.segment == 'mTAL':
            pos = (1*(cell.total-i))/(1*cell.total)
        elif cell.segment == 'SDL':
            pos = 0.3+0.7*i/(1*cell.total)
        elif cell.segment == 'OMCD':
            pos = (1*i)/(1*cell.total)
        if cell.species == 'rat':
            pHplasma = 7.323
            phpap = 7.0
        elif cell.species == 'mou':
            pHplasma = 7.323
            phpap = 7.0
        elif cell.species == 'hum':
            pHplasma = 7.4
            phpap = 7.3
        else:
            print('cell.species: ' + str(cell.species))
            raise Exception('what is species?')

        if cell.unx == 'Y':
            pHplasma = 7.4-0.1
            phpap = 7.3-0.1

        cell.pH[5] = pHplasma-(pHplasma-phpap)*2.0/7.0*pos

        facbic = np.exp(np.log(10)*(cell.pH[5]-pKHCO3))
        facpho = np.exp(np.log(10)*(cell.pH[5]-pKHPO4))
        facamm = np.exp(np.log(10)*(cell.pH[5]-pKNH3))
        fachco2 = np.exp(np.log(10)*(cell.pH[5]-pKHCO2))

        cell.conc[0,5] = cell.cm[0]+(cell.oi[0]-cell.cm[0])*pos
        cell.conc[1,5] = cell.cm[1]+(cell.oi[1]-cell.cm[1])*pos
        cell.conc[2,5] = cell.cm[2]+(cell.oi[2]-cell.cm[2])*pos
        cell.conc[3,5] = cell.cm[3]+(cell.oi[3]-cell.cm[3])*pos
        cell.conc[4,5] = cell.cm[4]+(cell.oi[4]-cell.cm[4])*pos
        cell.conc[5,5] = cell.cm[5]+(cell.oi[5]-cell.cm[5])*pos

        cell.conc[4,5] = cell.conc[5,5]/342.07
        cell.conc[3,5] = cell.conc[4,5]*10**(cell.pH[5]-pKHCO3)
        cell.conc[5,5] = cell.conc[4,5]*342.07

        Phototz = (cell.cm[6]+(cell.oi[6]-cell.cm[6])*pos)

        cell.conc[6,5] = Phototz*facpho/(1+facpho)
        cell.conc[7,5] = Phototz/(1+facpho)
        cell.conc[8,5] = cell.cm[8]+(cell.oi[8]-cell.cm[8])*pos

        TotAmmCM = 1.5
        Ammtotz = TotAmmCM+(cell.oi[9]-TotAmmCM)*pos
        cell.conc[9,5] = Ammtotz*facamm/(1+facamm)
        cell.conc[10,5] = Ammtotz/(1+facamm)
        cell.conc[11,5] = np.exp(-np.log(10)*cell.pH[5])*1000

        Hco2totz = cell.cm[12]+(cell.oi[12]-cell.cm[12])*pos
        cell.conc[12,5] = Hco2totz*fachco2/(1+fachco2)
        cell.conc[13,5] = Hco2totz/(1+fachco2)
        cell.conc[14,5] = cell.cm[14]+(cell.oi[14]-cell.cm[14])*pos

        elecS = 0.0
        for j in range(NS):
            elecS = elecS+zval[j]*cell.conc[j,5]

        cell.conc[2,5] = cell.conc[2,5]+elecS
 
        #  Concentrations of K, Cl are lower in female rat.
        if cell.sex=='female' and cell.species=='rat':
            female_conc(cell, i)
    
    elif cell.segment=='S3':
        if i == 1:
            if cell.species == 'rat':
                pHplasma = 7.323
                phpap = 7.0
            elif cell.species == 'mou':
                pHplasma = 7.323
                phpap = 7.0
            elif cell.species == 'hum':
                pHplasma = 7.4
                phpap = 7.3
            else:
                print('cell.species: ' + str(cell.species))
                raise Exception('what is species?')

            if cell.unx == 'Y':
                pHplasma = 7.4-0.1
                phpap = 7.3-0.1

            cell.pH[5] = pHplasma

            facbic = np.exp(np.log(10)*(cell.pH[5]-pKHCO3))
            facpho = np.exp(np.log(10)*(cell.pH[5]-pKHPO4))
            facamm = np.exp(np.log(10)*(cell.pH[5]-pKNH3))
            fachco2 = np.exp(np.log(10)*(cell.pH[5]-pKHCO2))

            cell.conc[0,5] = cell.cm[0]
            cell.conc[1,5] = cell.cm[1]
            cell.conc[2,5] = cell.cm[2]
            cell.conc[3,5] = cell.cm[3]
            cell.conc[4,5] = cell.cm[4]
            cell.conc[5,5] = cell.cm[5]

            cell.conc[4,5] = cell.conc[5,5]/342.07
            cell.conc[3,5] = cell.conc[4,5]*10**(cell.pH[5]-pKHCO3)
            cell.conc[6,5] = cell.cm[6]*facpho/(1+facpho)
            cell.conc[7,5] = cell.cm[6]/(1+facpho)

            cell.conc[8,5] = cell.cm[8]

            TotAmmCT = 0.1
            TotAmmCM = 1.5
            pos = 1

            Ammtotz = TotAmmCT+(TotAmmCM-TotAmmCT)*pos
            cell.conc[9,5] = Ammtotz*facamm/(1+facamm)
        
            cell.conc[10,5] = Ammtotz/(1+facamm)
            cell.conc[2,5] = cell.cm[2]+(Ammtotz/(1+facamm))

            cell.conc[11,5] = np.exp(-np.log(10)*cell.pH[5])*1000

            cell.conc[12,5] = cell.cm[12]*fachco2/(1+fachco2)
            cell.conc[13,5] = cell.cm[12]/(1+fachco2)
            cell.conc[14,5] = cell.cm[14]

            elecS = 0
            for j in range(NS):
                elecS = elecS+zval[j]*cell.conc[j,5]

            cell.conc[2,5] = cell.conc[2,5]+elecS
            
        else:
            if cell.species == 'hum':
                NPT=0.9*cell.total
                pHplasma = 7.4
                phpap = 7.3

            elif cell.species == 'rat':
                NPT=0.88*cell.total
                pHplasma = 7.323
                phpap = 7.0

            elif cell.species == 'mou':
                NPT=0.88*cell.total
                pHplasma = 7.323
                phpap = 7.0

            else:
                print('cell.species: ' + str(cell.species))
                raise Exception('what is species?')

            if cell.unx == 'Y':
                pHplasma = 7.4-0.1
                phpap = 7.3-0.1

            pos = (0.3*(i))/(1*cell.total-NPT)
            cell.pH[5] = pHplasma-(pHplasma-phpap)*2/7*pos

            facbic = np.exp(np.log(10)*(cell.pH[5]-pKHCO3))
            facpho = np.exp(np.log(10)*(cell.pH[5]-pKHPO4))
            facamm = np.exp(np.log(10)*(cell.pH[5]-pKNH3))
            fachco2 = np.exp(np.log(10)*(cell.pH[5]-pKHCO2))

            cell.conc[0,5] = cell.cm[0]+(cell.oi[0]-cell.cm[0])*pos
            cell.conc[1,5] = cell.cm[1]+(cell.oi[1]-cell.cm[1])*pos
            cell.conc[2,5] = cell.cm[2]+(cell.oi[2]-cell.cm[2])*pos
            cell.conc[3,5] = cell.cm[3]+(cell.oi[3]-cell.cm[3])*pos
            cell.conc[4,5] = cell.cm[4]+(cell.oi[4]-cell.cm[4])*pos
            cell.conc[5,5] = cell.cm[5]+(cell.oi[5]-cell.cm[5])*pos

            cell.conc[4,5] = cell.conc[5,5]/342.07
            cell.conc[3,5] = cell.conc[4,5]*10**(cell.pH[5]-pKHCO3)
            cell.conc[5,5] = cell.conc[4,5]*342.07

            Phototz = cell.cm[6]+(cell.oi[6]-cell.cm[6])*pos

            cell.conc[6,5] = Phototz*facpho/(1+facpho)
            cell.conc[7,5] = Phototz/(1+facpho)
            cell.conc[8,5] = cell.cm[8]+(cell.oi[8]-cell.cm[8])*pos

            TotAmmCM = 1.5
            Ammtotz = TotAmmCM+(cell.oi[9]-TotAmmCM)*pos
            cell.conc[9,5] = Ammtotz*facamm/(1+facamm)
            cell.conc[10,5] = Ammtotz/(1+facamm)
            cell.conc[11,5] = np.exp(-np.log(10)*cell.pH[5])*1000

            Hco2totz = cell.cm[12]+(cell.oi[12]-cell.cm[12])*pos
            cell.conc[12,5] = Hco2totz*fachco2/(1+fachco2)
            cell.conc[13,5] = Hco2totz/(1+fachco2)
            cell.conc[14,5] = cell.cm[14]+(cell.oi[14]-cell.cm[14])*pos

            elecS = 0.0
            for j in range(NS):
                elecS = elecS+zval[j]*cell.conc[j,5]

            cell.conc[2,5] = cell.conc[2,5]+elecS
            
        #  Concentrations of K, Cl are lower in female rat.
        if cell.sex=='female' and cell.species=='rat':
            female_conc(cell, i)

    elif cell.segment == 'IMCD' or cell.segment == 'LDL' or cell.segment == 'LAL':
        if cell.type == 'jux1':
            looplen = 0.2
        elif cell.type == 'jux2':
            looplen = 0.4
        elif cell.type == 'jux3':
            looplen = 0.6
        elif cell.type == 'jux4':
            looplen = 0.8
        elif cell.type == 'jux5':
            looplen = 1.0

        if cell.segment == 'IMCD':
            pos = i/cell.total
        elif cell.segment == 'LDL':
            pos = looplen*i/cell.total #looplen = 0.2 for jux1, 0.4 for jux2, 0.6 for jux3, 0.8 for jux4, 1.0 for jux5
        elif cell.segment == 'LAL':
            pos = looplen*(cell.total-i)/cell.total                
        
        if cell.species == 'rat':
            pHplasma = 7.323
            phpap = 7.0
        elif cell.species == 'mou':
            pHplasma = 7.323
            phpap = 7.0
        elif cell.species == 'hum':
            pHplasma = 7.4
            phpap = 7.3
        else:
            print('cell.species: ' + str(cell.species))
            raise Exception('what is species?')

        if cell.unx == 'Y':
            pHplasma = 7.4-0.1
            phpap = 7.3-0.1

        cell.ep[5]=-0.001e-3/EPref

        cell.pH[5] = pHplasma-(pHplasma-phpap)*((2/7)+(5/7)*pos)

        facbic = np.exp(np.log(10)*(cell.pH[5]-pKHCO3))
        facpho = np.exp(np.log(10)*(cell.pH[5]-pKHPO4))
        facamm = np.exp(np.log(10)*(cell.pH[5]-pKNH3))
        fachco2 = np.exp(np.log(10)*(cell.pH[5]-pKHCO2))

        cell.conc[0,5] = cell.oi[0]+(cell.pap[0]-cell.oi[0])*pos
        cell.conc[1,5] = cell.oi[1]+(cell.pap[1]-cell.oi[1])*pos
        cell.conc[2,5] = cell.oi[2]+(cell.pap[2]-cell.oi[2])*pos
        cell.conc[3,5] = cell.oi[3]+(cell.pap[3]-cell.oi[3])*pos
        cell.conc[4,5] = cell.oi[4]+(cell.pap[4]-cell.oi[4])*pos
        cell.conc[5,5] = cell.oi[5]+(cell.pap[5]-cell.oi[5])*pos

        cell.conc[4,5] = cell.conc[5,5]/342.07
        cell.conc[3,5] = cell.conc[4,5]*10.0**(cell.pH[5]-pKHCO3)
        cell.conc[5,5] = cell.conc[4,5]*342.07

        Phototz = cell.oi[6]+(cell.pap[6]-cell.oi[6])*pos

        cell.conc[6,5] = Phototz*facpho/(1+facpho)
        cell.conc[7,5] = Phototz/(1+facpho)
        cell.conc[8,5] = cell.oi[8]+(cell.pap[8]-cell.oi[8])*pos

        Ammtotz = cell.oi[9]+(cell.pap[9]-cell.oi[9])*pos
        cell.conc[9,5] = Ammtotz*facamm/(1+facamm)
        cell.conc[10,5] = Ammtotz/(1+facamm)
        cell.conc[11,5] = np.exp(-np.log(10)*cell.pH[5])*1000

        Hco2totz = cell.oi[12]+(cell.pap[12]-cell.oi[12])*pos
        cell.conc[12,5] = Hco2totz*fachco2/(1+fachco2)
        cell.conc[13,5] = Hco2totz/(1+fachco2)
        cell.conc[14,5] = cell.oi[14]+(cell.pap[14]-cell.oi[14])*pos

        elecS = 0.0
        for j in range(NS):
            elecS = elecS+zval[j]*cell.conc[j,5]

        cell.conc[2,5] = cell.conc[2,5]+elecS

        #  Concentrations of K, Cl are lower in female rat.
        if cell.sex=='female' and cell.species=='rat':
            female_conc(cell, i)

    else:
        cell.conc[:,5] = cell.conc[:,5]
