from driver import compute
from values import *
from defs import *
from set_params import set_torq_params
import re
import electrochemical 
import water
import glucose
import cotransport
import NHE3
import ATPase
import NKCC
import KCC
import NCC
import ENaC
import Pendrin
import AE1
import NHE1
import flux
import output

def compute_segment(sup_or_jux,sex,species,sup_or_multi,diabete,inhib,unx,preg,HT,file_to_save):
    solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
    compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
    cw=Vref*60e6
    #========================================================
    # Proximal convolute tubule
    #========================================================
    print('%s PCT start' %(sup_or_jux))
    if species == 'human':
        NPT = 181
    elif species == 'rat':
        NPT = 176
    elif species == 'mouse':
        NPT = 176
    else:
        print(str(species))
        raise Exception('what is species?')

    if sex == 'Male':
        filename = './datafiles/PTparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/PTparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/PTparams_F_'+species[0:3]+'.dat'

    pt=compute(NPT,filename,'Broyden',sup_or_jux,diabete,species,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx, preg = preg, HT=HT)

    Scaletorq = np.zeros(NPT)
    
    for j in range(NPT):
        if pt[j].segment == 'PT':
            TS = 1.3
            scaleT = 1.0
        elif pt[j].segment == 'S3':
            TS = 1.3
            scaleT = 0.5

        #torque-modulated effects                
        PM=pt[j].pres[0]                
        Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(pt[j].species,pt[j].sex,pt[j].preg)
                                    
        if pt[j].species == 'rat':
            fac1 = 8.0*visc*(pt[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif pt[j].species == 'mou':
            fac1 = 8.0*visc*(pt[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif pt[j].species == 'hum':
            fac1 = 8.0*visc*(pt[j].volref[0]*Vref)*torqL/(Radref**2)
        else:
            print('pt.species: ' + str(pt[j].species))
            raise Exception('what is species?')
        fac2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
        TM0= fac1*fac2
            
        RMtorq = torqR*(1.0e0+torqvm*(PM - PbloodPT))

        factor1 = 8.0*visc*(pt[j].vol[0]*Vref)*torqL/(RMtorq**2) 
        factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
        Torque = factor1*factor2
            
        Scaletorq[j] = 1.0 + TS*scaleT*(Torque/TM0-1.0)

    output.output_segment_results(pt,sup_or_jux,Scaletorq,file_to_save,NPT)

    print('%s PCT finished.' %(sup_or_jux))
    print('\n')
    
    #========================================================
    # S3
    #========================================================
    print('%s S3 start' %(sup_or_jux))
    if species == 'human':
        NS3 = 20
    elif species == 'rat':
        NS3 = 25
    elif species == 'mouse':
        NS3 = 25
    else:
        print('cell.species: ' + str(species))
        raise Exception('what is species?')
    
    if sex == 'Male':
        filename = './datafiles/S3params_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/S3params_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/S3params_F_'+species[0:3]+'.dat'
    s3=compute(NS3,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx,preg = preg, HT = HT)

    Scaletorq = np.zeros(NS3)
    
    for j in range(NS3):
        if s3[j].segment == 'PT':
            TS = 1.3
            scaleT = 1.0
        elif s3[j].segment == 'S3':
            TS = 1.3
            scaleT = 0.5

        #torque-modulated effects
                
        PM=s3[j].pres[0]

        Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(s3[j].species,s3[j].sex,s3[j].preg)
                    
        if s3[j].species == 'rat':
            fac1 = 8.0*visc*(s3[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif s3[j].species == 'mou':
            fac1 = 8.0*visc*(s3[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif s3[j].species == 'hum':
            fac1 = 8.0*visc*(s3[j].volref[0]*Vref)*torqL/(Radref**2) 
        else:
            print('s3.species: ' + str(s3[j].species))
            raise Exception('what is species?')
        fac2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
        TM0= fac1*fac2
            
        RMtorq = torqR*(1.0e0+torqvm*(PM - PbloodPT))

        factor1 = 8.0*visc*(s3[j].vol[0]*Vref)*torqL/(RMtorq**2) 
        factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
        Torque = factor1*factor2
        
        Scaletorq[j] = 1.0 + TS*scaleT*(Torque/TM0-1.0)

    output.output_segment_results(s3,sup_or_jux,Scaletorq,file_to_save,NS3)

    print('%s S3 finished.' %(sup_or_jux))
    print('\n')
    
    #========================================================
    # Short descending limb
    #========================================================
    print('%s SDL start' %(sup_or_jux))
    NSDL = 200
    if species == 'human':
        method = 'Newton'
    elif species == 'rat':
        method = 'Broyden'
    elif species == 'mouse':
        method = 'Broyden'
    else:
        print('species: ' + str(species))
        raise Exception('what is species?')
        
    if sex == 'Male':
        filename = './datafiles/SDLparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/SDLparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/SDLparams_F_'+species[0:3]+'.dat'
    #sdl=compute(NSDL,filename,'Broyden',diabete)
    sdl=compute(NSDL,filename,method,sup_or_jux,diabete,species,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx, preg = preg, HT = HT)
    
    Scaletorq = np.ones(NSDL)
    output.output_segment_results(sdl,sup_or_jux,Scaletorq,file_to_save,NSDL)

    print('%s SDL finished.' %(sup_or_jux))
    print('\n')
    
    #========================================================
    # Long descending limb
    #========================================================
    if sup_or_jux != 'sup':
        print('%s LDL start' %(sup_or_jux))
        NLDL = 200
        if sex == 'Male':
            filename = './datafiles/LDLparams_M_'+species[0:3]+'.dat'
        elif sex == 'Female':
            filename = './datafiles/LDLparams_F_'+species[0:3]+'.dat'
        else:
            filename ='./datafiles/LDLparams_F_'+species[0:3]+'.dat'
        ldl=compute(NLDL,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx, preg = preg, HT = HT)

        Scaletorq = np.ones(NLDL)
        output.output_segment_results(ldl,sup_or_jux,Scaletorq,file_to_save,NLDL)
        
        print('%s LDL finished.' %(sup_or_jux))
        print('\n')
        
    #========================================================
    # Long ascending limb
    #========================================================
        print('%s LAL start' %(sup_or_jux))
        NLAL = 200
        if sex == 'Male':
            filename = './datafiles/LALparams_M_rat.dat'
        elif sex == 'Female':
            filename = './datafiles/LALparams_F_rat.dat'
        else:
            filename ='./datafiles/LALparams_F_rat.dat'
        lal=compute(NLAL,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx, preg = preg, HT = HT)

        Scaletorq = np.ones(NLAL)
        output.output_segment_results(lal,sup_or_jux,Scaletorq,file_to_save,NLAL)

        print('%s LAL finished.' %(sup_or_jux))
        print('\n')

    #========================================================
    # Medulla thick ascending limb
    #========================================================
    print('%s mTAL start' %(sup_or_jux))
    NmTAL = 200
    if sex == 'Male':
        filename = './datafiles/mTALparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/mTALparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/mTALparams_F_'+species[0:3]+'.dat'
    mtal=compute(NmTAL,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi,inhib,unx = unx, preg = preg, HT = HT)

    Scaletorq = np.ones(NmTAL)
    
    output.output_segment_results(mtal,sup_or_jux,Scaletorq,file_to_save,NmTAL)

    print('%s mTAL finished.' %(sup_or_jux))
    print('\n')

    #========================================================
    # Cortex thick ascending limb
    #========================================================
    print('%s cTAL start' %(sup_or_jux))
    NcTAL = 200
    if sex == 'Male':
        filename = './datafiles/cTALparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/cTALparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/cTALparams_F_'+species[0:3]+'.dat'
    ctal=compute(NcTAL,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi,inhib,unx = unx, preg = preg, HT = HT)

    Scaletorq = np.ones(NcTAL)
    
    output.output_segment_results(ctal,sup_or_jux,Scaletorq,file_to_save,NcTAL)

    print('%s cTAL finished.' %(sup_or_jux))
    print('\n')
    
    #========================================================
    # Distal convoluted tubule
    #========================================================
    print('%s DCT start' %(sup_or_jux))
    NDCT = 200
    if sex == 'Male':
        filename = './datafiles/DCTparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/DCTparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/DCTparams_F_'+species[0:3]+'.dat'
    dct=compute(NDCT,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi,inhib,unx = unx, preg = preg, HT = HT)

    Scaletorq = np.ones(NDCT)
    
    output.output_segment_results(dct,sup_or_jux,Scaletorq,file_to_save,NDCT)

    print('%s DCT finished.'%(sup_or_jux))
    print('\n')
    
    #=======================================================
    # Connecting tubule
    #=======================================================
    print('%s CNT start' %(sup_or_jux))
    NCNT = 200
    if sex == 'Male':
        filename = './datafiles/CNTparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/CNTparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/CNTparams_F_'+species[0:3]+'.dat'
    cnt=compute(NCNT,filename,'Newton',sup_or_jux,diabete,species,sup_or_multi,inhib,unx = unx, preg = preg, HT = HT)

    Scaletorq = np.ones(NCNT)
    
    output.output_segment_results(cnt,sup_or_jux,Scaletorq,file_to_save,NCNT)

    print('%s CNT finished.'%(sup_or_jux))
    print(sup_or_jux+' finished.')
    print('\n')
