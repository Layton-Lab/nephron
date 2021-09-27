from driver import compute
from values import *
from defs import *
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
import os
import sys
import argparse
import multiprocessing
from set_params import set_torq_params
import output

solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
cw=Vref*60e6

parser = argparse.ArgumentParser()
# required input
parser.add_argument('--sex',choices=['Male','Female'],required = True,type = str,help = 'Sex')
parser.add_argument('--species',choices=['human','rat'],required = True,type = str, help = 'Human model or Rat model')
parser.add_argument('--type',choices = ['superficial','multiple'],required = True,type=str,help='superficial nephron or multiple nephrons?')
parser.add_argument('--file2save', required = True, type = str, help = 'where to save?')

# diabetic options
parser.add_argument('--diabetes',choices = ['Severe','Moderate'],default='Non',type=str,help='diabete status (Severe/Moderate)')
parser.add_argument('--inhibition',choices=['ACE','SGLT2','NHE3-50','NHE3-80','NKCC2-70','NKCC2-100','NCC-70','NCC-100','ENaC-70','ENaC-100','SNB-70','SNB-100'],default = None,type = str,help = 'any transporter inhibition?')
parser.add_argument('--unx',choices=['N','Y'],default = 'N',type = str,help = 'uninephrectomy status')
# pregnancy option
parser.add_argument('--pregnant', choices=['mid','late'], default='non', type=str, help='pregnant female? (mid/late)')

args = parser.parse_args()
sex = args.sex
humOrrat = args.species
sup_or_multi = args.type
diabete = args.diabetes
inhib = args.inhibition
unx = args.unx
preg = args.pregnant

if diabete != 'Non':
    if preg != 'non':
        raise Exception('pregnant diabetic not done')
    # if inhib != None:
    #     file_to_save = inhib+'_'+sex+'_'+humOrrat[0:3]+'_'+diabete+'_diab'+'_'+unx+'_unx'
    # else:
    #     file_to_save = sex+'_'+humOrrat[0:3]+'_'+diabete+'_diab'+'_'+unx+'_unx'
elif preg != 'non':
    if sex == 'Male':
        raise Exception('pregnant only for female')
    if humOrrat[0:3] == 'hum':
        raise Exception('pregnant model not set up for human yet')
    if inhib != None:
        raise Exception('pregnant model does not have inhibition set up yet')

    #file_to_save = preg+'pregnant_'+humOrrat[0:3]
# else:
#     file_to_save = sex + '_' + humOrrat[0:3] +'_normal'

file_to_save = args.file2save
    
if os.path.isdir(file_to_save) == False:
    os.makedirs(file_to_save)
    
if sup_or_multi == 'superficial':
    parts = ['sup']
else:
    parts = ['sup','jux1','jux2','jux3','jux4','jux5']

def compute_segmentPTS3(sup_or_jux,sex,humOrrat,sup_or_multi,diabete,inhib,unx,preg,file_to_save):
    solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
    compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
    cw=Vref*60e6
    #========================================================
    # Proximal convolute tubule
    #========================================================
    print('%s PCT start' %(sup_or_jux))
    if humOrrat == 'human':
        NPT = 181
    elif humOrrat == 'rat':
        NPT = 176
    if sex == 'Male':
        filename = './datafiles/PTparams_M_'+humOrrat[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/PTparams_F_'+humOrrat[0:3]+'.dat'
    else:
        filename ='./datafiles/PTparams_F_'+humOrrat[0:3]+'.dat'

    pt=compute(NPT,filename,'Broyden',sup_or_jux,diabete,humOrrat,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx, preg = preg)
    
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
        Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(pt[j].humOrrat,pt[j].sex,pt[j].preg)
                                    
        if pt[j].humOrrat == 'rat':
            fac1 = 8.0*visc*(pt[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif pt[j].humOrrat == 'mou':
            fac1 = 8.0*visc*(pt[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif pt[j].humOrrat == 'hum':
            fac1 = 8.0*visc*(pt[j].volref[0]*Vref)*torqL/(Radref**2)
        else:
            print('pt.humOrrat: ' + str(pt[j].humOrrat))
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
    if humOrrat == 'human':
        NS3 = 20
    elif humOrrat == 'rat':
        NS3 = 25
    if sex == 'Male':
        filename = './datafiles/S3params_M_'+humOrrat[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/S3params_F_'+humOrrat[0:3]+'.dat'
    else:
        filename ='./datafiles/S3params_F_'+humOrrat[0:3]+'.dat'
    s3=compute(NS3,filename,'Newton',sup_or_jux,diabete,humOrrat,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx,preg = preg)

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

        Radref,torqR,torqvm,PbloodPT,torqL,torqd = set_torq_params(s3[j].humOrrat,s3[j].sex,s3[j].preg)
                    
        if s3[j].humOrrat == 'rat':
            fac1 = 8.0*visc*(s3[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif s3[j].humOrrat == 'mou':
            fac1 = 8.0*visc*(s3[j].vol_init[0]*Vref)*torqL/(Radref**2)
        elif s3[j].humOrrat == 'hum':
            fac1 = 8.0*visc*(s3[j].volref[0]*Vref)*torqL/(Radref**2) 
        else:
            print('s3.humOrrat: ' + str(s3[j].humOrrat))
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
#=============================
# end compute_segmentPTS3
#=============================

def multiprocessing_funcPTS3(sup_or_jux):
    compute_segmentPTS3(sup_or_jux, sex, humOrrat, sup_or_multi, diabete, inhib, unx, preg, file_to_save)

if __name__ == '__main__':
    pool = multiprocessing.Pool()
    pool.map(multiprocessing_funcPTS3, parts)
    pool.close()
