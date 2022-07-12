from compute_sup_jux_segment import compute_segment 
from driver import compute
from values import *
from defs import *
import os
import argparse
import multiprocessing
import output
import re

solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
cw=Vref*60e6

parser = argparse.ArgumentParser()
# required input
parser.add_argument('--sex',choices=['Male','Female'],required = True,type = str,help = 'Sex')
parser.add_argument('--species',choices=['human','rat','mouse'],required = True,type = str, help = 'Human, Rat, or Mouse model')
parser.add_argument('--type',choices = ['superficial','multiple'],required = True,type=str,help='superficial nephron or multiple nephrons?')

# diabetic options
parser.add_argument('--diabetes',choices = ['Severe','Moderate'],default='Non',type=str,help='diabete status (Severe/Moderate)')
parser.add_argument('--inhibition',choices=['ACE','SGLT2','NHE3-50','NHE3-80','NKCC2-70','NKCC2-100','NCC-70','NCC-100','ENaC-70','ENaC-100','SNB-70','SNB-100'],default = None,type = str,help = 'any transporter inhibition?')
parser.add_argument('--unx',choices=['N','Y'],default = 'N',type = str,help = 'uninephrectomy status')

# pregnancy option
parser.add_argument('--pregnant', choices=['mid','late'], default='non', type=str, help='pregnant female? (mid/late)')
parser.add_argument('--HT',choices=['N','Y'],default = 'N',type = str,help = 'hypertension?')

args = parser.parse_args()
sex = args.sex
species = args.species
sup_or_multi = args.type
diabete = args.diabetes
inhib = args.inhibition
unx = args.unx
preg = args.pregnant
HT = args.HT

if diabete != 'Non':
    if preg != 'non':
        raise Exception('pregnant diabetic not done')
    if inhib != None:
        file_to_save = inhib+'_'+sex+'_'+species[0:3]+'_'+diabete+'_diab'+'_'+unx+'_unx'
    else:
        file_to_save = sex+'_'+species[0:3]+'_'+diabete+'_diab'+'_'+unx+'_unx'
elif preg != 'non':
    if sex == 'Male':
        raise Exception('pregnant only for female')
    if species[0:3] == 'hum' or species[0:3] == 'mou':
        raise Exception('pregnant model not set up for human or mouse yet')
    if HT == 'Y':
        file_to_save = preg+'pregnant_'+species[0:3]+'_HT'
    else:
        file_to_save = preg+'pregnant_'+species[0:3]
elif inhib != None:
    file_to_save = inhib+'_'+sex+'_'+species[0:3]
elif HT != 'N':
    file_to_save = sex + '_' + species[0:3]+'_HT'
else:
    file_to_save = sex + '_' + species[0:3] +'_normal'
    
if os.path.isdir(file_to_save) == False:
    os.makedirs(file_to_save)
    
if sup_or_multi == 'superficial':
    parts = ['sup']
else:
    parts = ['sup','jux1','jux2','jux3','jux4','jux5']

if os.path.isdir('outlets') == False:
    os.makedirs('outlets')
    
def multiprocessing_func(sup_or_jux):
    compute_segment(sup_or_jux, sex, species, sup_or_multi, diabete, inhib, unx, preg, HT, file_to_save)

if __name__ == '__main__':

    pool = multiprocessing.Pool()
    pool.map(multiprocessing_func,parts)
    pool.close()

    #========================================================
    # Cortical collecting duct
    #========================================================
    print('Collecting duct begin')
    print('CCD start')
    NCCD = 200
    if sex == 'Male':
        filename = './datafiles/CCDparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/CCDparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/CCDparams_F_'+species[0:3]+'.dat'
    ccd=compute(NCCD,filename,'Newton',diabete=diabete,species=species,sup_or_multi=sup_or_multi,inhibition = inhib,unx = unx, preg=preg, HT=HT)

    Scaletorq = np.ones(NCCD)

    output.output_segment_results(ccd,"",Scaletorq,file_to_save,NCCD)

    print('CCD finished.')
    print('\n')
    
    #========================================================
    # Outer medullary collecting duct
    #========================================================
    print('OMCD start')
    NOMCD = 200
    if sex == 'Male':
        filename = './datafiles/OMCDparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/OMCDparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/OMCDparams_F_'+species[0:3]+'.dat'
    if ccd[0].sex == 'male':
        omcd=compute(NOMCD,filename,'Newton',diabete=diabete,species=species,sup_or_multi=sup_or_multi,inhibition=inhib,unx=unx, preg=preg, HT=HT)
    elif ccd[0].sex == 'female':
        omcd=compute(NOMCD,filename,'Newton',diabete=diabete,species=species,sup_or_multi=sup_or_multi,inhibition=inhib,unx=unx, preg=preg, HT=HT)

    Scaletorq = np.ones(NOMCD)

    output.output_segment_results(omcd,"",Scaletorq,file_to_save,NOMCD)
    
    print('OMCD finished.')
    print('\n')

    #========================================================
    # Inner medullary collecting duct
    #========================================================
    print('IMCD start')
    NIMCD = 200
    if sex == 'Male':
        filename = './datafiles/IMCDparams_M_'+species[0:3]+'.dat'
    elif sex == 'Female':
        filename = './datafiles/IMCDparams_F_'+species[0:3]+'.dat'
    else:
        filename ='./datafiles/IMCDparams_F_'+species[0:3]+'.dat'
    imcd=compute(NIMCD,filename,'Newton',diabete=diabete,species=species,sup_or_multi=sup_or_multi,inhibition=inhib,unx=unx, preg=preg, HT=HT)

    Scaletorq = np.ones(NIMCD)

    jvol = output.output_segment_results(imcd,"",Scaletorq,file_to_save,NIMCD)
    
    print('IMCD finished.')
