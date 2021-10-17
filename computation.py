# This file is used to check individual segment. Type 'python3 computation.py' in the terminal to run.
# requires outlet files from previous simulation
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
import argparse
import re

solute = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu']
compart = ['Lumen','Cell','ICA','ICB','LIS','Bath']
cw=Vref*60e6

parser=argparse.ArgumentParser()
# required input
parser.add_argument('--sex',choices=['Male','Female'],required = True,type = str,help = 'Sex')
parser.add_argument('--species',choices=['human','rat', 'mouse'],required = True,type = str, help = 'Human model or Rat model')
parser.add_argument('--type',choices = ['superficial','multiple'],required = True,type=str,help='superficial nephron or multiple nephrons?')
parser.add_argument('--segment', choices = ['PT','S3','SDL', 'LDL', 'LAL', 'mTAL','cTAL','DCT', 'CNT', 'CCD', 'OMCD', 'IMCD'], required=True, type=str, help = 'choose segment')
parser.add_argument('--savefile', required=True, type=str, help = 'where to save?')
# optional input
parser.add_argument('--suporjux', choices=['sup','jux1','jux2','jux3','jux4','jux5', ''], default='', type=str, help = 'which nephron type? (sup/jux1/jux2/etc), '' is for collecting duct')
# diabetic options
parser.add_argument('--diabetes',choices = ['Severe','Moderate'],default='Non',type=str,help='diabete status (Severe/Moderate)')
parser.add_argument('--inhibition',choices=['ACE','SGLT2','NHE3-50','NHE3-80','NKCC2-70','NKCC2-100','NCC-70','NCC-100','ENaC-70','ENaC-100','SNB-70','SNB-100'],default = None,type = str,help = 'any transporter inhibition?')
parser.add_argument('--unx',choices=['N','Y'],default = 'N',type = str,help = 'uninephrectomy status')
# pregnancy option
parser.add_argument('--pregnant', choices=['mid','late'], default='non', type=str, help='pregnant female? (mid/late)')

args=parser.parse_args()

sex = args.sex
species = args.species
sup_or_multi = args.type
segment = args.segment
sup_or_jux = args.suporjux

if sup_or_jux == '':
	if sup_or_multi == 'multiple':
		if segment[-2:] != 'CD':
			print('segment: ' + segment)
			raise Exception('sup or jux required for multiple nephron segments')
	elif sup_or_multi == 'superficial':
		if segment[-2:] != 'CD':
			sup_or_jux = 'sup'
	else:
		print('sup_or_multi: ' + sup_or_multi)
		raise Exception('what is this sup_or_multi? ' + sup_or_multi)

diabete = args.diabetes
inhib = args.inhibition
unx = args.unx

preg = args.pregnant

file_to_save = args.savefile
if os.path.isdir(file_to_save) == False:
    os.makedirs(file_to_save)

if sex == 'Male':
	filename='./datafiles/'+segment+'params_M_'+species[0:3]+'.dat'
elif sex == 'Female':
	filename='./datafiles/'+segment+'params_F_'+species[0:3]+'.dat'
else:
	print('sex: ' + sex)
	raise Exception('must be male or female')

file = open(filename, 'r')
line = file.readline()
while (line):
	line = line.replace('\t', ' ')
	terms = line.split(' ')
	if ((line[0][0] != '#') and ('total' == terms[0])):
		first_space_pos = line.index(' ')
		num = re.findall(r'-?\d+\.?\d*[Ee]?[+-]?\d*', line[first_space_pos:len(line)])
		N = float(num[0])
		break
	else:
		line = file.readline()
file.close()

if segment == 'PT':
    N = 7*N / 8
elif segment == 'S3':
    N = N / 8

N = int(N)
method = 'Newton'
cell=compute(N,filename,method,sup_or_jux,diabete=diabete,species=species,sup_or_multi = sup_or_multi,inhibition=inhib,unx = unx, preg=preg)
if sup_or_jux != '':
	sup_or_jux = '_' + sup_or_jux

if cell[0].preg != 'non':
	sex_or_preg = cell[0].preg + 'pregnant'
else:
	sex_or_preg = cell[0].sex

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Lumen_potential'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].ep[0])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Cell_potential'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].ep[1])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_LIS_potential'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].ep[4])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_potential_gradient_Lumen_Cell'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].ep[0]-cell[j-1].ep[1])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_potential_gradient_Lumen_LIS'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].ep[0]-cell[j-1].ep[4])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_conc_gradient_Lumen_LIS'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].conc[0,0]-cell[j-1].conc[0,4])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_conc_gradient_Cell_LIS'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].conc[0,1]-cell[j-1].conc[0,4])+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_pressure_in_Lumen'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].pres[0])+'\n')
file.close()

for i in range(NS):
	file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_con_of_'+solute[i]+'_in_Lumen'+sup_or_jux+'.txt','w')
	for j in range(1,N):
		file.write(str(cell[j-1].conc[i,0])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_con_of_'+solute[i]+'_in_Cell'+sup_or_jux+'.txt','w')
	for j in range(1,N):
		file.write(str(cell[j-1].conc[i,1])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_con_of_'+solute[i]+'_in_LIS'+sup_or_jux+'.txt','w')
	for j in range(1,N):
		file.write(str(cell[j-1].conc[i,4])+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_con_of_'+solute[i]+'_in_Bath'+sup_or_jux+'.txt','w')
	for j in range(1,N):
		file.write(str(cell[j-1].conc[i,5])+'\n')
	file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_water_volume_in_Lumen'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].vol[0]*cw)+'\n')
file.close()
file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_water_volume_in_Cell'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(cell[j-1].vol[1]*cw)+'\n')
file.close()

file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_pH_in_Lumen'+sup_or_jux+'.txt','w')
for j in range(1,N):
	file.write(str(-np.log(abs(cell[j-1].conc[11,0])/1000)/np.log(10))+'\n')
file.close()

for i in range(NS):
	file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_flow_of_'+solute[i]+'_in_Lumen'+sup_or_jux+'.txt','w')
	for j in range(1,N):
		file.write(str(cell[j-1].conc[i,0]*cell[j-1].vol[0]*cw)+'\n')
	file.close()
for i in range(NS):
	file=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_flow_of_'+solute[i]+'_in_Cell'+sup_or_jux+'.txt','w')
	for j in range(1,N):
		file.write(str(cell[j-1].conc[i,1]*cell[j-1].vol[1]*cw)+'\n')
	file.close()

file_lumen = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_osmolality_in_Lumen'+sup_or_jux+'.txt','w')
file_cell = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_osmolality_in_Cell'+sup_or_jux+'.txt','w')
file_lis = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_osmolality_in_LIS'+sup_or_jux+'.txt','w')
file_bath = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_osmolality_in_Bath'+sup_or_jux+'.txt','w')
for j in range(N):
	osm_l = 0
	osm_c = 0
	osm_lis = 0
	osm_b = 0
	for i in range(NS):
		osm_l = osm_l +cell[j].conc[i,0]
		osm_c = osm_c +cell[j].conc[i,1]
		osm_lis = osm_lis+cell[j].conc[i,4]
		osm_b = osm_b +cell[j].conc[i,5]

	file_lumen.write(str(osm_l)+'\n')
	file_cell.write(str(osm_c)+'\n')
	file_lis.write(str(osm_lis)+'\n')
	file_bath.write(str(osm_b)+'\n')
file_lumen.close()
file_cell.close()
file_lis.close()
file_bath.close()

for j in range(1,N):
	cell[j].area[4][5] = 0.02*max(cell[j].vol[4]/cell[j].volref[4],1.0)
	cell[j].area[5][4] = cell[j].area[4][5]

	jvol = water.compute_water_fluxes(cell[j])
	jsol,delmu = electrochemical.compute_ecd_fluxes(cell[j],jvol)
	file_na_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_Na_apical_ecd.txt','a')
	file_na_apical.write(str(jsol[0,0,1])+'\n')
	file_k_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_K_apical_ecd.txt','a')
	file_k_apical.write(str(jsol[1,0,1])+'\n')
	file_nh4_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_NH4_apical_ecd.txt','a')
	file_nh4_apical.write(str(jsol[10,0,1])+'\n')
	file_h_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_H_apical_ecd.txt','a')
	file_h_apical.write(str(jsol[11,0,1])+'\n')
	file_cl_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_Cl_apical_ecd.txt','a')
	file_cl_apical.write(str(jsol[2,0,1])+'\n')
	file_hco3_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_HCO3_apical_ecd.txt','a')
	file_hco3_apical.write(str(jsol[3,0,1])+'\n')
	file_hpo4_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_HPO4_apical_ecd.txt','a')
	file_hpo4_apical.write(str(jsol[6,0,1])+'\n')
	file_h2po4_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_H2PO4_apical_ecd.txt','a')
	file_h2po4_apical.write(str(jsol[7,0,1])+'\n')
	file_hco2_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_HCO2_apical_ecd.txt','a')
	file_hco2_apical.write(str(jsol[12,0,1])+'\n')
	file_nh3_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_NH3_apical_ecd.txt','a')
	file_nh3_apical.write(str(jsol[9,0,1])+'\n')
	file_glu_apical=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_glu_apical_ecd.txt','a')
	file_glu_apical.write(str(jsol[14,0,1])+'\n')

	file_na_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_Na_para_ecd.txt','a')
	file_na_para.write(str(jsol[0,0,4])+'\n')
	file_k_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_K_para_ecd.txt','a')
	file_k_para.write(str(jsol[1,0,4])+'\n')
	file_nh4_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_NH4_para_ecd.txt','a')
	file_nh4_para.write(str(jsol[10,0,4])+'\n')
	file_h_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_H_para_ecd.txt','a')
	file_h_para.write(str(jsol[11,0,4])+'\n')
	file_cl_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_Cl_para_ecd.txt','a')
	file_cl_para.write(str(jsol[2,0,4])+'\n')
	file_hco3_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_HCO3_para_ecd.txt','a')
	file_hco3_para.write(str(jsol[3,0,4])+'\n')
	file_hpo4_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_HPO4_para_ecd.txt','a')
	file_hpo4_para.write(str(jsol[6,0,4])+'\n')
	file_h2po4_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_H2PO4_para_ecd.txt','a')
	file_h2po4_para.write(str(jsol[7,0,4])+'\n')
	file_hco2_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_HCO2_para_ecd.txt','a')
	file_hco2_para.write(str(jsol[12,0,4])+'\n')
	file_nh3_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_NH3_para_ecd.txt','a')
	file_nh3_para.write(str(jsol[9,0,4])+'\n')
	file_glu_para=open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_glu_para_ecd.txt','a')
	file_glu_para.write(str(jsol[14,0,4])+'\n')
	for i in range(len(cell[j].trans)):
		transporter_type = cell[j].trans[i].type
		memb_id = cell[j].trans[i].membrane_id

		if transporter_type == 'SGLT1':
			solute_id,fluxs = glucose.sglt1(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'SGLT2':
			solute_id,fluxs = glucose.sglt2(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'GLUT1':
			solute_id,fluxs=glucose.glut1(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')
		elif transporter_type == 'GLUT2':
			solute_id,fluxs=glucose.glut2(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len([solute_id])):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs)+'\n')			
		elif transporter_type == 'NHE3':
			solute_id,fluxs=NHE3.nhe3(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NaKATPase':
			solute_id,fluxs=ATPase.nakatpase(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs[k])+'\n')

		elif transporter_type == 'HATPase':
			solute_id,fluxs=ATPase.hatpase(cell[j],cell[j].ep,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2A':
			solute_id,fluxs=NKCC.nkcc2(cell[j],memb_id,cell[j].trans[i].act,cell[j].area,'A')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2B':
			solute_id,fluxs=NKCC.nkcc2(cell[j],memb_id,cell[j].trans[i].act,cell[j].area,'B')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC2F':
			solute_id,fluxs=NKCC.nkcc2(cell[j],memb_id,cell[j].trans[i].act,cell[j].area,'F')
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')       
		elif transporter_type == 'KCC4':
			solute_id,fluxs=KCC.kcc4(cell[j].conc,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+str(memb_id[0])+str(memb_id[1])+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'ENaC':
			solute_id,fluxs=ENaC.ENaC(cell[j],j,memb_id,cell[j].trans[i].act,cell[j].area,jvol)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NCC':
			solute_id,fluxs=NCC.NCC(cell[j],j,memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'Pendrin':
			solute_id,fluxs=Pendrin.Pendrin(cell[j],memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type =='AE1':
			solute_id,fluxs=AE1.AE1(cell[j],memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'HKATPase':
			solute_id,fluxs=ATPase.hkatpase(cell[j],memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NHE1':
			solute_id,fluxs=NHE1.NHE1(cell[j],memb_id,cell[j].trans[i].act,cell[j].area)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		elif transporter_type == 'NKCC1':
			solute_id,fluxs=NKCC.nkcc1(cell[j],memb_id,cell[j].trans[i].act,delmu)
			for k in range(len(solute_id)):
				file = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[j].segment+'_'+transporter_type+'_'+solute[solute_id[k]]+'.txt','a')
				file.write(str(fluxs[k])+'\n')
		else:
			print('What is this?',transporter_type)
# 	if cell[j].segment == 'mTAL' or cell[j].segment == 'cTAL' or cell[j].segment == 'DCT':
# 		NaH_1_4 = cell[j].area[1,4]*cell[j].dLA[2].perm*(delmu[0,1,4]-delmu[11,1,4])
# 		NaH_1_5 = cell[j].area[1,5]*cell[j].dLA[3].perm*(delmu[0,1,5]-delmu[11,1,5])
# 		file_cell_lis = open(sex_or_preg+'_'+cell[j].segment+'_'+'NaHexchanger_Cell_LIS.txt','a')
# 		file_cell_lis.write(str(NaH_1_4)+'\n')
# 		file_cell_bath = open(sex_or_preg+'_'+cell[j].segment+'_'+'NaHexchanger_Cell_Bath.txt','a')
# 		file_cell_bath.write(str(NaH_1_5)+'\n')

file_Na_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_apical_flux'+sup_or_jux+'.txt','w')
file_K_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_K_apical_flux'+sup_or_jux+'.txt','w')
file_NH4_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_NH4_apical_flux'+sup_or_jux+'.txt','w')
file_H_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H_apical_flux'+sup_or_jux+'.txt','w')
file_Cl_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Cl_apical_flux'+sup_or_jux+'.txt','w')
file_HCO3_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO3_apical_flux'+sup_or_jux+'.txt','w')
file_HPO4_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HPO4_apical_flux'+sup_or_jux+'.txt','w')
file_H2PO4_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H2PO4_apical_flux'+sup_or_jux+'.txt','w')
file_HCO2_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO2_apical_flux'+sup_or_jux+'.txt','w')

file_Na_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_para_flux'+sup_or_jux+'.txt','w')
file_K_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_K_para_flux'+sup_or_jux+'.txt','w')
file_NH4_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_NH4_para_flux'+sup_or_jux+'.txt','w')
file_H_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H_para_flux'+sup_or_jux+'.txt','w')
file_Cl_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Cl_para_flux'+sup_or_jux+'.txt','w')
file_HCO3_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO3_para_flux'+sup_or_jux+'.txt','w')
file_HPO4_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HPO4_para_flux'+sup_or_jux+'.txt','w')
file_H2PO4_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H2PO4_para_flux'+sup_or_jux+'.txt','w')
file_HCO2_flux_para = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO2_para_flux'+sup_or_jux+'.txt','w')

file_Na_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_alpha_flux'+sup_or_jux+'.txt','w')
file_K_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_K_alpha_flux'+sup_or_jux+'.txt','w')
file_NH4_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_NH4_alpha_flux'+sup_or_jux+'.txt','w')
file_H_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H_alpha_flux'+sup_or_jux+'.txt','w')
file_Cl_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Cl_alpha_flux'+sup_or_jux+'.txt','w')
file_HCO3_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO3_alpha_flux'+sup_or_jux+'.txt','w')
file_HPO4_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HPO4_alpha_flux'+sup_or_jux+'.txt','w')
file_H2PO4_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H2PO4_alpha_flux'+sup_or_jux+'.txt','w')
file_HCO2_flux_alpha = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO2_alpha_flux'+sup_or_jux+'.txt','w')
file_Na_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_beta_flux'+sup_or_jux+'.txt','w')
file_K_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_K_beta_flux'+sup_or_jux+'.txt','w')
file_NH4_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_NH4_beta_flux'+sup_or_jux+'.txt','w')
file_H_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H_beta_flux'+sup_or_jux+'.txt','w')
file_Cl_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Cl_beta_flux'+sup_or_jux+'.txt','w')
file_HCO3_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO3_beta_flux'+sup_or_jux+'.txt','w')
file_HPO4_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HPO4_beta_flux'+sup_or_jux+'.txt','w')
file_H2PO4_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_H2PO4_beta_flux'+sup_or_jux+'.txt','w')
file_HCO2_flux_beta = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_HCO2_beta_flux'+sup_or_jux+'.txt','w')

file_K_flux_total = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_K_total_flux'+sup_or_jux+'.txt','w')
file_Na_flux_total = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_total_flux'+sup_or_jux+'.txt','w')
file_Na_flux_baso = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_baso_flux'+sup_or_jux+'.txt','w')
file_Na_flux_Cell_LIS = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_Na_Cell_LIS_flux'+sup_or_jux+'.txt','w')
file_water_flux = open('./'+file_to_save+'/'+sex_or_preg+'_'+species[0:3]+'_'+cell[0].segment+'_water_flux_apical'+sup_or_jux+'.txt','w')

for j in range(1,N):
	jsol=np.zeros([15,6,6])
	jvol,jsol = flux.compute_fluxes(cell[j],j)

	file_water_flux.write('{} {} {} {} \n'.format(jvol[0,1],jvol[0,4],jvol[1,5],jvol[4,5]))

	file_Na_flux.write(str(jsol[0,0,1])+'\n')
	file_K_flux.write(str(jsol[1,0,1])+'\n')
	file_NH4_flux.write(str(jsol[10,0,1])+'\n')
	file_H_flux.write(str(jsol[11,0,1])+'\n')
	file_Cl_flux.write(str(jsol[2,0,1])+'\n')
	file_HCO3_flux.write(str(jsol[3,0,1])+'\n')
	file_HPO4_flux.write(str(jsol[6,0,1])+'\n')
	file_H2PO4_flux.write(str(jsol[7,0,1])+'\n')
	file_HCO2_flux.write(str(jsol[12,0,1])+'\n')
	file_Na_flux_para.write(str(jsol[0,0,4])+'\n')
	file_K_flux_para.write(str(jsol[1,0,4])+'\n')
	file_NH4_flux_para.write(str(jsol[10,0,4])+'\n')
	file_H_flux_para.write(str(jsol[11,0,4])+'\n')
	file_Cl_flux_para.write(str(jsol[2,0,4])+'\n')
	file_HCO3_flux_para.write(str(jsol[3,0,4])+'\n')
	file_HPO4_flux_para.write(str(jsol[6,0,4])+'\n')
	file_H2PO4_flux_para.write(str(jsol[7,0,4])+'\n')
	file_HCO2_flux_para.write(str(jsol[12,0,4])+'\n')
	file_Na_flux_alpha.write(str(jsol[0,0,2])+'\n')
	file_K_flux_alpha.write(str(jsol[1,0,2])+'\n')
	file_NH4_flux_alpha.write(str(jsol[10,0,2])+'\n')
	file_H_flux_alpha.write(str(jsol[11,0,2])+'\n')
	file_Cl_flux_alpha.write(str(jsol[2,0,2])+'\n')
	file_HCO3_flux_alpha.write(str(jsol[3,0,2])+'\n')
	file_HPO4_flux_alpha.write(str(jsol[6,0,2])+'\n')
	file_H2PO4_flux_alpha.write(str(jsol[7,0,2])+'\n')
	file_HCO2_flux_alpha.write(str(jsol[12,0,2])+'\n')
	file_Na_flux_beta.write(str(jsol[0,0,3])+'\n')
	file_K_flux_beta.write(str(jsol[1,0,3])+'\n')
	file_NH4_flux_beta.write(str(jsol[10,0,3])+'\n')
	file_H_flux_beta.write(str(jsol[11,0,3])+'\n')
	file_Cl_flux_beta.write(str(jsol[2,0,3])+'\n')
	file_HCO3_flux_beta.write(str(jsol[3,0,3])+'\n')
	file_HPO4_flux_beta.write(str(jsol[6,0,3])+'\n')
	file_H2PO4_flux_beta.write(str(jsol[7,0,3])+'\n')
	file_HCO2_flux_beta.write(str(jsol[12,0,3])+'\n')
	file_K_flux_total.write(str(jsol[1,0,1]+jsol[1,0,2]+jsol[1,0,3]+jsol[1,0,4])+'\n')
	file_Na_flux_total.write(str(jsol[0,0,1]+jsol[0,0,2]+jsol[0,0,3]+jsol[0,0,4])+'\n')
	file_Na_flux_baso.write(str(jsol[0,4,5])+'\n')
	file_Na_flux_Cell_LIS.write(str(jsol[0,1,4])+'\n')
file_Na_flux_baso.close()
file_Na_flux_Cell_LIS.close()
file_Na_flux.close()
file_K_flux.close()
file_NH4_flux.close()
file_H_flux.close()
file_Cl_flux.close()
file_HCO3_flux.close()
file_HPO4_flux.close()
file_H2PO4_flux.close()
file_HCO2_flux.close()
file_Na_flux_para.close()
file_K_flux_para.close()
file_NH4_flux_para.close()
file_H_flux_para.close()
file_Cl_flux_para.close()
file_HCO3_flux_para.close()
file_HPO4_flux_para.close()
file_H2PO4_flux_para.close()
file_HCO2_flux_para.close()
file_Na_flux_alpha.close()
file_K_flux_alpha.close()
file_NH4_flux_alpha.close()
file_H_flux_alpha.close()
file_Cl_flux_alpha.close()
file_HCO3_flux_alpha.close()
file_HPO4_flux_alpha.close()
file_H2PO4_flux_alpha.close()
file_HCO2_flux_alpha.close()
file_Na_flux_beta.close()
file_K_flux_beta.close()
file_NH4_flux_beta.close()
file_H_flux_beta.close()
file_Cl_flux_beta.close()
file_HCO3_flux_beta.close()
file_HPO4_flux_beta.close()
file_H2PO4_flux_beta.close()
file_HCO2_flux_beta.close()
file_K_flux_total.close()
file_Na_flux_total.close()
file_water_flux.close()
