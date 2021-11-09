# this is to set the parameters for the pregnancy model for the rat
# note that diabetes and unx has not been characterized so then there are no diabetes/unx options in this file
import re
from defs import *
from values import *
import time

# pick out compartment IDs. E.g. given label=Area_Lumen_Cell, it would return
# (0,1) for (lumen, cell)
def get_interface_id(label):
    tmp = (label).split('_')
    if len(tmp)==3:  # only have two compartment IDs
        ind1,ind2 = compart_id[tmp[1]],compart_id[tmp[2]]
        return ind1, ind2
    elif len(tmp)==4:  # solute ID, followed by two compartment IDs
        sid,ind1,ind2 = solute_id[tmp[1]],compart_id[tmp[2]],compart_id[tmp[3]]
        return sid,ind1,ind2

# for coupled transpoters
# pick out solute IDs, compartment IDs, and coefficients
def get_coupled_id(label):
    tmp = (label).split('_')
#    print(tmp)
    sid1,sid2 = solute_id[tmp[3]],solute_id[tmp[4]]
    ind1,ind2 = compart_id[tmp[1]],compart_id[tmp[2]]
    if len(tmp)==5:  # only involves 2 solutes
        return ind1,ind2,sid1,sid2
    elif len(tmp)==6:  # involves 3 solutes
        sid3 = solute_id[tmp[5]]
        return ind1,ind2,sid1,sid2,sid3
    else:
        print("Wrong label",tmp)

# is str_short at the beginning of str_long?
# case insensitive
def compare_string_prefix(str_long,str_short):
    return (str_long.lower())[:len(str_short)] == str_short.lower()

def compare_sex(sex, cell):
    if sex.lower() == 'male' or sex.lower() == 'female':
        return cell.sex.lower() == sex.lower()
    else:
        return True

def read_params_preg(cell,filename,j):

    file = open(filename,'r')
    
    cell.segment=filename[12:-16]

    # error messages
    if cell.preg == 'non':
        print('cell.preg: ' + cell.preg)
        raise Exception('read_params_preg is for pregnant model')
    if cell.diabete != 'Non':
        print('cell.diabete: ' + cell.diabete)
        raise Exception('diabetes during pregnancy has not been done yet')
    if cell.unx == 'Y':
        print('cell.unx:' + cell.unx)
        raise Exception('unx has not been done for pregnancy yet')
    if cell.species != 'rat':
        print('cell.species: ' + cell.species)
        raise Exception('pregnancy only done for rat, human model not set up yet')
    
    line = file.readline()
    while (line):
        line = line.replace('\t',' ')
        terms = line.split(' ')
        if line[0][0]!='#':
            id = terms[0] #re.findall(r'[A-Za-z_]+', line)
            sex = id.split('_')[-1].lower()
            if compare_sex(sex, cell):
                id = id.replace('_male', '')
                id = id.replace('_female', '')
            else:
                line = file.readline()
                continue;
                
            # Skip over the label, which may contain numbers like in HCO3
            first_space_pos = line.index(' ')
            num = re.findall(r'-?\d+\.?\d*[Ee]?[+-]?\d*', line[first_space_pos:len(line)])
            if num: # if this line is numerical parameter
                value = float(num[0])


            if id.lower() == "Sex".lower():
                find = re.findall("Female".lower(), terms[-1].lower())
                if find:
                    cell.sex = "Female".lower()
            
            # Diameter:
            elif compare_string_prefix(id,"Diameter"):
                # pregnant diameter
                preg_rat = 1.0 # reset preg_rat
                if cell.segment == 'PT' or cell.segment == 'S3':
                    if cell.preg == 'mid':
                        preg_rat = 1.07
                    elif cell.preg == 'late':
                        preg_rat = 1.07
                else:
                    if cell.preg == 'mid':
                        preg_rat = 1.05
                    elif cell.preg == 'late':
                        preg_rat = 1.05

                if cell.type != 'sup':
                    if cell.preg == 'mid':
                        preg_rat = preg_rat*1.03
                    elif cell.preg == 'late':
                        preg_rat = preg_rat*1.03

                cell.diam = value*preg_rat

            # Length:
            elif compare_string_prefix(id,"Length"):
                # pregnant PT length
                if cell.segment == 'PT' or cell.segment == 'S3':
                    if cell.preg == 'mid':
                        cell.len = value*1.1425
                    elif cell.preg == 'late':
                        cell.len = value*1.165
                # juxtamedullary segments lengths
                elif cell.segment == 'LDL' or cell.segment == 'LAL':
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

                    cell.len = value*looplen
                else:
                    cell.len = value

                if cell.type != 'sup' and cell.species == 'rat':
                    if cell.segment == 'cTAL':
                        if cell.sex == 'male':
                            cell.len = 0.05
                        elif cell.sex == 'female':
                            cell.len = 0.05*0.9 #updated female
                    elif cell.segment == 'CNT':
                        if cell.sex == 'male':
                            cell.len = 0.3
                        elif cell.sex == 'female':
                            cell.len = 0.3*0.9 #updated female

                if cell.type != 'sup' and cell.species == 'hum':
                    if cell.segment == 'cTAL':
                        if cell.sex == 'male':
                            cell.len = 0.125
                        elif cell.sex == 'female':
                            cell.len = 0.125
                    elif cell.segment == 'CNT':
                        if cell.sex == 'male':
                            cell.len = 0.6
                        elif cell.sex == 'female':
                            cell.len = 0.6

            # Total number of cells:
            elif compare_string_prefix(id,"Total"):
                cell.total = value
                                      
            # Luminal pressure:
            elif compare_string_prefix(id,"Pressure"):
                cell.pres[0] = value
                if cell.type !='sup' and cell.segment == 'PT' and cell.species == 'rat':
                    if cell.preg == 'mid':
                        cell.pres[0] = 12.75
                    elif cell.preg == 'late':
                        cell.pres[0] = 12.75

            # pH:
            elif compare_string_prefix(id,"pH"):
                for i in range(6):
                    cell.pH[i] = num[i]
           
            # Surface area multiplication factor:
            elif compare_string_prefix(id,"Area"):
                ind1,ind2 = get_interface_id(id)
                cell.area[ind1][ind2] = value
                cell.area[ind2][ind1] = value  # symmetry
                cell.area_init[ind1][ind2] = value
                cell.area_init[ind2][ind1] = value
                if cell.type != 'sup' and cell.species == 'rat':
                    if cell.segment == 'PT' or cell.segment == 'S3':
                        cell.area[ind1][ind2] = 1.75*cell.area[ind1][ind2]
                        cell.area[ind2][ind1] = cell.area[ind1][ind2]

            # Water permeabilities:
            elif compare_string_prefix(id,"Pf"):
                ind1,ind2 = get_interface_id(id)
                # Units of dimensional water flux (in 'value'): cm3/s/cm2 epith
                # Non-dimensional factor for water flux: (Pfref)*Vwbar*Cref
                # Calculate non-dimensional dLPV = Pf*Vwbar*Cref / (Pfref*Vwbar*Cref)
                # dLPV = Pf/Pfref

                preg_rat = 1.0 #reset to 1.0, will change if needed

                if cell.preg != 'non':
                    # pregnancy water perm (transcellular)
                    if ind1 == 0 and ind2 == 1:
                        if cell.segment == 'PT' or cell.segment == 'S3':
                            if cell.preg == 'mid':
                                preg_rat = 1.0
                            elif cell.preg == 'late':
                                preg_rat = 1.0
                        elif cell.segment == 'SDL':
                            if cell.preg == 'mid':
                                preg_rat = 1.1
                            elif cell.preg == 'late':
                                preg_rat = 1.5
                        elif cell.segment == 'LDL':
                            if cell.preg == 'mid':
                                preg_rat = 1.1
                            elif cell.preg == 'late':
                                preg_rat = 1.5
                        elif cell.segment == 'CCD':
                            if cell.preg == 'mid':
                                preg_rat = 1.4
                            elif cell.preg == 'late':
                                preg_rat = 1.4
                        elif cell.segment == 'OMCD':
                            if cell.preg == 'mid':
                                preg_rat = 1.4
                            elif cell.preg == 'late':
                                preg_rat = 1.4
                        elif cell.segment == 'IMCD':
                            if cell.preg == 'mid':
                                preg_rat = 1.8
                            elif cell.preg == 'late':
                                preg_rat = 1.8
                    elif ind1 == 1:
                        if ind2 == 4 or ind2 == 5:
                            if cell.segment == 'PT' or cell.segment == 'S3':
                                if cell.preg == 'mid':
                                    preg_rat = 1.0
                                elif cell.preg == 'late':
                                    preg_rat = 1.0
                            elif cell.segment == 'SDL':
                                if cell.preg == 'mid':
                                    preg_rat = 1.1
                                elif cell.preg == 'late':
                                    preg_rat = 1.5
                            elif cell.segment == 'LDL':
                                if cell.preg == 'mid':
                                    preg_rat = 1.1
                                elif cell.preg == 'late':
                                    preg_rat = 1.5
                            elif cell.segment == 'CCD':
                                if cell.preg == 'mid':
                                    preg_rat = 1.4
                                elif cell.preg == 'late':
                                    preg_rat = 1.4
                            elif cell.segment == 'OMCD':
                                if cell.preg == 'mid':
                                    preg_rat = 1.4
                                elif cell.preg == 'late':
                                    preg_rat = 1.4
                            elif cell.segment == 'IMCD':
                                if cell.preg == 'mid':
                                    preg_rat = 1.8
                                elif cell.preg == 'late':
                                    preg_rat = 1.8

                cell.dLPV[ind1][ind2] = value/Pfref*preg_rat
                #print('water permeability')
                #print(ind1, ind2, cell.dLPV[ind1][ind2],preg_rat)
                # symmetry
                cell.dLPV[ind2][ind1] = value/Pfref*preg_rat

                if cell.segment == 'SDL' and cell.type == 'sup':
                    if j>=0.46*cell.total:
                        cell.dLPV[0,1]=0.00*cell.dLPV[0,1]
                        cell.dLPV[0,4]=0.00*cell.dLPV[0,4]
                elif cell.segment == 'LDL':
                    if cell.sex == 'male':
                        if j>=0.4*cell.total:
                            cell.dLPV[0,1]=0.00*cell.dLPV[0,1]
                            cell.dLPV[0,4]=0.00*cell.dLPV[0,4]
                    elif cell.sex == 'female':
                        if j>=0.5*cell.total:
                            cell.dLPV[0,1]=0.00*cell.dLPV[0,1]
                            cell.dLPV[0,4]=0.00*cell.dLPV[0,4]
                
                if cell.segment == 'CNT' and cell.type !='sup':
                    if cell.sex == 'female':
                        cell.dLPV = cell.dLPV*4/3
                                
            # Reflection coefficients:
            elif compare_string_prefix(id,"sig"):
                sid,ind1,ind2 = get_interface_id(id)
                cell.sig[sid][ind1][ind2] = value
                # symmetry
                cell.sig[sid][ind2][ind1] = value
                cell.sig[sid][4][5] = 0.0
                cell.sig[sid][5][4] = 0.0
                
            # Membrane solute permeabilities:
            elif compare_string_prefix(id,"perm"):
                sid,ind1,ind2 = get_interface_id(id)
                cell.h[sid][ind1][ind2] = value*1.0e-5/href
                # Symmetry:
                cell.h[sid][ind2][ind1] = value*1.0e-5/href
                # Same permeability on basolateral membrane (around bath or LIS):
                if ind1==1 and ind2==5:
                    cell.h[sid][ind1][4] = value*1.0e-5/href
                    cell.h[sid][4][ind1] = value*1.0e-5/href
                elif ind1==1 and ind2==4:
                    cell.h[sid][ind1][5] = value*1.0e-5/href
                    cell.h[sid][5][ind1] = value*1.0e-5/href
                elif ind1==2 and ind2==4:
                    cell.h[sid][ind1][5] = value*1.0e-5/href
                elif ind1==3 and ind2==4:
                    cell.h[sid][ind1][5] = value*1.0e-5/href
                if cell.segment == 'OMCD':
                    cell.h[sid][0][3] = 0.0
                    cell.h[sid][3][4] = 0.0
                    cell.h[sid][3][5] = 0.0
                if cell.segment == 'IMCD':
                    if cell.sex == 'male':
                        if j>3*cell.total/4-1:
                            cell.h[8,0,1] = 300.0*1.0e-5/href
                    elif cell.sex == 'female':
                        if j>2*cell.total/3-1:
                            cell.h[8,0,1] = 300.0*1.0e-5/href
                if cell.segment == 'LDL':
                    if cell.sex == 'male':
                        if j>=0.4*cell.total:
                            cell.h[0,0,1]=80.0
                            cell.h[0,0,4]=80.0
                            cell.h[1,0,1]=100.0
                            cell.h[1,0,4]=100.0
                            cell.h[2,0,1]=80.0
                            cell.h[2,0,4]=80.0
                            cell.h[3,0,1]=20.0
                            cell.h[3,0,4]=20.0
                            cell.h[10,0,1]=20.0
                            cell.h[10,0,4]=20.0
                            cell.h[8,0,1]=80.0 #80
                            cell.h[8,0,4]=80.0 #80
                    elif cell.sex == 'female':
                        if j>=0.5*cell.total:
                            cell.h[0,0,1]=40.0
                            cell.h[0,0,4]=80.0
                            cell.h[1,0,1]=100.0
                            cell.h[1,0,4]=100.0
                            cell.h[2,0,1]=40.0
                            cell.h[2,0,4]=80.0
                            cell.h[3,0,1]=20.0
                            cell.h[3,0,4]=20.0
                            cell.h[10,0,1]=20.0
                            cell.h[10,0,4]=20.0
                            cell.h[8,0,1]=80.0
                            cell.h[8,0,4]=80.0

                # K secretion
                if cell.segment == 'DCT':
                    if j>0.66*cell.total:
                        #DCT2
                        if cell.preg == 'late':
                            preg_rat = 0.45
                            cell.h[1,0,1] = 0.6*preg_rat
                        elif cell.preg == 'mid':
                            preg_rat = 0.5
                            cell.h[1,0,1] = 0.6*preg_rat
                elif cell.segment == 'CNT':
                    if cell.preg == 'late':
                        preg_rat = 0.45
                        cell.h[1,0,1] = 8.0*preg_rat
                    elif cell.preg == 'mid':
                        preg_rat = 0.5
                        cell.h[1,0,1] = 8.0*preg_rat
                elif cell.segment == 'CCD':
                    if cell.preg == 'late':
                        preg_rat = 0.7
                        cell.h[1,0,1] = 2.8*preg_rat
                    elif cell.preg == 'mid':
                        preg_rat = 0.75
                        cell.h[1,0,1] = 2.8*preg_rat
                elif cell.segment == 'OMCD':
                    if cell.preg == 'late':
                        preg_rat = 0.7
                        cell.h[1,0,1] = 2.4*preg_rat
                    elif cell.preg == 'mid':
                        preg_rat = 0.75 
                        cell.h[1,0,1] = 2.4*preg_rat
                    
                            
            # Coupled transporters:
            elif compare_string_prefix(id,"coupled"):
                # retrieve interface and solute id
                vals = get_coupled_id(id)
                newdLA = coupled_transport()
                newdLA.perm = value / (href*Cref)
                newdLA.membrane_id = [vals[0],vals[1]]
                
                coef = []  # retrieve coupling coefficients

                for i in range(1,len(num)):
                    coef.append(int(num[i]))
                newdLA.coef = coef
                newdLA.solute_id = vals[2:len(vals)]
                
                # NaPi2
                if newdLA.solute_id == (0,7):
                    if cell.segment == 'PT' or cell.segment == 'S3':
                        if cell.preg == 'mid':
                            newdLA.perm = 0.9*newdLA.perm 
                        elif cell.preg == 'late':
                            newdLA.perm = 0.85*newdLA.perm
                    else:
                        print('segment: ' + cell.segment)
                        raise Exception('NaPi2 in pregnancy not characterized for this segment')
                # K-Cl cotransporter
                elif newdLA.solute_id == (1,2):
                    if cell.segment == 'PT' or cell.segment == 'S3':
                        if cell.preg == 'mid':
                            newdLA.perm = 1.35*newdLA.perm
                        elif cell.preg == 'late':
                            newdLA.perm = 1.3*newdLA.perm
                    elif cell.segment == 'DCT':
                        if cell.preg == 'mid':
                            newdLA.perm = 1.4*newdLA.perm
                        elif cell.preg == 'late':
                            newdLA.perm = 1.3*newdLA.perm
                    elif cell.segment == 'IMCD':
                        if cell.preg == 'mid':
                            newdLA.perm = 1.4*newdLA.perm
                        elif cell.preg == 'late':
                            newdLA.perm = 1.3*newdLA.perm
                    else:
                        print('segment: '+cell.segment)
                        raise Exception('K-Cl coupled transporter not characterized for pregnancy in this segment')
                #Na-Cl cotransporter
                elif newdLA.solute_id == (0,2):
                    if cell.segment == 'OMCD':
                        if cell.preg == 'mid':
                            newdLA.perm = 1.15*newdLA.perm
                        elif cell.preg == 'late':
                            newdLA.perm = 1.05*newdLA.perm
                    elif cell.segment == 'IMCD':
                        if cell.preg == 'mid':
                            newdLA.perm = 1.15*newdLA.perm
                        elif cell.preg == 'late':
                            newdLA.perm = 1.05*newdLA.perm
                    else:
                        print('segment: ' + cell.segment)
                        raise Exception('Na-Cl coupled transporter not characterized for pregnancy in this segment')

                cell.dLA.append(newdLA)


            # Specific transporters:
            elif compare_string_prefix(id,"transport"):
                tmp = (id).split('_')
                ind1,ind2 = compart_id[tmp[1]],compart_id[tmp[2]]
                newTransp = transporter()
                newTransp.membrane_id = [ind1,ind2]
                newTransp.type = tmp[3]
                newTransp.act = value/(href*Cref)
                #print('transporter')
                #print(newTransp.membrane_id,newTransp.type,newTransp.act)
                if cell.type != 'sup' and cell.sex == 'female' and cell.species == 'rat':
                    if cell.segment == 'mTAL' or cell.segment == 'cTAL':
                        if newTransp.type == 'NKCC2A' or newTransp.type == 'NKCC2B' or newTransp.type == 'NKCC2F' or newTransp.type == 'NaKATPase':
                            newTransp.act = 1.5*value/(href*Cref)
                        elif newTransp.type == 'KCC4':
                            newTransp.act = 2.0*value/(href*Cref)

                # pregnant model values
                if newTransp.type == 'NHE3':
                    # PCT, S3, mTAL, cTAL, DCT
                    if cell.segment == 'PT' or cell.segment == 'S3':
                        if cell.preg == 'mid':
                            preg_rat = 1.3
                        elif cell.preg == 'late':
                            preg_rat = 1.175
                    elif cell.segment == 'mTAL' or cell.segment == 'cTAL' or cell.segment == 'DCT':
                        if cell.preg == 'mid':
                            preg_rat = 1.3
                        elif cell.preg == 'late':
                            preg_rat = 1.175
                    else:
                        print('segment: ' + cell.segment)
                        raise Exception('NHE3 activity not done for pregnancy in this segment')
                elif newTransp.type == 'NaKATPase':
                    if cell.segment == 'PT' or cell.segment == 'S3' or cell.segment == 'cTAL':
                        if cell.preg == 'mid':
                            preg_rat = 0.725
                        elif cell.preg == 'late':
                            preg_rat = 0.65
                    elif cell.segment == 'DCT':
                        if cell.preg == 'mid':
                            preg_rat = 0.75
                        elif cell.preg == 'late':
                            preg_rat = 0.75
                    elif cell.segment == 'CNT':
                        if cell.preg == 'mid':
                            preg_rat = 0.75
                        elif cell.preg == 'late':
                            preg_rat = 0.75
                    elif cell.segment == 'mTAL':
                        if cell.preg == 'mid':
                            preg_rat = 1.15
                        elif cell.preg == 'late':
                            preg_rat = 1.0
                    elif cell.segment == 'CCD':
                        if cell.preg == 'mid':
                            preg_rat = 0.75
                        elif cell.preg == 'late':
                            preg_rat = 0.7
                    elif cell.segment == 'IMCD' or cell.segment == 'OMCD':
                        if cell.preg == 'mid':
                            preg_rat = 1.125
                        elif cell.preg == 'late':
                            preg_rat = 1.0
                    else:
                        print('segment: ' + cell.segment)
                        raise Exception('NaKATPase activity not done for pregnancy in segment')
                elif newTransp.type == 'NKCC2A' or newTransp.type == 'NKCC2B' or newTransp.type == 'NKCC2F':
                    if cell.preg == 'mid':
                        preg_rat = 1.15
                    elif cell.preg == 'late':
                        preg_rat = 1.5
                elif newTransp.type == 'KCC4':
                    if cell.preg == 'mid':
                        preg_rat = 1.4
                    elif cell.preg == 'late':
                        preg_rat = 1.35 # 1.35
                elif newTransp.type == 'NCC':
                    if cell.preg == 'mid':
                        preg_rat = 1.0
                    elif cell.preg == 'late':
                        preg_rat = 0.9
                elif newTransp.type == 'ENaC':
                    if cell.preg == 'mid':
                        preg_rat = 1.85
                    elif cell.preg == 'late':
                        preg_rat = 2.15
                elif newTransp.type == 'HKATPase':
                    if cell.preg == 'mid':
                        preg_rat = 2.5
                    elif cell.preg == 'late':
                        preg_rat = 2.75
                elif newTransp.type == 'HATPase':
                    if cell.preg == 'mid':
                        preg_rat = 1.0
                    elif cell.preg == 'late':
                        preg_rat = 1.0
                elif newTransp.type == 'Pendrin':
                    if cell.preg == 'mid':
                        preg_rat = 1.65
                    elif cell.preg == 'late':
                        preg_rat = 1.7
                elif newTransp.type == 'NHE1':
                    if cell.preg == 'mid':
                        preg_rat = 0.9
                    elif cell.preg == 'late':
                        preg_rat = 0.9
                elif newTransp.type == 'AE1':
                    if cell.preg == 'mid':
                        preg_rat = 1.0 
                    elif cell.preg == 'late':
                        preg_rat = 1.0
                else:
                    preg_rat = 1.0
                newTransp.act = preg_rat*newTransp.act
                #print(cell.preg + 'pregnant transporter activity: ' + str(newTransp.act))
                cell.trans.append(newTransp)

            # Solute concentrations:
            elif compare_string_prefix(id,"conc"):
                tmp = (id).split('_')
                cell.conc[solute_id[tmp[1]]][0] = float(num[0])
                cell.conc[solute_id[tmp[1]]][1] = float(num[1])
                cell.conc[solute_id[tmp[1]]][4] = float(num[2])
                cell.conc[solute_id[tmp[1]]][5] = float(num[3])
                if len(num) > 4:
                    cell.conc[solute_id[tmp[1]]][2] = float(num[4])
                    if len(num) > 5:
                        cell.conc[solute_id[tmp[1]]][3] = float(num[5])

            # Reference impermeat concentration (for cell)
            # or oncotic pressure for lumen/LIS/bath
            elif compare_string_prefix(id,"cimpref"):
                tmp = (id).split('_')
                cell.cimpref[compart_id[tmp[1]]] = float(num[0])

            # Impermeant properties
            elif compare_string_prefix(id,"zimp"):
                tmp = (id).split('_')
                cell.zimp[compart_id[tmp[1]]] = float(num[0])
                
            # Reference buffer concentrations (for cell):
            elif compare_string_prefix(id,"cbuftot"):
                tmp = (id).split('_')
                cell.cbuftot[compart_id[tmp[1]]] = float(num[0])
                
            # Rates used in HCO3/H2CO3 reaction:
            elif compare_string_prefix(id,"dkd"):
                tmp = (id).split('_')
                cell.dkd[compart_id[tmp[1]]] = float(num[0])
            elif compare_string_prefix(id,"dkh"):
                tmp = (id).split('_')
                cell.dkh[compart_id[tmp[1]]] = float(num[0])
            
            # Reference volume flows:
            elif compare_string_prefix(id,"volref"):
                tmp = (id).split('_')
                cell.volref[compart_id[tmp[1]]] = float(num[0])

            # Actual volume flows:
            elif compare_string_prefix(id,"vol"):
                tmp = (id).split('_')  
                if cell.segment == 'PT' and cell.type == 'sup':
                    # SNGFR for sup nephrons
                    if compart_id[tmp[1]] == 0:
                        if cell.preg == 'mid':
                            cell.vol[0] = 0.0052 #0.004*1.3
                        elif cell.preg == 'late':
                            cell.vol[0] = 0.0048 #0.004*1.2
                        cell.vol_init[0] = cell.vol[0]
                    else:
                        cell.vol[compart_id[tmp[1]]] = float(num[0]) 
                        cell.vol_init[compart_id[tmp[1]]] = float(num[0])
                elif cell.segment == 'PT' and cell.type != 'sup' and cell.species == 'rat':
                    # SNGFR for jux nephrons
                    if compart_id[tmp[1]] == 0:
                        if cell.preg == 'mid':
                            cell.vol[0] = 0.00728 #0.0056*1.3
                        elif cell.preg == 'late':
                            cell.vol[0] = 0.00672 #0.0056*1.2
                        cell.vol_init[0] = cell.vol[0]
                    else:
                        cell.vol[compart_id[tmp[1]]] = float(num[0])
                        cell.vol_init[compart_id[tmp[1]]] = float(num[0])
                else:
                    cell.vol[compart_id[tmp[1]]] = float(num[0])
                    cell.vol_init[compart_id[tmp[1]]] = float(num[0])

            # Membrane potential:
            elif compare_string_prefix(id,"ep"):
                tmp = (id).split('_')
                cell.ep[compart_id[tmp[1]]] = float(num[0])
            
            # Interstitial concentration parameters
            # CM: cortico-medullary boundary; OI: outer-inner stripe boundary; Pap: papillary tip
            # changes for altered concentrations (i.e., in pregnancy) are made in boundaryBath
            elif compare_string_prefix(id,'cm'):
                cell.cm = []
                for i in num:
                    cell.cm.append(float(i))
            elif compare_string_prefix(id,'oi'):
                cell.oi = []
                for i in num:
                    cell.oi.append(float(i))
            elif compare_string_prefix(id,'pap'):
                cell.pap = []
                for i in num:
                    cell.pap.append(float(i))
                        
            # invalue keyword
#            else:
#                print("Wrong id",id)
                # if not a comment line
        
        line = file.readline()

        # Updating the reflective coefficients at the LIS-Bath interface:
        for i in range(NS):
            cell.sig[i,4,5] = 0.0
            cell.sig[i,5,4] = 0.0  
    file.close()
