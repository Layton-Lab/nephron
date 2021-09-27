# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 10:41:12 2021

@author: melis
"""

import os
import matplotlib.pyplot as plt
import numpy as np

#===================
# settings
#===================
compare = 1# 1,2,or 3 comparisons

segs_early = ['PT', 'S3', 'SDL']
segs_jux = ['LDL', 'LAL']
segs_later = ['mTAL', 'cTAL', 'DCT', 'CNT']

sup_list = segs_early + segs_later
jux_list = segs_early + segs_jux + segs_later

segs_cd = ['CCD', 'OMCD', 'IMCD']

direct1 = 'male-nephron'
sex1 = 'male'

direct2 = 'female-original'
sex2 = 'female'

direct3 =''
sex3 = ''




label1 = direct1
label2 = direct2
label3 = direct3

solute_list = ['Na', 'water', 'K', 'Cl']
#solute_list = ['Na', 'K', 'Cl', 'glu', 'HCO3', 'NH4','TA', 'water'] 


# shift for the horizontal labeling lines
shift = -1
#======================
#======================

def get_data(direct, sex, solute, segment, sup_or_jux):
    os.chdir(direct)
    y = []
    if solute == 'TA':
        y = get_data_TA(direct, sex, segment, sup_or_jux)
    else:
        if solute == 'water':
            filename = sex +'_rat_'+segment+'_water_volume_in_Lumen'+sup_or_jux+'.txt'
        elif solute == 'pH':
            filename = sex+'_rat_'+segment+'_pH_in_Lumen.txt'
        else:
            filename = sex + '_rat_'+segment+'_flow_of_'+solute+'_in_Lumen'+sup_or_jux+'.txt'
        file = open(filename, 'r')
        y.extend(np.loadtxt(filename, delimiter = '\n', unpack = True))
        file.close()
    os.chdir('..')
    return y

# [TA] = (10**(7.4-6.8)*[H2PO4]-[HPO4])/(1+10**(7.4-6.8))
def get_data_TA(direct, sex, segment, sup_or_jux):
    os.chdir(direct)
    y = []
    fname1 = sex + '_rat_'+segment+'_flow_of_H2PO4_in_Lumen'+sup_or_jux+'.txt'
    fname2 = sex + '_rat_'+segment+'_flow_of_HPO4_in_Lumen'+sup_or_jux+'.txt'
    file1 = open(fname1, 'r')
    file2 = open(fname2, 'r')
    H2PO4_vals = np.loadtxt(fname1, delimiter = '\n', unpack = True)
    HPO4_vals = np.loadtxt(fname2, delimiter = '\n', unpack = True)
    file1.close()
    file2.close()
    TA_vals = (10**(7.4-6.8) * H2PO4_vals - HPO4_vals)/(1 + 10**(7.4-6.8))
    y.extend(TA_vals)
    os.chdir('..')
    return y

def list_data(seg_list, solute, sup_or_jux):
    data1 = []
    data2 = []
    data3 = []
    for seg in range(len(seg_list)):
        segment = seg_list[seg]
        print(segment)
        data1.extend(get_data(direct1, sex1, solute, segment, sup_or_jux))
        if compare > 1:
            data2.extend(get_data(direct2, sex2, solute, segment, sup_or_jux))
            if compare > 2:
                data3.extend(get_data(direct3, sex3, solute, segment, sup_or_jux))
    return data1, data2, data3
#=============
# get data
#=============
for sol in range(len(solute_list)):
    solute = solute_list[sol]
    print('\n')
    print(solute)
    
    print('sup data')
    sup_data1, sup_data2, sup_data3 = list_data(sup_list, solute, '_sup')
    
    print('jux data')
    jux_data1, jux_data2, jux_data3 = list_data(jux_list, solute, '_jux5')
    
    print('cd data')
    cd_data1, cd_data2, cd_data3 = list_data(segs_cd, solute, '')

#====================================
# position/length of each segment
#====================================
    sex = sex1
    if sex == 'male':
        lens_sup = {'PT': 1.1, 'SDL': 0.14, 'mTAL': 0.2, 'cTAL': 0.2,
                    'DCT': 0.1, 'CNT': 0.2}
        # jux lengths that are different from sup nephron
        lens_jux = {'LDL': 0.5, 'LAL': 0.5, 'cTAL': 0.05, 'CNT': 0.3}
        lens_cd = {'CCD':0.2, 'OMCD':0.2, 'IMCD':0.425}
    else:
        print('sex: '+ sex)
        raise Exception('sex segment lengths not set up yet')
        
    pos_PT = [lens_sup['PT']*i/199 for i in range(200)] #this is PCT/S3
    pos_SDL = [pos_PT[-1]+lens_sup['SDL']*i/199 for i in range(200)]
    pos_LDL = [pos_SDL[-1]+lens_jux['LDL']*i/199 for i in range(200)]
    pos_LAL = [pos_LDL[-1]+lens_jux['LAL']*i/199 for i in range(200)]
    pos_mTAL = [pos_LAL[-1]+lens_sup['mTAL']*i/199 for i in range(200)]
    # cTAL length for sup and jux different, sup longer
    pos_cTAL_sup = [pos_mTAL[-1]+lens_sup['cTAL']*i/199 for i in range(200)]
    pos_cTAL_jux = [pos_mTAL[-1]+lens_jux['cTAL']*i/199 for i in range(200)]
    pos_DCT = [pos_cTAL_sup[-1]+lens_sup['DCT']*i/199 for i in range(200)]
    # CNT length for sup and jux different, jux longer
    pos_CNT_sup = [pos_DCT[-1]+lens_sup['CNT']*i/199 for i in range(200)]
    pos_CNT_jux = [pos_DCT[-1]+lens_jux['CNT']*i/199 for i in range(200)]
    pos_CCD = [pos_CNT_jux[-1]+lens_cd['CCD']*i/199 for i in range(200)]
    pos_OMCD = [pos_CCD[-1]+lens_cd['OMCD']*i/199 for i in range(200)]
    pos_IMCD = [pos_OMCD[-1]+lens_cd['IMCD']*i/199 for i in range(200)]
    
    
    sup_early_pos = pos_PT + pos_SDL
    sup_late_pos = pos_mTAL + pos_cTAL_sup + pos_DCT + pos_CNT_sup
    jux_pos_early = pos_PT + pos_SDL + pos_LDL + pos_LAL + pos_mTAL + pos_cTAL_jux
    jux_pos_late = pos_DCT + pos_CNT_jux
    cd_pos = pos_CCD + pos_OMCD + pos_IMCD
    
    #==========================
    # make plots
    #=========================
    # note that the plots are for jux5 and sup nephron
    
    
    fig, ax = plt.subplots()
    plt.rcParams['axes.xmargin'] =0
    # fig settings
    #fig.set_figheight(12)
    #fig.set_figwidth(12)
       
    # colors
    c1 = 'c'
    c2 = 'mediumvioletred'
    c3 = 'green'
    
    # fontsizes
    xlab_size = 14
    xticklab_size = 14
    ylab_size = 14
    yticklab_size = 14
    title_size = 16
    leg_size = 14
    text_size = 12
    
    
    sup_early1 = ax.plot(sup_early_pos, sup_data1[0:400], color = c1, label = label1)
    sup_late1 = ax.plot(sup_late_pos, sup_data1[401:1400], color = c1)
    jux_early1= ax.plot(jux_pos_early, jux_data1[0:1200], color = c1, linestyle = 'dashed')
    jux_late1 = ax.plot(jux_pos_late, jux_data1[1200:1600], color = c1, linestyle = 'dashed')
    cd1 = ax.plot(cd_pos, cd_data1, color = c1)
    
    if compare > 1:
        sup_early2 = ax.plot(sup_early_pos, sup_data2[0:400], color = c2, label = label2)
        sup_late2 = ax.plot(sup_late_pos, sup_data2[401:1400], color = c2)
        jux_early2= ax.plot(jux_pos_early, jux_data2[0:1200], color = c2, linestyle = 'dashed')
        jux_late2 = ax.plot(jux_pos_late, jux_data2[1200:1600], color = c2, linestyle = 'dashed')
        cd2 = ax.plot(cd_pos, cd_data2, color = c2)
    
    print('dashed is jux5 values')
    
    plt.yticks(fontsize=yticklab_size)
    plt.xticks(fontsize=xticklab_size)
    plt.title(solute + ' flow', fontsize = title_size)
    plt.legend(fontsize = leg_size)
    ax = plt.gca()
    if solute == 'water':
        ax.set_ylabel('volume flow (nl/min)', fontsize = ylab_size)
    else:
        ax.set_ylabel(solute + ' flow (pmol/min)', fontsize = ylab_size)
    ax.set_xlabel('nephron length (normalized)', fontsize = xlab_size)
    ax.axes.xaxis.set_ticks([])
    
    height = sup_data1[0]+shift
    plt.hlines(height, pos_PT[0], pos_PT[-1],colors = 'k', linestyles = 'dotted')
    plt.text(0.4, height+0.1, 'PT', fontsize = text_size)
    
    height = sup_data1[200]-1
    plt.hlines(height, pos_SDL[0], pos_SDL[-1], colors = 'k', linestyles = 'dotted')
    plt.text(1.0, height+0.15, 'SDL', fontsize = text_size)
    
    height = jux_data1[400]+shift+1.5
    plt.hlines(height, pos_LDL[0], pos_LDL[-1], colors = 'k', linestyles = 'dotted')
    plt.text(1.35, height+0.15, 'LDL', fontsize = text_size)
    
    height = jux_data1[600]+shift
    plt.hlines(height, pos_LAL[0], pos_LAL[-1], colors = 'k', linestyles = 'dotted')
    plt.text(1.87, height+0.15, 'LAL')
    
    height = sup_data1[400]+shift
    plt.hlines(height, pos_mTAL[0], pos_cTAL_sup[-1], colors = 'k', linestyles = 'dotted')
    plt.text(2.3, height+0.15, 'TAL')
    
    height = sup_data1[600]+shift-1
    plt.hlines(height, pos_DCT[0], pos_DCT[-1], colors = 'k', linestyles = 'dotted')
    plt.text(2.6, height-0.5, 'DCT')
    
    height = sup_data1[800]+shift+1
    plt.hlines(height, pos_CNT_jux[0], pos_CNT_jux[-1], colors = 'k', linestyles = 'dotted')
    plt.text(2.8, height+0.15, 'CNT')
    
    height = cd_data1[200]+shift
    plt.hlines(height, cd_pos[0], cd_pos[-1], colors = 'k', linestyles = 'dotted')
    plt.text(3.3, height+0.15, 'CD')
