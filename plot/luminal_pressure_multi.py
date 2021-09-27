# -*- coding: utf-8 -*-
"""
This is for plotting the luminal pressure of the superficial
and each of the juxtamedullary nephrons for a given model
output
"""

import numpy as np
import matplotlib.pyplot as plt
import os

#===========================
# settings
#===========================
direct = 'male-torqR'
sex = 'male'

segs_early = ['PT', 'S3', 'SDL']
segs_jux = ['LDL', 'LAL']
segs_later = ['mTAL', 'cTAL', 'DCT', 'CNT']

sup_list = segs_early + segs_later
jux_list = segs_early + segs_jux + segs_later

segs_cd = ['CCD', 'OMCD', 'IMCD']

shift = -1

#============================
# functions used
#============================
def get_data(direct, sex, segment, sup_or_jux):
    os.chdir(direct)
    y = []
    fname = sex + '_rat_' + segment + '_pressure_in_Lumen'+sup_or_jux+'.txt'
    file = open(fname, 'r')
    y.extend(np.loadtxt(fname, delimiter = '\n', unpack = True))
    file.close()
    os.chdir('..')
    return y

#==================
# get data
#==================
print('sup nephron start')
sup_data = []
for seg in range(len(sup_list)):
    segment = sup_list[seg]
    print(segment)
    sup_data.extend(get_data(direct, sex, segment, '_sup'))

print('jux nephrons start')
jux1_data = []
jux2_data = []
jux3_data = []
jux4_data = []
jux5_data = []
for seg in range(len(jux_list)):
    segment = jux_list[seg]
    print(segment)
    jux1_data.extend(get_data(direct, sex, segment, '_jux1'))
    jux2_data.extend(get_data(direct,sex,segment, '_jux2'))
    jux3_data.extend(get_data(direct,sex,segment, '_jux3'))
    jux4_data.extend(get_data(direct,sex,segment, '_jux4'))
    jux5_data.extend(get_data(direct,sex,segment, '_jux5'))
    
print('cd start')
cd_data = []
for seg in range(len(segs_cd)):
    segment = segs_cd[seg]
    print(segment)
    cd_data.extend(get_data(direct,sex,segment,''))
    
print('sup CNT end LP: ' + str(sup_data[-1]))
print('jux1 CNT end LP: '+str(jux1_data[-1])) 
print('jux2 CNT end LP: '+str(jux2_data[-1]))
print('jux3 CNT end LP: '+str(jux3_data[-1]))   
print('jux4 CNT end LP: '+str(jux4_data[-1]))
print('jux5 CNT end LP: '+str(jux5_data[-1]))
print('cd end LP: '+str(cd_data[-1]))
#==================================
# position/length of each segment
#==================================
#if sex == 'male':
lens_sup = {'PT': 1.1, 'SDL': 0.14, 'mTAL': 0.2, 'cTAL': 0.2,
            'DCT': 0.1, 'CNT': 0.2}
# jux lengths that are different from sup nephron
lens_jux = {'LDL': 0.5, 'LAL': 0.5, 'cTAL': 0.05, 'CNT': 0.3}
lens_cd = {'CCD':0.2, 'OMCD':0.2, 'IMCD':0.425}
# else:
#     print('sex: '+ sex)
#     raise Exception('sex segment lengths not set up yet')

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


sup_early = ax.plot(sup_early_pos, sup_data[0:400], color = c1, label = 'sup')
sup_late = ax.plot(sup_late_pos, sup_data[401:1400], color = c1)
jux_early= ax.plot(jux_pos_early, jux5_data[0:1200], color = c2, linestyle = 'dashed', label = 'jux5')
jux_late = ax.plot(jux_pos_late, jux5_data[1200:1600], color = c2, linestyle = 'dashed')
cd = ax.plot(cd_pos[0:len(cd_data)], cd_data, color = c3, label = 'cd')

plt.yticks(fontsize=yticklab_size)
plt.xticks(fontsize=xticklab_size)
plt.title(direct + ' luminal pressure', fontsize = title_size)
plt.legend(fontsize = leg_size)
ax = plt.gca()
ax.set_ylabel('luminal pressure (mmHg)', fontsize = ylab_size)
ax.set_xlabel('nephron length (cm)', fontsize = xlab_size)

height = sup_data[0]+shift
plt.hlines(height, pos_PT[0], pos_PT[-1],colors = 'k', linestyles = 'dotted')
plt.text(0.4, height+0.1, 'PT', fontsize = text_size)

height = sup_data[200]+shift
plt.hlines(height, pos_SDL[0], pos_SDL[-1], colors = 'k', linestyles = 'dotted')
plt.text(1.0, height+0.15, 'SDL', fontsize = text_size)

height = jux5_data[400]+shift+1.5
plt.hlines(height, pos_LDL[0], pos_LDL[-1], colors = 'k', linestyles = 'dotted')
plt.text(1.35, height+0.15, 'LDL', fontsize = text_size)

height = jux5_data[600]+shift
plt.hlines(height, pos_LAL[0], pos_LAL[-1], colors = 'k', linestyles = 'dotted')
plt.text(1.87, height+0.15, 'LAL')

height = sup_data[400]+shift
plt.hlines(height, pos_mTAL[0], pos_cTAL_sup[-1], colors = 'k', linestyles = 'dotted')
plt.text(2.3, height+0.15, 'TAL')

height = sup_data[600]+shift-1
plt.hlines(height, pos_DCT[0], pos_DCT[-1], colors = 'k', linestyles = 'dotted')
plt.text(2.6, height-0.5, 'DCT')

height = sup_data[800]+shift+1
plt.hlines(height, pos_CNT_jux[0], pos_CNT_jux[-1], colors = 'k', linestyles = 'dotted')
plt.text(2.8, height+0.15, 'CNT')

height = cd_data[200]+shift
plt.hlines(height, cd_pos[0], cd_pos[-1], colors = 'k', linestyles = 'dotted')
plt.text(3.3, height+0.15, 'CD')

# height = cd_data[0]+shift
# plt.hlines(height, pos_CCD[0], pos_CCD[-1], colors = 'k', linestyles = 'dotted')
# plt.text(3.0, height+0.15, 'CCD')

# height = cd_data[200]+shift
# plt.hlines(height, pos_OMCD[0], pos_OMCD[-1], colors = 'k', linestyles = 'dotted')
# plt.text(3.3, height+0.15, 'OMCD')

# height = cd_data[400]+shift
# plt.hlines(height, pos_IMCD[0], pos_IMCD[-1], colors = 'k', linestyles = 'dotted')
# plt.text(3.5, height-0.5, 'IMCD')