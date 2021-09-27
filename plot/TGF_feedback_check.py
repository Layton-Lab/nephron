# -*- coding: utf-8 -*-
"""
This is for looking at the luminal pressure and needed measurements for the TGF
mechanism part of the model using model output
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from defs import*

#======================================
# Settings
#======================================
direct = 'female-tgf3'
sex = 'female'

compare = 0 # set to 1 if want to compare to another model
direct2 = 'female-tgf1'
sex2 = 'female'

neph_segments = ['PT', 'S3', 'SDL', 'mTAL', 'cTAL', 'DCT', 'CNT']
cd_segments = ['CCD', 'OMCD', 'IMCD']

Ncells = 200*np.ones(len(neph_segments)+len(cd_segments))
if neph_segments[0] == 'PT':
    Ncells[0] = 176
    Ncells[1] = 25
elif neph_segments[0] == 'S3':
    Ncells[0] = 25

title = direct + ' luminal pressure'

shift = -1

cd_flag = 1 # set to 1 if want cd information
if cd_flag == 0:
    cd_segments = []

#========================================
#========================================
def get_data(direct, sex, segment, sup_or_jux):
    os.chdir(direct)
    y = []
    fname = sex + '_rat_' + segment + '_pressure_in_Lumen'+sup_or_jux+'.txt'
    file = open(fname, 'r')
    y.extend(np.loadtxt(fname, delimiter = '\n', unpack = True))
    file.close()
    os.chdir('..')
    return y

def get_CNT_Cl(direct, sex, sup_or_jux):
    os.chdir(direct)
    fname = sex + '_rat_cTAL_con_of_Cl_in_Lumen_'+sup_or_jux+'.txt'
    file = open(fname, 'r')
    cl_con = np.loadtxt(fname, delimiter = '\n', unpack = True)[-1]
    file.close()
    os.chdir('..')
    return cl_con

print('nephrons start')
sup_data = []
jux1_data = []
jux2_data = []
jux3_data = []
jux4_data = []
jux5_data = []
for seg in range(len(neph_segments)):
    segment = neph_segments[seg]
    print(segment)
    sup_data.extend(get_data(direct, sex, segment, '_sup'))
    jux1_data.extend(get_data(direct, sex, segment, '_jux1'))
    jux2_data.extend(get_data(direct,sex,segment, '_jux2'))
    jux3_data.extend(get_data(direct,sex,segment, '_jux3'))
    jux4_data.extend(get_data(direct,sex,segment, '_jux4'))
    jux5_data.extend(get_data(direct,sex,segment, '_jux5'))

if cd_flag:
    print('cd start')
    cd_data = []
    for seg in range(len(cd_segments)):
        segment = cd_segments[seg]
        print(segment)
        cd_data.extend(get_data(direct,sex,segment,''))
    
print('sup CNT end LP: ' + str(sup_data[-1]))
print('jux1 CNT end LP: '+str(jux1_data[-1])) 
print('jux2 CNT end LP: '+str(jux2_data[-1]))
print('jux3 CNT end LP: '+str(jux3_data[-1]))   
print('jux4 CNT end LP: '+str(jux4_data[-1]))
print('jux5 CNT end LP: '+str(jux5_data[-1]))

    
if compare:
    print('\n')
    print(direct2 + ' start')
    print('nephrons start')
    sup_data2 = []
    jux1_data2 = []
    jux2_data2 = []
    jux3_data2 = []
    jux4_data2 = []
    jux5_data2 = []
    for seg in range(len(neph_segments)):
        segment = neph_segments[seg]
        print(segment)
        sup_data2.extend(get_data(direct2, sex2, segment, '_sup'))
        jux1_data2.extend(get_data(direct2, sex2, segment, '_jux1'))
        jux2_data2.extend(get_data(direct2,sex2,segment, '_jux2'))
        jux3_data2.extend(get_data(direct2,sex2,segment, '_jux3'))
        jux4_data2.extend(get_data(direct2,sex2,segment, '_jux4'))
        jux5_data2.extend(get_data(direct2,sex2,segment, '_jux5'))
    if cd_flag:
        print('cd start')
        cd_data2 = []
        for seg in range(len(cd_segments)):
            segment = cd_segments[seg]
            print(segment)
            cd_data2.extend(get_data(direct2,sex2,segment,''))



#=================
# make figures
#=================
fig = plt.figure()
plt.rcParams['axes.xmargin']=0
fig.set_figheight(10)
fig.set_figwidth(12)

# colors
c1 = 'blue'
c2 = 'orange'
c3 = 'green'
c4 = 'blueviolet'
c5 = 'red'
c6 = 'black'
c7 = 'hotpink'

# fontsizes
xlab_size = 18
xticklab_size = 20
ylab_size = 16
yticklab_size = 16
title_size = 16
leg_size = 16

# plot nephron luminal presssure
plt.plot(sup_data, label='sup', color = c1)
plt.plot(jux1_data, label = 'jux1', color = c2)
plt.plot(jux2_data,label='jux2', color=c3)
plt.plot(jux3_data, label='jux3', color=c4)
plt.plot(jux4_data, label='jux4', color=c5)
plt.plot(jux5_data, label='jux5', color=c6)

if cd_flag:
    # plot cd luminal pressure
    cd_pos = np.arange(sum(Ncells[0:len(neph_segments)]), len(sup_data)+len(cd_data))
    plt.plot(cd_pos, cd_data, label = 'cd', color=c7)
    print('end of cd luminal pressure: '+str(cd_data[-1]))

if compare:
    # plot nephron luminal presssure
    plt.plot(sup_data2, color = c1, linestyle = 'dashed')
    print(direct2 + ' sup CNT end LP: ' + str(sup_data[-1]))
    plt.plot(jux1_data2, color = c2, linestyle = 'dashed')
    print(direct2 +' jux1 CNT end LP: '+str(jux1_data[-1]))
    plt.plot(jux2_data2, color=c3, linestyle = 'dashed')
    print(direct2+' jux2 CNT end LP: '+str(jux2_data[-1]))
    plt.plot(jux3_data2, color=c4, linestyle = 'dashed')
    print(direct2+' jux3 CNT end LP: '+str(jux3_data[-1]))
    plt.plot(jux4_data2, color=c5, linestyle = 'dashed')
    print(direct2+' jux4 CNT end LP: '+str(jux4_data[-1]))
    plt.plot(jux5_data2, color=c6, linestyle = 'dashed')
    
    if cd_flag:
        # plot cd luminal pressure
        cd_pos2 = np.arange(sum(Ncells[0:len(neph_segments)]), len(sup_data)+len(cd_data2))
        plt.plot(cd_pos2, cd_data2, color = c7, linestyle='dashed')
        print(direct2+' end of cd luminal pressure: '+str(cd_data[-1]))
    
    print(direct2+ ' values are the dashed lines')


start = 0
seg_list = neph_segments + cd_segments
for seg in range(len(seg_list)):
    end = Ncells[seg] + start
    if int(start) < len(sup_data):
        height = sup_data[int(start)] + shift
    else:
        height = cd_data[int(start)-len(sup_data)] + shift
    plt.hlines(height, start, end, color = 'k', linestyle = 'dotted')
    temp = start + (end - start)/2.0
    plt.text(temp, height, seg_list[seg])
    start = end
plt.yticks(fontsize=yticklab_size)
plt.title(title, fontsize = title_size)
plt.legend(fontsize=leg_size)
ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.set_ylabel('luminal pressure (mmHg)', fontsize=ylab_size)

#========================
# SNGFR computation
#========================
print('\n')
clMD_sup = get_CNT_Cl(direct, sex, 'sup')
if sex == 'male':
    #SNGFR0_sup = 30
    SNGFR0_sup = 29.6 #tgf4 value 
    Cop_sup = 40
elif sex == 'female':
    SNGFR0_sup = 27.0 #tgf3 value
    Cop_sup = 30
DQ = 10
m_sup = 20
SNGFR_sup = SNGFR0_sup - DQ*np.tanh((clMD_sup - Cop_sup)/m_sup)
print('SNGFR_sup: '+str(SNGFR_sup))

if sex == 'male':
    #SNGFR0_jux = 45
    SNGFR0_jux = 40.0 #tgf4 value
    Cop_jux = 80
elif sex == 'female':
    #SNGFR0_jux = 36
    SNGFR0_jux = 35.7 #tgf3 value
    Cop_jux = 60
m_jux = 40
jux_nephs = ['jux1', 'jux2', 'jux3', 'jux4', 'jux5']
SNGFR_jux_list = np.zeros(len(jux_nephs))
clMD_jux_list = np.zeros(len(jux_nephs))
for j in range(len(jux_nephs)):
    jux = jux_nephs[j]
    print(jux)
    clMD_jux = get_CNT_Cl(direct,sex,jux)
    clMD_jux_list[j] = clMD_jux
    SNGFR_jux = SNGFR0_jux - DQ*np.tanh((clMD_jux-Cop_jux)/m_jux)
    print(jux + ' SNGFR_jux: ' + str(SNGFR_jux))
    SNGFR_jux_list[j] = SNGFR_jux
print('\n')
print('mean(SNGFR_jux_list): '+str(np.mean(SNGFR_jux_list)))

# conversion factors
vol_lumen = 0.004
Vref = 1e-4 
cw=Vref*60e6

SNGFR = cw*vol_lumen

SNGFR_try = 22.5
vol_lumen_try = SNGFR_try/cw