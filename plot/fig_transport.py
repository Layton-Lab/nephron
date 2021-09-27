import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

#========================================
# user input
#========================================
compare = 1#2, 3
save_figs = 0 # if want to save figs set to 1

solute_list = ['Na','K','Cl','HCO3','urea','NH4','TA', 'Volume']
#solute_list = ['Na']
#solute_list = ['K']

direct1 = 'latepregnant_rat2021-08-11'
sex1 = 'female'

direct2 = 'female-multi2021-08-11'
sex2 = 'female'

direct3 = '2021-07-23latepregnant_rat'
sex3 = 'female'

label1 = direct1
label2 = direct2
label3 = direct3

segment_early = ['PT', 'S3', 'SDL', 'mTAL', 'cTAL', 'DCT', 'CNT']
segment_jux = ['SDL', 'LDL', 'LAL']
segment_cd = ['CCD','OMCD','IMCD']



humOrrat = 'hum' # set to 'hum' for human model

# conversion factors
if humOrrat == 'rat':
    sup_ratio = 2.0/3.0
elif humOrrat == 'hum':
    sup_ratio = 0.85
    
jux_ratio = 1-sup_ratio
neph_weight = [sup_ratio, 0.4*jux_ratio, 0.3*jux_ratio, 0.15*jux_ratio, 0.1*jux_ratio, 0.05*jux_ratio ]

if humOrrat == 'rat':
    neph_per_kidney = 36000 #number of nephrons per kidney
elif humOrrat == 'hum':
    neph_per_kidney = 1000000
p_to_mu = 1e-6 #convert pmol to micromole
cf = neph_per_kidney * p_to_mu

#==========================================================================
# save figures/comments options (note: requires save_figs == 1)
#==========================================================================
if save_figs:
    plot_folder = input('where to save plots? ')
    comments = input('any comments? ')
    
    if os.path.isdir(plot_folder) == False:
        os.makedirs(plot_folder)
        
    comments_file = open('./'+plot_folder+'/comments.txt', 'w')
    comments_file.write(comments)
    comments_file.close()

#========================================
# functions used
#========================================

def get_transport(direct, sex, solute, segment, supOrjux):
    os.chdir(direct)
    if solute == 'TA':
        fname1 = sex+'_'+humOrrat+'_'+segment+'_flow_of_H2PO4_in_Lumen'+supOrjux+'.txt'
        fname2 = sex+'_'+humOrrat+'_'+segment+'_flow_of_HPO4_in_Lumen'+supOrjux+'.txt'
        file1 = open(fname1, 'r')
        file2 = open(fname2, 'r')
        H2PO4_start = float(file1.readline())
        HPO4_start = float(file2.readline())
        TA_start = (10**(7.4-6.8) * H2PO4_start - HPO4_start)/(1 + 10**(7.4-6.8))
        H2PO4_end = float(np.loadtxt(fname1, delimiter = '\n', unpack = True)[-1])
        HPO4_end = float(np.loadtxt(fname2, delimiter = '\n', unpack = True)[-1])
        TA_end = (10**(7.4-6.8) * H2PO4_end - HPO4_end)/(1 + 10**(7.4-6.8))
        file1.close()
        file2.close()
        transport = TA_start - TA_end
    elif solute == 'Volume':
        fname = sex+'_'+humOrrat+'_'+segment+'_water_volume_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        start = float(file.readline())
        end = float(np.loadtxt(fname, delimiter = '\n', unpack = True)[-1])
        file.close()
        transport = start - end
    else:
        fname = sex+'_'+humOrrat+'_'+segment+'_flow_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        start = float(file.readline())
        end = float(np.loadtxt(fname, delimiter = '\n', unpack = True)[-1])
        file.close()
        transport = start - end
    os.chdir('..')
    transport_cf = transport * cf
    return transport_cf

def get_data(direct, sex, solute, segments):
    # get transport data for nephron segments (i.e., not collecting duct)
    sup_trans = np.zeros(len(segments))
    jux1_trans = np.zeros(len(segments))
    jux2_trans = np.zeros(len(segments))
    jux3_trans = np.zeros(len(segments))
    jux4_trans = np.zeros(len(segments))
    jux5_trans = np.zeros(len(segments))
    
    for s in range(len(segments)):
        seg = segments[s]
        if seg[-2:].lower() == 'cd':
            print('segment: ' + seg)
            raise Exception('not for collecting duct, use get_cd_data for cd')
        if seg.lower() != 'ldl' and seg.lower() != 'lal':
            sup_trans[s] = get_transport(direct, sex, solute, seg, '_sup')
        jux1_trans[s] = get_transport(direct, sex, solute, seg, '_jux1')
        jux2_trans[s] = get_transport(direct, sex, solute, seg, '_jux2')
        jux3_trans[s] = get_transport(direct, sex, solute, seg, '_jux3')
        jux4_trans[s] = get_transport(direct, sex, solute, seg, '_jux4')
        jux5_trans[s] = get_transport(direct, sex, solute, seg, '_jux5')
        
    trans_vals_weighted = np.matrix([sup_trans*neph_weight[0], jux1_trans*neph_weight[1],
                                     jux2_trans*neph_weight[2], jux3_trans*neph_weight[3],
                                     jux4_trans*neph_weight[4], jux5_trans*neph_weight[5]])
        
    sup_vals = neph_weight[0]*sup_trans
    jux_vals = neph_weight[1]*jux1_trans + neph_weight[2]*jux2_trans + \
        neph_weight[3]*jux3_trans + neph_weight[4]*jux4_trans + neph_weight[5]*jux5_trans
    
    return sup_vals, jux_vals, trans_vals_weighted

def get_cd_data(direct, sex, solute, segments):
    trans = np.zeros(len(segments))
    
    for s in range(len(segments)):
        seg = segments[s]
        if seg[-2:].lower() != 'cd':
            print('segment: ' + seg)
            raise Exception('only for collecting duct segments')
        trans[s] = get_transport(direct, sex, solute, seg, '')
    trans = trans
    return trans
#=============================================================
# retrieve data
#==============================================================
for solute in solute_list:
    print(solute)
    sup_vals1, jux_vals1, trans_vals1 = get_data(direct1, sex1, solute, segment_early)
    if compare>1:
        sup_vals2, jux_vals2, trans_vals2 = get_data(direct2, sex2, solute, segment_early)
    if compare > 2:
        sup_vals3, jux_vals3, trans_vals3 = get_data(direct3, sex3, solute, segment_early)
    
    jux_dl_vals1 = get_data(direct1, sex1, solute, segment_jux)[1]
    if compare>1:
        jux_dl_vals2 = get_data(direct2, sex2, solute, segment_jux)[1]
    if compare > 2:
        jux_dl_vals3 = get_data(direct3, sex3, solute, segment_jux)[1]
    
    cd_vals1 = get_cd_data(direct1, sex1, solute, segment_cd)
    if compare>1:
        cd_vals2 = get_cd_data(direct2, sex2, solute, segment_cd)
    if compare > 2:
        cd_vals3 = get_cd_data(direct3, sex3, solute, segment_cd)
    
    #print relevant values 
    print(direct1)
    print('sup vals: ' + str(sup_vals1))
    print('jux vals: ' + str(jux_vals1))
    print('jux dl vals: '+str(jux_dl_vals1))
    print('cd vals: '+str(cd_vals1))
    print('\n')
    
    if compare>1:
        print(direct2)
        print('sup vals: ' + str(sup_vals2))
        print('jux vals: ' + str(jux_vals2))
        print('jux dl vals: '+str(jux_dl_vals2))
        print('cd vals: '+str(cd_vals2))
        print('\n')
    
    if compare > 2:
        print(direct3)
        print('sup vals: ' + str(sup_vals3))
        print('jux vals: ' + str(jux_vals3))
        print('jux dl vals: '+str(jux_dl_vals3))
        print('cd vals: '+str(cd_vals3))
        print('\n')
    #==================================================
    # compute total transport for the reported values
    #=================================================
    # these are segments and labels we want
    segment_transport = ['PT','DL','LAL','TAL','DCT','CNT','CD']
    
    # cd considered separately
    sup_final1 = np.zeros(len(segment_transport)-1)
    if compare>1:
        sup_final2 = np.zeros(len(segment_transport)-1)
    if compare > 2:
        sup_final3 = np.zeros(len(segment_transport)-1)
    
    jux_final1 = np.zeros(len(segment_transport)-1)
    if compare > 1:
        jux_final2 = np.zeros(len(segment_transport)-1)
    if compare > 2:
        jux_final3 = np.zeros(len(segment_transport)-1)
    
    # PT = pct & s3
    sup_final1[0] = sum(sup_vals1[0:1+1]) 
    if compare>1:
        sup_final2[0] = sum(sup_vals2[0:1+1])
    if compare>2:
        sup_final3[0] = sum(sup_vals3[0:1+1])
    
    jux_final1[0] = sum(jux_vals1[0:1+1])
    if compare>1:
        jux_final2[0] = sum(jux_vals2[0:1+1])
    if compare > 2:
        jux_final3[0] = sum(jux_vals3[0:1+1])
    
    # DL = sdl for sup, sdl + ldl for jux
    sup_final1[1] = sup_vals1[2]
    if compare>1:
        sup_final2[1] = sup_vals2[2]
    if compare > 2:
        sup_final3[1] = sup_vals3[2]
    #sdl + ldl
    jux_final1[1] = sum(jux_dl_vals1[0:1+1])
    if compare>1:
        jux_final2[1] = sum(jux_dl_vals2[0:1+1])
    if compare > 2:
        jux_final3[1] = sum(jux_dl_vals3[0:1+1])
    
    # LAL
    jux_final1[2] = jux_dl_vals1[-1]
    if compare>1:
        jux_final2[2] = jux_dl_vals2[-1]
    if compare > 2:
        jux_final3[2] = jux_dl_vals3[-1]
    
    # TAL = mTAL + cTAL
    sup_final1[3] = sum(sup_vals1[3:4+1])
    if compare>1:
        sup_final2[3] = sum(sup_vals2[3:4+1])
    if compare > 2:
        sup_final3[3] = sum(sup_vals3[3:4+1])
    
    jux_final1[3] = sum(jux_vals1[3:4+1])
    if compare>1:
        jux_final2[3] = sum(jux_vals2[3:4+1])
    if compare > 2:
        jux_final3[3] = sum(jux_vals3[3:4+1])
    
    # DCT, CNT
    sup_final1[4:5+1] = sup_vals1[5:6+1]
    if compare>1:
        sup_final2[4:5+1] = sup_vals2[5:6+1]
    if compare > 2:
        sup_final3[4:5+1] = sup_vals3[5:6+1]
    
    jux_final1[4:5+1] = jux_vals1[5:6+1]
    if compare>1:
        jux_final2[4:5+1] = jux_vals2[5:6+1]
    if compare > 2:
        jux_final3[4:5+1] = jux_vals3[5:6+1]
    
    # CD = ccd + omcd + imcd
    cd_final1 = sum(cd_vals1)
    if compare>1:
        cd_final2 = sum(cd_vals2)
    if compare > 2:
        cd_final3 = sum(cd_vals3)
    
    #=================================================================
    # make figure
    #==================================================================
    
    fig, ax = plt.subplots()
    # fig settings
    fig.set_figheight(10)
    fig.set_figwidth(12)
    
    bar_width = 0.25
    
    # colors
    c1 = 'c'
    c2 = 'mediumvioletred'
    c3 = 'green'
    
    # fontsizes
    xlab_size = 18
    xticklab_size = 20
    ylab_size = 22
    yticklab_size = 20
    title_size = 25
    leg_size = 18
    
    # positions
    neph_pos = np.arange(len(segment_transport)-1)
    cd_pos = np.arange(len(segment_transport)-1, len(segment_transport))
    
    # bar1
    sup1 = ax.bar(neph_pos, sup_final1, bar_width, align = 'center', edgecolor = 'black',
                  color = c1, label = label1)
    jux1 = ax.bar(neph_pos, jux_final1, bar_width, bottom = sup_final1, align = 'center', 
                  edgecolor = 'black', color = 'white')
    cd1 = ax.bar(cd_pos, cd_final1, bar_width, align = 'center', edgecolor = 'black',
                 color = c1)
    
    if compare>1:
        # bar2
        sup2 = ax.bar(neph_pos+bar_width, sup_final2, bar_width, align = 'center', edgecolor = 'black',
                      color = c2, label=label2)
        jux2 = ax.bar(neph_pos+bar_width, jux_final2, bar_width, bottom = sup_final2, align = 'center',
                      edgecolor = 'black', color = 'white')
        cd2 = ax.bar(cd_pos+bar_width, cd_final2, bar_width, align = 'center', edgecolor = 'black',
                     color = c2)
    
    if compare > 2:
        # bar3
        sup3 = ax.bar(neph_pos + 2*bar_width, sup_final3, bar_width, align = 'center', edgecolor = 'black',
                      color = c3, label=label3)
        jux3 = ax.bar(neph_pos + 2*bar_width, jux_final3, bar_width, bottom = sup_final3, align='center',
                      edgecolor='black', color='white')
        cd3 = ax.bar(cd_pos + 2*bar_width, cd_final3, bar_width, align ='center', edgecolor = 'black',
                     color = c3)
    
    
    ax.set_xticks(np.arange(len(segment_transport))+1*bar_width)
    ax.set_xticklabels(segment_transport, fontsize=xticklab_size)
    ax.legend(fontsize=leg_size)
    plt.yticks(fontsize=yticklab_size)
    plt.axhline(0, color = 'k')
    if solute == 'Volume':
        ax.set_ylabel('volume transport (ml/min)', fontsize = ylab_size)
    else:
        ax.set_ylabel(solute + ' transport ($\mu$mol/min)', fontsize = ylab_size)
    ax.set_title(solute + ' transport', fontsize = title_size)
    
    if save_figs:
        plt.savefig('./'+plot_folder+'/'+solute+' transport', bbox_inches = 'tight')