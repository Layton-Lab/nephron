import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

compare = 3 #2 or 3
save_figs = 0 # if want to save figs set to 1
#solute_list = ['Na','K','Cl','HCO3','H2CO3','CO2','HPO4','H2PO4','urea','NH3','NH4','H','HCO2','H2CO2','glu', 'TA', 'osmolality']
solute_list = ['Na','K','Cl','HCO3','urea','NH4','TA', 'osmolality', 'pH']

direct1 = 'Female_rat_Non_diab-multiple'
sup_or_jux1 = 'jux5'

direct2 = 'nephron'
sup_or_jux2 = 'jux5'

direct3 = 'nonpregnant_rat-final'
sup_or_jux3 = 'jux5'


sex1 = 'female'
sex2 = 'female'
sex3 = 'nonpregnant'

label1 = direct1
label2 = direct2
label3 = direct3

# segs_sup_early = ['pt', 's3', 'sdl']
# segs_sup_late = ['mtal', 'ctal', 'dct', 'cnt']
# segs_jux = ['pt', 's3', 'sdl', 'ldl','lal','mtal','ctal', 'dct', 'cnt']

segs_sup_early = ['sdl']
segs_sup_late = []
segs_jux = ['sdl']

plot_cd = 0
if plot_cd:
    segs_cd = ['ccd', 'omcd', 'imcd']



humOrrat = 'rat'

if plot_cd:
    segs_total = segs_jux + segs_cd
else:
    segs_total = segs_jux
    
Ncells = 200*np.ones(len(segs_total))
if segs_total[0].upper() == 'PT':
    Ncells[0] = 176
    Ncells[1] = 25
elif segs_total[0].upper() == 'S3':
    Ncells[0] = 25

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


#=======================================================

def get_data(direct, filename):
    os.chdir(direct)
    file = open(filename, 'r')
    y = np.loadtxt(filename, delimiter = '\n', unpack = True)
    file.close()
    os.chdir('..')
    return y

def get_data_TA(direct, sex, segment, supOrjux, humOrrat):
    os.chdir(direct)
    y = []
    fname1 = sex + '_'+ humOrrat +'_' + segment+'_con_of_H2PO4_in_Lumen'+supOrjux+'.txt'
    fname2 = sex + '_'+humOrrat +'_' + segment+'_con_of_HPO4_in_Lumen'+supOrjux+'.txt'
    file1 = open(fname1, 'r')
    file2 = open(fname2, 'r')
    H2PO4_vals = np.loadtxt(fname1, delimiter = '\n', unpack = True)
    HPO4_vals = np.loadtxt(fname2, delimiter = '\n', unpack = True)
    file1.close()
    file2.close()
    TA_vals = (10**(7.4-6.8) * H2PO4_vals - HPO4_vals)/(1 + 10**(7.4-6.8))
    y=TA_vals
    os.chdir('..')
    return y

def get_data_pH(direct, sex, segment, supOrjux, humOrrat):
    fname = sex + '_' + humOrrat+'_'+ segment + '_con_of_H_in_Lumen'+supOrjux+'.txt'
    y = get_data(direct, fname)
    y_ph = [-np.log(i/1000)/np.log(10) for i in y]
    return y_ph

def segment_data(segments, supOrjux, solute):
    y1 = []
    y2 = []
    y3 = []
    for seg in range(len(segments)):
        if solute == 'TA':
            y1.extend(get_data_TA(direct1, sex1, segments[seg], supOrjux, humOrrat))
            y2.extend(get_data_TA(direct2, sex2, segments[seg], supOrjux, humOrrat))
            if compare >2:
                y3.extend(get_data_TA(direct3, sex3, segments[seg], supOrjux, humOrrat))
        elif solute == 'pH':
            y1.extend(get_data_pH(direct1, sex1, segments[seg], supOrjux, humOrrat))
            y2.extend(get_data_pH(direct2, sex2, segments[seg], supOrjux, humOrrat))
            if compare > 2:
                y3.extend(get_data_pH(direct3, sex3, segments[seg], supOrjux, humOrrat))
        else:
            if solute == 'osmolality':
                fname1 = sex1 + '_' + humOrrat+'_'+ segments[seg].upper() + '_osmolality_in_Lumen'+supOrjux+'.txt'
                fname2 = sex2 + '_' + humOrrat+'_'+ segments[seg].upper() + '_osmolality_in_Lumen'+supOrjux+'.txt'
                fname3 = sex3 + '_' + humOrrat+'_'+segments[seg].upper() + '_osmolality_in_Lumen'+supOrjux+'.txt'
            else:
                fname1 = sex1 + '_' + humOrrat+'_'+ segments[seg] + '_con_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
                fname2 = sex2 + '_' + humOrrat+'_'+ segments[seg] + '_con_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
                fname3 = sex3 + '_' + humOrrat+'_'+segments[seg] + '_con_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
            y1.extend(get_data(direct1, fname1))
            y2.extend(get_data(direct2, fname2))
            if compare > 2:
                y3.extend(get_data(direct3, fname3))
    return y1, y2, y3

for sol in range(len(solute_list)):
    solute = solute_list[sol]    
    print('\n' + solute)
    data1_sup_early, data2_sup_early, data3_sup_early = segment_data(segs_sup_early, '_sup', solute)
    data1_sup_late, data2_sup_late, data3_sup_late = segment_data(segs_sup_late, '_sup', solute)
    data1_jux1, data2_jux1, data3_jux1 = segment_data(segs_jux, '_jux1', solute)
    data1_jux2, data2_jux2, data3_jux2 = segment_data(segs_jux, '_jux2', solute)
    data1_jux3, data2_jux3, data3_jux3 = segment_data(segs_jux, '_jux3', solute)
    data1_jux4, data2_jux4, data3_jux4 = segment_data(segs_jux, '_jux4', solute)
    data1_jux5, data2_jux5, data3_jux5 = segment_data(segs_jux, '_jux5', solute)
    if plot_cd:
        data1_cd, data2_cd, data3_cd = segment_data(segs_cd, '', solute)


#======================================================================
# make the figures
#======================================================================

    fig, ax = plt.subplots()
    plt.rcParams['axes.xmargin'] =0
    # fig settings
    fig.set_figheight(12)
    fig.set_figwidth(12)
       
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
    text_size = 18
    
    
    sup_early_pos = np.arange(sum(Ncells[0:len(segs_sup_early)]))
    sup_late_pos = np.arange(sum(Ncells[0:len(segs_sup_early)+2]), sum(Ncells[0:len(segs_jux)]))
    jux_pos= np.arange(sum(Ncells[0:len(segs_jux)]))
    cd_pos = np.arange(sum(Ncells[0:len(segs_jux)]), sum(Ncells))
    
    # shift
    if solute == 'pH':
        shift = -0.1
    else:
        shift = -3
    
    # only plotting the jux5 for the concentrations
    sup1_early = ax.plot(sup_early_pos, data1_sup_early, color = c1, linewidth=2, label = direct1)
    sup1_late = ax.plot(sup_late_pos, data1_sup_late, color = c1, linewidth=2)
    if sup_or_jux1 == 'jux1':
        jux1 = ax.plot(jux_pos,data1_jux1, color = c1, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux1 == 'jux2':
        jux1 = ax.plot(jux_pos,data1_jux2, color = c1, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux1 == 'jux3':
        jux1 = ax.plot(jux_pos,data1_jux3, color = c1, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux1 == 'jux4':
        jux1 = ax.plot(jux_pos,data1_jux4, color = c1, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux1 == 'jux5':
        jux1 = ax.plot(jux_pos,data1_jux5, color = c1, linestyle = 'dashed', linewidth = 2)
    else:
        raise Exception('which juxtamedullary nephron to plot?')
    if plot_cd:
        cd1 = ax.plot(cd_pos, data1_cd, color = c1, linewidth = 2)
    
    sup2_early = ax.plot(sup_early_pos, data2_sup_early, color = c2, linewidth=2, label = direct2)
    sup2_late = ax.plot(sup_late_pos, data2_sup_late, color = c2, linewidth=2)
    if sup_or_jux2 == 'jux1':
        jux2 = ax.plot(jux_pos,data2_jux1, color = c2, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux2 == 'jux2':
        jux2 = ax.plot(jux_pos,data2_jux2, color = c2, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux2 == 'jux3':
        jux2 = ax.plot(jux_pos,data2_jux3, color = c2, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux2 == 'jux4':
        jux2 = ax.plot(jux_pos,data2_jux4, color = c2, linestyle = 'dashed', linewidth = 2)
    elif sup_or_jux2 == 'jux5':
        jux2 = ax.plot(jux_pos,data2_jux5, color = c2, linestyle = 'dashed', linewidth = 2)
    else:
        raise Exception('which juxt nephron to plot?')
    if plot_cd:
        cd2 = ax.plot(cd_pos, data2_cd, color = c2, linewidth = 2)
    
    if compare >2:
        sup3_early = ax.plot(sup_early_pos, data3_sup_early, color = c3, linewidth=2, label = direct3)
        sup3_late = ax.plot(sup_late_pos, data3_sup_late, color = c3, linewidth=2)
        if sup_or_jux3 == 'jux1':
            jux3 = ax.plot(jux_pos,data3_jux1, color = c3, linestyle = 'dashed', linewidth = 2)
        elif sup_or_jux3 == 'jux2':
            jux3 = ax.plot(jux_pos,data3_jux2, color = c3, linestyle = 'dashed', linewidth = 2)
        elif sup_or_jux3 == 'jux3':
            jux3 = ax.plot(jux_pos,data3_jux3, color = c3, linestyle = 'dashed', linewidth = 2)
        elif sup_or_jux3 == 'jux4':
            jux3 = ax.plot(jux_pos,data3_jux4, color = c3, linestyle = 'dashed', linewidth = 2)
        elif sup_or_jux3 == 'jux5':
            jux3 = ax.plot(jux_pos,data3_jux5, color = c3, linestyle = 'dashed', linewidth = 2)
        else:
            raise Exception('which juxt nephron to plot?')
        if plot_cd:
            cd3 = ax.plot(cd_pos, data3_cd, color = c3, linewidth = 2)
    
    start = 0
    if plot_cd:
        jux_full = data1_jux5 + data1_cd
    else:
        jux_full = data1_jux5
    for seg in range(len(segs_total)):
        end = Ncells[seg] + start
        height = jux_full[int(start)] + shift
        plt.hlines(height, start, end, colors = 'k', linestyles = 'dashed')
        temp = start + (end - start)/2.0
        plt.text(temp, height, segs_total[seg].upper(), fontsize = text_size)
        start = end
    
    plt.xticks(fontsize=xticklab_size)
    plt.yticks(fontsize=yticklab_size)
    ax.legend(fontsize=leg_size)
    ax.axes.xaxis.set_ticks([])
    if solute == 'osmolality':
        ax.set_ylabel('osmolality (mosm/kg H$_2$O)', fontsize = ylab_size)
        ax.set_title('osmolality', fontsize = title_size)
    elif solute == 'pH':
        ax.set_ylabel('pH', fontsize = ylab_size)
        ax.set_title('pH', fontsize = title_size)
    else:
        ax.set_ylabel('['+solute+'] (mM)', fontsize = ylab_size)
        ax.set_title(solute + ' concentration', fontsize = title_size)
    
   
    if save_figs:
        if solute == 'pH':
            plt.savefig('./'+plot_folder+'/pH',bbox_inches = 'tight')
        elif solute == 'osmolality':
            plt.savefig('./'+plot_folder+'/osmolality',bbox_inches = 'tight')
        else:
            plt.savefig('./'+plot_folder+'/'+solute+' concentration', bbox_inches = 'tight')
            
print('dashed line is the juxtamedullary concentration, solid is superficial')
            
