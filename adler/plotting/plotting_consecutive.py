import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from math import ceil
import multiprocessing as mp
from functools import partial 
import seaborn as sns

from adler.data_managing_functions import download_data,check_file,time
from adler.plotting.plotting_main import set_scale,mask_arr,load_activity,compute_st_values,download_quantifiers,load_activity_dist
from adler.plotting.dyncode_main import get_consecutive_data_dyncode,get_exp_N_total_isolated_consecutive,get_activity_data_dyncode,get_conc_data

#%%
import matplotlib
matplotlib.rcParams['lines.markeredgecolor'] = 'black'
matplotlib.rcParams['lines.markeredgewidth'] = 0
matplotlib.rcParams['savefig.transparent'] = True
matplotlib.rc('pdf', fonttype=42)
matplotlib.rcParams['savefig.transparent'] = True


# change font
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


sns.despine()
sns.set(context='paper', style='ticks')
plt.grid(0)
plt.rc('axes.spines', top=True, bottom=True, left=True, right=True); 

def von_mises(theta,k=10):
    return np.exp(k*np.sin(theta))

#%%
def get_mean_value_place(trials,no_std = False):
    #por cada lista vacia, es como si hubiera una lista de ceros! 
    arr_aux = []
    if len([len(i) for i in trials]) > 0: 
        for j in range(np.max([len(i) for i in trials])):
            aux = []
           
            for i in trials:
                if len(i) > j : 
                    aux.append(i[j])
                else:
                    aux.append(0)
                
            arr_aux.append(aux)
    if no_std:
        return np.array([np.mean(k) for k in arr_aux])
    else:
        return(np.array([np.mean(k) for k in arr_aux]),np.array([np.std(k) for k in arr_aux]))


def load_consecutive_statistics_realizations(dataset,save_data_arr,T):
    mean_trains_cons_trials,total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials = [],[],[],[]
    
    for data_folder in save_data_arr:
        #para cada trial
        mean_trains_cons,total_pulses,isolated_pulses,consecutive_pulses = load_consecutive_statistics(dataset,data_folder,T)
        mean_trains_cons_trials.append(mean_trains_cons)
        total_pulses_trials.append(total_pulses/T)
        isolated_pulses_trials.append(isolated_pulses/T)
        consecutive_pulses_trials.append(consecutive_pulses/T)
        
    return get_mean_value_place(mean_trains_cons_trials,False),total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials


def load_consecutive_statistics(dataset,data_folder,T):
        ''' le pasas un experimento  y te devuelve la estadistica de pulsos cons
        
        FALTAAA VER CONSECUTIVE TRIAL

        para total, consec, isolated te devuelve la mediana del pulse rate (para un solo trial que le pasaste)
        
        '''
        isolated_pulses_dataset = []
        total_pulses_dataset = []
        consecutive_pulses_dataset = []
        consecutive_trains_dataset = []
        n_cell = 0
        
        for (order,row) in dataset.groupby(['order']):
            if n_cell < 68:
                number      = int(row.number)
                file_name   =  str(number)+'_'+str(order)+'.pkl'
                n_cell = n_cell+1
                
                if (check_file('i_'+file_name,data_folder)): 
                    #esto sucede solamente si hay pulsos en la serie temporal
                    
                    isolated_pulses = download_data(data_folder+'i_'+file_name)
                    isolated_pulses_dataset.append(isolated_pulses)
                
                    consecutive_trial = download_data(data_folder+'c_'+file_name) #esto es sin normalizar
                    #consecutive_trains_dataset.append(consecutive_trial)
                    total_pulses_dataset.append(consecutive_trial[0])
                    
                    consecutive_pulses = consecutive_trial[0]-isolated_pulses
                    #print(consecutive_trial[0],isolated_pulses,consecutive_pulses)
                    consecutive_pulses_dataset.append(consecutive_pulses)
                    
                    
                else:
                    isolated_pulses_dataset.append(0)
                    total_pulses_dataset.append(0)
                    consecutive_pulses_dataset.append(0)
                    
                if (check_file('exp_c_'+file_name,data_folder)):
                    #esto sucede sólo si hay máximos
                    consecutive_trial_exp = download_data(data_folder+'exp_c_'+file_name)
                    consecutive_trains_dataset.append(consecutive_trial_exp)
                else:
                    consecutive_trains_dataset.append([0])
                    #para get_mean_value_place, es lo mismo lista vacía que lista con cerois
        return get_mean_value_place(consecutive_trains_dataset,True),np.median(total_pulses_dataset),np.median(isolated_pulses_dataset),np.median(consecutive_pulses_dataset)

# =============================================================================
# Esto es para dist
# =============================================================================


def load_consecutive_statistics_realizations_dist(dataset,save_data_arr,T):
    mean_trains_cons_trials,total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials = [],[],[],[]
    
    for data_folder in save_data_arr:
        mean_trains_cons,total_pulses,isolated_pulses,consecutive_pulses = load_consecutive_statistics_dist(dataset,data_folder,T)
        
        mean_trains_cons_trials.append(mean_trains_cons)
        total_pulses_trials.append(total_pulses/T)
        isolated_pulses_trials.append(isolated_pulses/T)
        consecutive_pulses_trials.append(consecutive_pulses/T)
        
    return get_mean_value_place(mean_trains_cons_trials,False),total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials


def load_consecutive_statistics_dist(ref_,data_folder,T):
        ''' le pasas un experimento  y te devuelve la estadistica de pulsos cons
        Con distinto number! y order puede ser (osea distintas alphas)'''
        
        isolated_pulses_dataset = []
        total_pulses_dataset = []
        consecutive_pulses_dataset = []
        consecutive_trains_dataset = []
        cell_n = 0
        for (number,order),row in ref_.groupby(['number','order']):
            file_name   =  str(number)+'_'+str(order)+'.pkl'
            
            if cell_n < 69:
                if (check_file('i_'+file_name,data_folder)): 
                    cell_n = cell_n + 1
                    isolated_pulses = download_data(data_folder+'i_'+file_name)
                    isolated_pulses_dataset.append(isolated_pulses)
                
                    consecutive_trial = download_data(data_folder+'c_'+file_name)
                    #consecutive_trains_dataset.append(consecutive_trial)
                    total_pulses_dataset.append(consecutive_trial[0])
                
                    consecutive_pulses = consecutive_trial[0]-isolated_pulses
                    consecutive_pulses_dataset.append(consecutive_pulses)
                else:
                    isolated_pulses_dataset.append(0)
                    total_pulses_dataset.append(0)
                    consecutive_pulses_dataset.append(0)
                    
                if (check_file('exp_c_'+file_name,data_folder)):
                    consecutive_trial_exp = download_data(data_folder+'exp_c_'+file_name)
                    consecutive_trains_dataset.append(consecutive_trial_exp)
                else:
                    consecutive_trains_dataset.append([0])
                    
        #print('len total_pulses_dataset' , len(total_pulses_dataset))
        return get_mean_value_place(consecutive_trains_dataset,True),np.median(total_pulses_dataset),np.median(isolated_pulses_dataset),np.median(consecutive_pulses_dataset)
        #return get_mean_value_place(consecutive_trains_dataset,True),sum(total_pulses_dataset),sum(isolated_pulses_dataset),sum(consecutive_pulses_dataset)


#%%
        
def plot_consecutiveness_activity(dt,beg,T,d,N,Delta,description_file,data_folder,save_folder,dyncode_filename,save_data_arr):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['omega', 'alpha','D','number'])
    plot_consecutiveness_activity__ = partial(plot_consecutiveness_activity_,dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr)
    pool.map(plot_consecutiveness_activity__,tuple_)
    
    plot_time_series_square_dataset_ = partial(plot_time_series_square_dataset,dt,beg,T,d,N,Delta,data_folder,save_folder)
    pool.map(plot_time_series_square_dataset_,tuple_)
    
    pool.close()
    pool.join()
    return (2)

def plot_consecutiveness_activity_(dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr,tuple_):
    
    green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]

    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    (omega,alpha,D,number),dataset = tuple_[0],tuple_[1]
    delta = np.round(alpha/omega,4)  
# =============================================================================
#     quantifiers hist plot
# =============================================================================
    gs_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main[0:2])
    dyncode_df = get_conc_data(dyncode_filename)['an_WT_ESL']
    DT,IPI,joint_duration,dm,pulse_rate = download_quantifiers(dataset,data_folder,T,dt,d,False)
    
    ax1 = plt.subplot(gs_row_1[1,0])
    ax1_dc = plt.subplot(gs_row_1[0,0])
    bins_dc = ax1_dc.hist(dyncode_df.dt_peaks.dropna().values/3,bins=np.linspace(0,20,21),density=True,color=green,alpha=1,linewidth=0); 
    print("dyncode")
    compute_st_values(ax1_dc,dyncode_df.dt_peaks.dropna().values/3,bins_dc,1,10)   
    
    if len(DT) > 0:
        ax1.axvspan(6, 8.33, color=green, alpha=0.3, lw=0)
        bins = ax1.hist(DT,bins=np.linspace(0,20,21),density=True,alpha=1,linewidth=0); 
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        print("TS")
        compute_st_values(ax1,DT,bins,1,10)   
    else:
        print(delta,D,"no data")
    
    for ax in [ax1,ax1_dc]:    
        ax.set_ylim([0,0.2]);
        ax.set_xlim([0,20])
        set_scale(ax,[0,5,10,15,20], [0,0.2])
        ax.set_yticklabels([0,0.2])
        ax.tick_params(labelsize=10)
    ax1.set_xticklabels([0,5,10,15,20])
    ax1_dc.set_xticklabels([])
    ax1.set_xlabel('duración' ,fontsize=10); 
    
    ax2 = plt.subplot(gs_row_1[1,1])
    ax2_dc = plt.subplot(gs_row_1[0,1])
    bins_dc = ax2_dc.hist(dyncode_df.IPI.dropna().values/3,np.linspace(0,40,21),density=True,color = green,alpha=1,linewidth=0); 
    print("dyncode")
    compute_st_values(ax2_dc,dyncode_df.IPI.dropna().values/3,bins_dc,1,10)   
    if len(DT) > 0:
        ax2.axvspan(8, 18.67, color=green, alpha=0.3, lw=0)
        bins = ax2.hist(IPI,bins=np.linspace(0,40,21),density=True,alpha=1,linewidth=0); 
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        print("TS")
        compute_st_values(ax2,IPI,bins,1,10)   
    else:
        print(delta,D,"no data")
    
    for ax in [ax2,ax2_dc]:
        ax.set_ylim([0,0.1]);
        ax.set_xlim([0,40])
        set_scale(ax,[0,10,20,30,40], [0,0.1])
        ax.set_yticklabels([0,0.1])
        ax.tick_params(labelsize=10)
    ax2.set_xticklabels([0,10,20,30,40])
    ax2_dc.set_xticklabels([])
    ax2.set_xlabel('IPI' ,fontsize=10); 

    ax3 = plt.subplot(gs_row_1[1,2])
    ax3_dc = plt.subplot(gs_row_1[0,2])
    dyncode_pr = dyncode_df.groupby(level="cell").amp_peaks.count()/(dyncode_df.groupby(level="cell").FRAME.count()/3)
    bins_dc = ax3_dc.hist(dyncode_pr,bins=np.linspace(0,0.08,10),density=True,color = green,alpha=1,linewidth=0); 
    print("dyncode")
    compute_st_values(ax3_dc,dyncode_pr,bins_dc,1,10)   
    if len(DT) > 0:
        ax3.axvspan( 0.02, 0.06, color=green, alpha=0.3, lw=0)
        bins = ax3.hist(pulse_rate,bins=np.linspace(0,0.08,10),density=True,alpha=1,linewidth=0); 
        #ax3.axvspan(0.0067, 0.02, color=green, alpha=0.3, lw=0) #old,creo que es en frames
        print("TS")
        compute_st_values(ax3,pulse_rate,bins,1,10)   
    else:
        print(delta,D,"no data")
    
    for ax in [ax3_dc,ax3]:
        ax.set_ylim([0,40]);
        ax.set_xlim([0,0.08])
        set_scale(ax,[0,0.02,0.04,0.06,0.08], [0,40])
        ax.set_yticklabels([0,40])
        ax.tick_params(labelsize=10)    

    ax3.set_xticklabels([0,0.02,0.04,0.06,0.08])
    ax3_dc.set_xticklabels([])
    ax3.set_xlabel('tasa de pulsado' ,fontsize=10); 
    
# =============================================================================
#     consecutiveness plot
# =============================================================================
    gs_row_3 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[3])

    (mean_trains_cons,std_trains_cons),total_pulses_median,isolated_pulses_median,consecutive_pulses_median = load_consecutive_statistics_realizations(dataset,save_data_arr,T)
    #print('total,is,con',total_pulses_median,isolated_pulses_median,consecutive_pulses_median )
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    ax5 = plt.subplot(gs_row_3[0])
     
    ax5.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , markersize=7, alpha=1,label = "model")
    ax5.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,alpha = 0.2,linewidth=0)
    x,y = get_consecutive_data_dyncode(dyncode_filename)
    ax5.plot(x,y,linewidth=0.5, marker = "." , markersize=7, alpha=1,color = green,label = "serum + LIF")

    ax5.set_xlim([0,10.5]);
    ax5.set_yscale('log')
    ax5.set_ylim([10**(-2),10**1])
    ax5.set_ylabel('counts x trace',fontsize=10); ax5.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax5.xaxis.set_label_coords(0.5, -0.08);ax5.yaxis.set_label_coords(-0.2,0.5);
    ax5.set_xticks([1,4,7,10])  
    ax5.legend()      

# =============================================================================
#     consecutiveness boxplot
# =============================================================================
    
    ax6 = plt.subplot(gs_row_3[1])
    ax6.axhline(y = 1,color = green,linewidth=0.5,alpha = 1,linestyle = 'dashed')

    #hay que hacer varios trials para tener este plot
    
    total_median,isolated_median,consecutive_median = get_exp_N_total_isolated_consecutive(dyncode_filename) 
    #print( 'dyncode',total_median,isolated_median,consecutive_median)

    total_pulses_normed = [i/total_median for i in total_pulses_median]
    isolated_pulses_normed = [i/isolated_median for i in isolated_pulses_median]
    consecutive_pulses_normed = [i/consecutive_median for i in consecutive_pulses_median]
   # print( total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed)
    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
    #print('len total pulses normed' , len(total_pulses_normed))
    
    X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
    bp1 = ax6.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp1['boxes']):
         box_.set( color=colors[0], linewidth=0.0,facecolor=colors[0],alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp1['whiskers']):
        whisker.set(color=colors[0],linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp1['caps']):
        cap.set(color=colors[0],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp1['medians']):
        median.set(color=colors[0],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp1['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians
    
    for i in range(len(X1)):
        xA = np.random.normal(0, 0.1, len(arr[i])), 
        ax6.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

#    total_N,isolated_N,consecutive_N = get_exp_N_total_isolated_consecutive(dyncode_filename) 
#    total_pulses_normed = sum([i/total_N for i in total_pulses])
#    isolated_pulses_normed = sum([i/isolated_N for i in isolated_pulses])
#    consecutive_pulses = [t-i for t,i in zip(total_pulses,isolated_pulses)]
#    consecutive_pulses_normed = sum([i/consecutive_N for i in consecutive_pulses])
#    
#    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
#    print(total_pulses_normed,total_pulses,total_N)
    
#    ax6.plot([1,2,3],arr,' o ')
    ax6.tick_params(axis='x', labelsize=8,length=2); 
    ax6.tick_params(axis='y', labelsize=8,length=2)
    ax6.set_xlabel(' ',fontsize=8)
    ax6.set_ylabel('cuentas normalizadas',fontsize=8)
    ax6.set_ylim([0.0,2.5])
    ax6.xaxis.set_label_coords(0.5, -0.12);ax6.yaxis.set_label_coords(-0.05,0.5)
    ax6.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax6.set_xticklabels([' total' ,'isolated','consecutive'],rotation = 0)


# =============================================================================
#     activity population and mean plot
# =============================================================================
    gs_row_2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[2])

    ax3 = plt.subplot(gs_row_2[0:2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    activity,silent,n_cell = load_activity(dataset,data_folder,dt,T,d)
    
    #population activity
    if len(activity) > 0:
        p1 = ax3.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0)
        p2 = ax3.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=0.9,alpha=0.8,linewidth=0.0)
        x , y , silent_experiment = get_activity_data_dyncode(dyncode_filename)
        p3 = ax3.bar(x,y,bottom=silent_experiment,width=0.9,alpha=0.3,linewidth=0.0,color = green)
        
    ax3.set_xlim([0,n_cell]);ax3.set_ylim([0,100])
    ax3.set_xlabel( ' trazas ',fontsize=8); 
    ax3.set_xticks([1,n_cell + 1])
    ax3.set_yticks([0,50,100])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax3.xaxis.set_label_coords(0.5,-0.06)
    
    #mean activity
    ax4 = plt.subplot(gs_row_2[2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    
    if len(activity) > 0:
        p1 = ax4.barh(1,width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
        p2 = ax4.barh(1,width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),alpha=0.8,linewidth=0.0,height=0.6)
        p3 = ax4.barh(1,width = np.mean(x),left=np.mean(silent_experiment),alpha=0.3,linewidth=0.0,height=0.6,color = green)

    ax4.set_xticks([0,50,100])
    ax4.set_xlim([0,100])
    #ax4.set_yticks([1,n_cell + 1])
    ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax4.invert_yaxis()
    
    
    ax3.set_ylabel('fraction of cell track' ,fontsize=8); 
    ax4.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_folder+ 'consecutiveness_activity_'+str(delta)+'_'+str(D)+'.pdf', format='pdf')
    


def plot_time_series_square_dataset(dt,beg,T,d,N,Delta,data_folder,save_path_name,tuple_):
    '''
    
    '''
    (omega,alpha,D,number),dataset = tuple_[0],tuple_[1]
    delta = np.round(alpha/omega,4)  
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ;# ylim = [-0.2,2.2] ;   
    ylim = [0,22500] ;        

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(10, 7, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    axs = axs.ravel(); 
        
    for ix, (order,data)  in  enumerate(dataset.groupby(['order'])):
        if ix < len(axs):
            ax = axs[ix]
          #  print(delta,D)
            
            ################################################
            #### download data
            ################################################
            file_name =  str(number)+'_'+str(order)+'.pkl'
            if check_file(file_name,data_folder):            
                    
                theta = download_data(data_folder + file_name) 
                t = time(dt,T,d)
                end = len(t)
                beg_ = int(beg/(dt*d))
                assert len(t) == len(theta), (len(t),len(theta))
                
                #ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2)
                ax.plot(t[beg_:end:Delta],von_mises(theta)[beg_:end:Delta],linewidth=2)
                
            if check_file('max_'+file_name,data_folder):            
                    
                MAX          = mask_arr(beg_,end, download_data(data_folder + 'max_'+file_name))
                left_minima  = mask_arr(beg_,end, download_data(data_folder + 'left_minima_'+ file_name) )
                right_minima = mask_arr(beg_,end, download_data(data_folder + 'right_minima_'+ file_name) )
                
                if len(MAX) > 0:
                   # ax.plot(t[beg_:end][MAX],(1+ np.sin(theta))[beg_:end][MAX],'o',color = 'red',markersize = 8)
                   # ax.plot(t[beg_:end][left_minima],(1+ np.sin(theta))[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                   # ax.plot(t[beg_:end][right_minima],(1+ np.sin(theta))[beg_:end][right_minima],'>',color='black',markersize = 8)
                    ax.plot(t[beg_:end][MAX],von_mises(theta)[beg_:end][MAX],'o',color = 'red',markersize = 8)
                    ax.plot(t[beg_:end][left_minima],von_mises(theta)[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                    ax.plot(t[beg_:end][right_minima],von_mises(theta)[beg_:end][right_minima],'>',color='black',markersize = 8)


            ax.set_ylim(ylim);
            ax.set_xlim(xlim)
            
            ###############################################
            #### Plotting
            ################################################
            text = str(ix)
            ax.text(0.9, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            
            if (ix == len(axs)-6):
                ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                ax.set_xlabel('time (min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
            
            #set_scale(ax,[beg,T], [0,2]) 
            set_scale(ax,[beg,T], [0,22500])
            ax.set_xticklabels([beg,T])
            #ax.set_yticklabels([0,2])
            ax.set_yticklabels([0,22500])
            ax.tick_params(labelsize=20)
    

    
    plt.savefig(save_path_name + 'time_series'+str(delta)+'_'+str(D)+'.pdf', format='pdf')
    return(0)

#%%
# =============================================================================
# for dist daata!!
# =============================================================================
#%%
        
def plot_consecutiveness_activity_dist(dt,beg,T,d,N,Delta,description_file,data_folder,save_folder,dyncode_filename,save_data_arr):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['omega','D'])
    plot_consecutiveness_activity_dist__ = partial(plot_consecutiveness_activity_dist_,dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr)
    pool.map(plot_consecutiveness_activity_dist__,tuple_)
    
    plot_time_series_square_dataset_dist_ = partial(plot_time_series_square_dataset_dist,dt,beg,T,d,N,Delta,data_folder,save_folder)
    pool.map(plot_time_series_square_dataset_dist_,tuple_)
    
    pool.close()
    pool.join()
    return (2)

def plot_consecutiveness_activity_dist_(dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr,tuple_):
    
    green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]
    orange = sns.color_palette("deep",6)[1]
    
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    (omega,D),ref_ = tuple_[0],tuple_[1]

# =============================================================================
#     quantifiers hist plot
# =============================================================================
    dyncode_df = get_conc_data(dyncode_filename)['an_WT_ESL']
    gs_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main[0:2])
    
    
    DT,IPI,pulse_rate= [],[],[]
    for (alpha,number),dataset in ref_.groupby(['alpha','number']):
         DT_aux,IPI_aux,joint_duration_aux,dm_aux,pulse_rate_aux = download_quantifiers(dataset,data_folder,T,dt,d,False)
         DT= DT + list(DT_aux); IPI = IPI + list(IPI_aux); pulse_rate = pulse_rate + list(pulse_rate_aux)
         
    ax1 = plt.subplot(gs_row_1[1,0])
    ax1_dc = plt.subplot(gs_row_1[0,0])
    bins_dc = ax1_dc.hist(dyncode_df.dt_peaks.dropna().values/3,bins=np.linspace(0,20,21),density=True,color=green,alpha=1,linewidth=0); 
    print("dyncode")
    compute_st_values(ax1_dc,dyncode_df.dt_peaks.dropna().values/3,bins_dc,1,10)   
    
    
    if len(DT) > 0:
        ax1.axvspan(6, 8.33, color=green, alpha=0.3, lw=0)
        bins = ax1.hist(DT,bins=np.linspace(0,20,21),density=True,color = orange,alpha=1,linewidth=0); 
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        compute_st_values(ax1,DT,bins,1,10)
    else:
        pass
  
    
    for ax in [ax1,ax1_dc]:
        ax.set_ylim([0,0.2]);
        ax.set_xlim([0,20])
        set_scale(ax,[0,5,10,15,20], [0,0.2])
        ax.set_yticklabels([0,0.2])
        ax.tick_params(labelsize=10)
    ax1.set_xticklabels([0,5,10,15,20])
    ax1_dc.set_xticklabels([])
    ax1.set_xlabel('duración' ,fontsize=10);
    
    ax2 = plt.subplot(gs_row_1[1,1])
    ax2_dc = plt.subplot(gs_row_1[0,1])
    bins_dc = ax2_dc.hist(dyncode_df.IPI.dropna().values/3,np.linspace(0,40,21),density=True,color = green,alpha=1,linewidth=0); 
    compute_st_values(ax2_dc,dyncode_df.IPI.dropna().values/3,bins_dc,1,10)   
    
    
    if len(DT) > 0:
        ax2.axvspan(8, 18.67, color=green, alpha=0.3, lw=0)
        bins = ax2.hist(IPI,bins=np.linspace(0,40,21),density=True,color = orange,alpha=1,linewidth=0); 
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        compute_st_values(ax2,IPI,bins,1,10)   
    else:
        pass
    
    for ax in [ax2,ax2_dc]:
        ax.set_ylim([0,0.1]);
        ax.set_xlim([0,40])
        set_scale(ax,[0,10,20,30,40], [0,0.1])
        ax.set_yticklabels([0,0.1])
        ax.tick_params(labelsize=10)
    ax2.set_xticklabels([0,10,20,30,40])
    ax2_dc.set_xticklabels([])
    ax2.set_xlabel('IPI' ,fontsize=10); 

    ax3 = plt.subplot(gs_row_1[1,2])
    ax3_dc = plt.subplot(gs_row_1[0,2])
    dyncode_pr = dyncode_df.groupby(level="cell").amp_peaks.count()/(dyncode_df.groupby(level="cell").FRAME.count()/3)
    bins_dc = ax3_dc.hist(dyncode_pr,bins=np.linspace(0,0.08,10),density=True,color = green,alpha=1,linewidth=0); 
    compute_st_values(ax3_dc,dyncode_pr,bins_dc,1,10)   
    
    if len(DT) > 0:
        ax3.axvspan( 0.02, 0.06, color=green, alpha=0.3, lw=0)
        bins = ax3.hist(pulse_rate,bins=np.linspace(0,0.08,10),density=True,color = orange,alpha=1,linewidth=0); 
        #ax3.axvspan(0.0067, 0.02, color=green, alpha=0.3, lw=0) #old,creo que es en frames
        compute_st_values(ax3,pulse_rate,bins,1,10)   
    else:
        pass
    
    
    for ax in [ax3_dc,ax3]:
        ax.set_ylim([0,40]);
        ax.set_xlim([0,0.08])
        set_scale(ax,[0,0.02,0.04,0.06,0.08], [0,40])
        ax.set_yticklabels([0,40])
        ax.tick_params(labelsize=10)    

    ax3.set_xticklabels([0,0.02,0.04,0.06,0.08])
    ax3_dc.set_xticklabels([])
    ax3.set_xlabel('tasa de pulsado' ,fontsize=10); 
# =============================================================================
#     consecutiveness plot
# =============================================================================
    (mean_trains_cons,std_trains_cons),total_pulses,isolated_pulses,consecutive_pulses = load_consecutive_statistics_realizations_dist(ref_,save_data_arr,T)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']


    gs_row_3 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[3])
    ax5 = plt.subplot(gs_row_3[0])
     
    ax5.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, color = orange, marker = "." , markersize=7, alpha=1,label = "model")
    ax5.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,color = orange,alpha = 0.2,linewidth=0)
    x,y = get_consecutive_data_dyncode(dyncode_filename)
    ax5.plot(x,y,linewidth=0.5, marker = "." , markersize=7, alpha=1,color = green,label = "experiment")

    ax5.set_xlim([0,10.5]);
    ax5.set_yscale('log')
    ax5.set_ylim([10**(-2),10**1])
    ax5.set_ylabel('counts x trace',fontsize=10); ax5.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax5.xaxis.set_label_coords(0.5, -0.08);ax5.yaxis.set_label_coords(-0.2,0.5);
    ax5.set_xticks([1,4,7,10])  
    ax5.legend()      
 

# =============================================================================
#     consecutiveness boxplot
# =============================================================================
    
    ax6 = plt.subplot(gs_row_3[1])
    ax6.axhline(y = 1,color = green,linewidth=1,alpha = 1,linestyle = 'dashed')

#hay que hacer varios trials para tener este plot
    
    total_N,isolated_N,consecutive_N = get_exp_N_total_isolated_consecutive(dyncode_filename) 
    

    total_pulses_normed = [i/total_N for i in total_pulses]
    isolated_pulses_normed = [i/isolated_N for i in isolated_pulses]
    consecutive_pulses_normed = [i/consecutive_N for i in consecutive_pulses]
    #print(consecutive_N,consecutive_pulses)
    
    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
    #print(total_pulses_normed,total_pulses,total_N)
    
    X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
    bp1 = ax6.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp1['boxes']):
         box_.set( color=orange, linewidth=0.0,facecolor=orange,alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp1['whiskers']):
        whisker.set(color=orange,linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp1['caps']):
        cap.set(color=orange,linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp1['medians']):
        median.set(color=orange,linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp1['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians
    
    for i in range(len(X1)):
        xA = np.random.normal(0, 0.1, len(arr[i])), 
        ax6.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
        

    ax6.tick_params(axis='x', labelsize=8,length=2); 
    ax6.tick_params(axis='y', labelsize=8,length=2)
    ax6.set_xlabel(' ',fontsize=8)
    ax6.set_ylabel('cuentas normalizadas',fontsize=8)
    ax6.set_ylim([0.0,2.5])
    ax6.xaxis.set_label_coords(0.5, -0.12);ax6.yaxis.set_label_coords(-0.05,0.5)
    ax6.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax6.set_xticklabels([' total' ,'isolated','consecutive'],rotation = 0)

# =============================================================================
#     activity population and mean plot
# =============================================================================
    gs_row_2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[2])
    ax3 = plt.subplot(gs_row_2[0:2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    #activity,silent,n_cell = load_activity(dataset,data_folder,dt,T,d)
    activity,silent,n_cell = load_activity_dist(ref_,data_folder,dt,T,d)
    
    #population activity
    if len(activity) > 0:
        p1 = ax3.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0)
        p2 = ax3.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=0.9,color=orange,alpha=0.8,linewidth=0.0)
        x , y , silent_experiment = get_activity_data_dyncode(dyncode_filename)
        p3 = ax3.bar(x,y,bottom=silent_experiment,width=0.9,alpha=0.3,linewidth=0.0,color = green)
        
    ax3.set_xlim([0,n_cell]);ax3.set_ylim([0,100])
    ax3.set_xlabel( ' trazas ',fontsize=8); 
    ax3.set_xticks([1,n_cell + 1])
    ax3.set_yticks([0,50,100])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax3.xaxis.set_label_coords(0.5,-0.06)
    
    #mean activity
    ax4 = plt.subplot(gs_row_2[2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    
    if len(activity) > 0:
        p1 = ax4.barh(1,width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
        p2 = ax4.barh(1,width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),alpha=0.8,color=orange,linewidth=0.0,height=0.6)
        p3 = ax4.barh(1,width = np.mean(x),left=np.mean(silent_experiment),alpha=0.3,linewidth=0.0,height=0.6,color = green)

    ax4.set_xticks([0,50,100])
    ax4.set_xlim([0,100])
    #ax4.set_yticks([1,n_cell + 1])
    ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax4.invert_yaxis()
    
    
    ax3.set_ylabel('fraction of cell track' ,fontsize=8); 
    ax4.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_folder+ 'consecutiveness_activity_'+str(D)+'.pdf', format='pdf')

    




def plot_time_series_square_dataset_dist(dt,beg,T,d,N,Delta,data_folder,save_path_name,tuple_):
    '''
    
    '''
    orange = sns.color_palette("deep",6)[1]
    (omega,D),ref_ = tuple_[0],tuple_[1]
    #delta = np.round(alpha/omega,4)  
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.2,2.2] ;         

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(10, 7, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    axs = axs.ravel(); 
        
    for ix, ((alpha,number,order),data)  in  enumerate(ref_.groupby(['alpha','number','order'])):
        if ix < len(axs):
            ax = axs[ix]
            delta = np.round(alpha/omega,4) 
           # print(delta,D)
            
            ################################################
            #### download data
            ################################################
            file_name =  str(number)+'_'+str(order)+'.pkl'
            if check_file(file_name,data_folder):            
                    
                theta = download_data(data_folder + file_name) 
                t = time(dt,T,d)
                end = len(t)
                beg_ = int(beg/(dt*d))
                assert len(t) == len(theta), (len(t),len(theta))
                
                ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color=orange)
            
            if check_file('max_'+file_name,data_folder):            
                    
                MAX          = mask_arr(beg_,end, download_data(data_folder + 'max_'+file_name))
                left_minima  = mask_arr(beg_,end, download_data(data_folder + 'left_minima_'+ file_name) )
                right_minima = mask_arr(beg_,end, download_data(data_folder + 'right_minima_'+ file_name) )
                
                if len(MAX) > 0:
                    ax.plot(t[beg_:end][MAX],(1+ np.sin(theta))[beg_:end][MAX],'o',color = 'blue',markersize = 8)
                    ax.plot(t[beg_:end][left_minima],(1+ np.sin(theta))[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                    ax.plot(t[beg_:end][right_minima],(1+ np.sin(theta))[beg_:end][right_minima],'>',color='black',markersize = 8)


            ax.set_ylim(ylim);
            ax.set_xlim(xlim)
            
            ###############################################
            #### Plotting
            ################################################
            text = str(ix)
            ax.text(0.9, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            
            if (ix == len(axs)-6):
                ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                ax.set_xlabel('time (min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
            
            set_scale(ax,[beg,T], [0,2])
            ax.set_xticklabels([beg,T])
            ax.set_yticklabels([0,2])
            ax.tick_params(labelsize=20)
    

    
    plt.savefig(save_path_name + 'time_series_'+str(D)+'.pdf', format='pdf')
    return(0)

#%% FOR OU DATA
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# Plotting ou
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================


def plot_consecutiveness_activity_ou(dt,beg,T,d,N,Delta,description_file,data_folder,save_folder,dyncode_filename,save_data_arr):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['omega', 'alpha0','sigma','tau','number'])
    plot_consecutiveness_activity__ = partial(plot_consecutiveness_activity_ou_,dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr)
    pool.map(plot_consecutiveness_activity__,tuple_)
    
    plot_time_series_square_dataset_ = partial(plot_time_series_square_dataset_ou,dt,beg,T,d,N,Delta,data_folder,save_folder)
    pool.map(plot_time_series_square_dataset_,tuple_)
    
    pool.close()
    pool.join()
    return (2)

def plot_consecutiveness_activity_ou_(dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr,tuple_):
    
    red = sns.color_palette("deep",6)[3]
    green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]

    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    (omega,alpha0,sigma,tau,number),dataset = tuple_[0],tuple_[1]
    delta0 = np.round(alpha0/omega,4)  
# =============================================================================
#     quantifiers hist plot
# =============================================================================
    gs_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main[0:2])
    DT,IPI,joint_duration,dm,pulse_rate = download_quantifiers(dataset,data_folder,T,dt,d,False)
    dyncode_df = get_conc_data(dyncode_filename)['an_WT_ESL']
    
    
    ax1 = plt.subplot(gs_row_1[1,0])
    ax1_dc = plt.subplot(gs_row_1[0,0])
    bins_dc = ax1_dc.hist(dyncode_df.dt_peaks.dropna().values/3,bins=np.linspace(0,20,21),density=True,color=green,alpha=1,linewidth=0); 
    compute_st_values(ax1_dc,dyncode_df.dt_peaks.dropna().values/3,bins_dc,1,10)   
    
    if len(DT) > 0:
        
        bins = ax1.hist(DT,bins=np.linspace(0,20,21),density=True,color = red,alpha=1,linewidth=0); 
        ax1.axvspan(6, 8.33, color=green, alpha=0.3, lw=0)
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        compute_st_values(ax1,DT,bins,1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax1,ax1_dc]:
        ax.set_ylim([0,0.2]);
        ax.set_xlim([0,20])
        set_scale(ax,[0,5,10,15,20], [0,0.2])
        ax.set_yticklabels([0,0.2])
        ax.tick_params(labelsize=10)
    ax1.set_xticklabels([0,5,10,15,20])
    ax1_dc.set_xticklabels([])
    ax1.set_xlabel('duración' ,fontsize=10); 



    ax2 = plt.subplot(gs_row_1[1,1])
    ax2_dc = plt.subplot(gs_row_1[0,1])
    bins_dc = ax2_dc.hist(dyncode_df.IPI.dropna().values/3,np.linspace(0,40,21),density=True,color = green,alpha=1,linewidth=0); 
    compute_st_values(ax2_dc,dyncode_df.IPI.dropna().values/3,bins_dc,1,10)   

    if len(DT) > 0:
        ax2.axvspan(8, 18.67, color=green, alpha=0.3, lw=0)
        bins = ax2.hist(IPI,bins=np.linspace(0,40,21),density=True,color = red,alpha=1,linewidth=0); 
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        compute_st_values(ax2,IPI,bins,1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax2,ax2_dc]:
        ax.set_ylim([0,0.1]);
        ax.set_xlim([0,40])
        set_scale(ax,[0,10,20,30,40], [0,0.1])
        ax.set_yticklabels([0,0.1])
        ax.tick_params(labelsize=10)
    ax2.set_xticklabels([0,10,20,30,40])
    ax2_dc.set_xticklabels([])
    ax2.set_xlabel('IPI' ,fontsize=10); 

    ax3 = plt.subplot(gs_row_1[1,2])
    ax3_dc = plt.subplot(gs_row_1[0,2])
    dyncode_pr = dyncode_df.groupby(level="cell").amp_peaks.count()/(dyncode_df.groupby(level="cell").FRAME.count()/3)
    bins_dc = ax3_dc.hist(dyncode_pr,bins=np.linspace(0,0.08,10),density=True,color = green,alpha=1,linewidth=0); 
    compute_st_values(ax3_dc,dyncode_pr,bins_dc,1,10)   

    if len(DT) > 0:
        ax3.axvspan( 0.02, 0.06, color=green, alpha=0.3, lw=0)
        bins = ax3.hist(pulse_rate,bins=np.linspace(0,0.08,10),density=True,color = red,alpha=1,linewidth=0); 
        #ax3.axvspan(0.0067, 0.02, color=green, alpha=0.3, lw=0) #old,creo que es en frames
        compute_st_values(ax3,pulse_rate,bins,1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax3,ax3_dc]:
        ax.set_ylim([0,40]);
        ax.set_xlim([0,0.08])
        set_scale(ax,[0,0.08], [0,40])       
        ax.set_yticklabels([0,40])
        ax.tick_params(labelsize=10)    
    ax3.set_xticklabels([0,0.08])
    ax3_dc.set_xticklabels([])
    ax3.set_xlabel('pulse rate' ,fontsize=10); 

    
    
# =============================================================================
#     consecutiveness plot
# =============================================================================
    gs_row_3 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[3])

    (mean_trains_cons,std_trains_cons),total_pulses_median,isolated_pulses_median,consecutive_pulses_median = load_consecutive_statistics_realizations(dataset,save_data_arr,T)
    #print('total,is,con',total_pulses_median,isolated_pulses_median,consecutive_pulses_median )
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    ax5 = plt.subplot(gs_row_3[0])
     
    ax5.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , color = red,markersize=7, alpha=1,label = "model")
    ax5.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,color = red,alpha = 0.2,linewidth=0)
    x,y = get_consecutive_data_dyncode(dyncode_filename)
    ax5.plot(x,y,linewidth=0.5, marker = "." , markersize=7, alpha=1,color = green,label = "exp")

    ax5.set_xlim([0,10.5]);
    ax5.set_yscale('log')
    ax5.set_ylim([10**(-2),10**1])
    ax5.set_ylabel('counts x trace',fontsize=10); ax5.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax5.xaxis.set_label_coords(0.5, -0.08);ax5.yaxis.set_label_coords(-0.2,0.5);
    ax5.set_xticks([1,4,7,10])  
    ax5.legend()       

# =============================================================================
#     consecutiveness boxplot
# =============================================================================
    
    ax6 = plt.subplot(gs_row_3[1])
    ax6.axhline(y = 1,color = green,linewidth=0.5,alpha = 1,linestyle = 'dashed')

    #hay que hacer varios trials para tener este plot
    
    total_median,isolated_median,consecutive_median = get_exp_N_total_isolated_consecutive(dyncode_filename) 
    #print( 'dyncode',total_median,isolated_median,consecutive_median)

    total_pulses_normed = [i/total_median for i in total_pulses_median]
    isolated_pulses_normed = [i/isolated_median for i in isolated_pulses_median]
    consecutive_pulses_normed = [i/consecutive_median for i in consecutive_pulses_median]
    #print( total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed)
    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
    #print('len total pulses normed' , len(total_pulses_normed))
    
    X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
    bp1 = ax6.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp1['boxes']):
         box_.set( color=red, linewidth=0.0,facecolor=red,alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp1['whiskers']):
        whisker.set(color=red,linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp1['caps']):
        cap.set(color=red,linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp1['medians']):
        median.set(color=red,linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp1['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians
    
    for i in range(len(X1)):
        xA = np.random.normal(0, 0.1, len(arr[i])), 
        ax6.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

#    total_N,isolated_N,consecutive_N = get_exp_N_total_isolated_consecutive(dyncode_filename) 
#    total_pulses_normed = sum([i/total_N for i in total_pulses])
#    isolated_pulses_normed = sum([i/isolated_N for i in isolated_pulses])
#    consecutive_pulses = [t-i for t,i in zip(total_pulses,isolated_pulses)]
#    consecutive_pulses_normed = sum([i/consecutive_N for i in consecutive_pulses])
#    
#    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
#    print(total_pulses_normed,total_pulses,total_N)
    
#    ax6.plot([1,2,3],arr,' o ')
    ax6.tick_params(axis='x', labelsize=8,length=2); 
    ax6.tick_params(axis='y', labelsize=8,length=2)
    ax6.set_xlabel(' ',fontsize=8)
    ax6.set_ylabel('cuentas normalizadas',fontsize=8)
    ax6.set_ylim([0.0,2.5])
    ax6.xaxis.set_label_coords(0.5, -0.12);ax6.yaxis.set_label_coords(-0.05,0.5)
    ax6.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax6.set_xticklabels([' total' ,'isolated','consecutive'],rotation = 0)


# =============================================================================
#     activity population and mean plot
# =============================================================================
    gs_row_2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[2])

    ax3 = plt.subplot(gs_row_2[0:2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    activity,silent,n_cell = load_activity(dataset,data_folder,dt,T,d)
    
    #population activity
    if len(activity) > 0:
        p1 = ax3.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0)
        p2 = ax3.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=0.9,alpha=0.8,color = red,linewidth=0.0)
        x , y , silent_experiment = get_activity_data_dyncode(dyncode_filename)
        p3 = ax3.bar(x,y,bottom=silent_experiment,width=0.9,alpha=0.3,linewidth=0.0,color = green)
        
    ax3.set_xlim([0,n_cell]);ax3.set_ylim([0,100])
    ax3.set_xlabel( ' trazas ',fontsize=8); 
    ax3.set_xticks([1,n_cell + 1])
    ax3.set_yticks([0,50,100])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax3.xaxis.set_label_coords(0.5,-0.06)
    
    #mean activity
    ax4 = plt.subplot(gs_row_2[2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    
    if len(activity) > 0:
        p1 = ax4.barh(1,width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
        p2 = ax4.barh(1,width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),color = red,alpha=0.8,linewidth=0.0,height=0.6)
        p3 = ax4.barh(1,width = np.mean(x),left=np.mean(silent_experiment),alpha=0.3,linewidth=0.0,height=0.6,color = green)

    ax4.set_xticks([0,50,100])
    ax4.set_xlim([0,100])
    #ax4.set_yticks([1,n_cell + 1])
    ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax4.invert_yaxis()
    
    
    ax3.set_ylabel('fraction of cell track' ,fontsize=8); 
    ax4.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_folder+ 'consecutiveness_activity_'+str(delta0)+'_'+str(sigma)+'_'+str(tau)+'.pdf', format='pdf')
    


def plot_time_series_square_dataset_ou(dt,beg,T,d,N,Delta,data_folder,save_path_name,tuple_):
    '''
    
    '''
    (omega,alpha0,sigma,tau,number),dataset = tuple_[0],tuple_[1]
    delta0 = np.round(alpha0/omega,4)  
    red = sns.color_palette("deep",6)[3]
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.2,2.2] ;         

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(10, 7, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    axs = axs.ravel(); 
        
    for ix, (order,data)  in  enumerate(dataset.groupby(['order'])):
        if ix < len(axs):
            ax = axs[ix]
          #  print(delta,D)
            
            ################################################
            #### download data
            ################################################
            file_name =  str(number)+'_'+str(order)+'.pkl'
            if check_file(file_name,data_folder):            
                    
                theta = download_data(data_folder + file_name) 
                t = time(dt,T,d)
                end = len(t)
                beg_ = int(beg/(dt*d))
                assert len(t) == len(theta), (len(t),len(theta))
                
                ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color = red)
            
            if check_file('max_'+file_name,data_folder):            
                    
                MAX          = mask_arr(beg_,end, download_data(data_folder + 'max_'+file_name))
                left_minima  = mask_arr(beg_,end, download_data(data_folder + 'left_minima_'+ file_name) )
                right_minima = mask_arr(beg_,end, download_data(data_folder + 'right_minima_'+ file_name) )
                
                if len(MAX) > 0:
                    ax.plot(t[beg_:end][MAX],(1+ np.sin(theta))[beg_:end][MAX],'o',color = 'blue',markersize = 8)
                    ax.plot(t[beg_:end][left_minima],(1+ np.sin(theta))[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                    ax.plot(t[beg_:end][right_minima],(1+ np.sin(theta))[beg_:end][right_minima],'>',color='black',markersize = 8)


            ax.set_ylim(ylim);
            ax.set_xlim(xlim)
            
            ###############################################
            #### Plotting
            ################################################
            text = str(ix)
            ax.text(0.9, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            
            if (ix == len(axs)-6):
                ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                ax.set_xlabel('time (min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
            
            set_scale(ax,[beg,T], [0,2])
            ax.set_xticklabels([beg,T])
            ax.set_yticklabels([0,2])
            ax.tick_params(labelsize=20)
    

    
    plt.savefig(save_path_name + 'time_series'+str(delta0)+'_'+str(sigma)+'_'+str(tau)+'.pdf', format='pdf')
    return(0)
    
#%%
#%% FOR OU DATA
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# Plotting ou
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================


def plot_consecutiveness_activity_ou_D(dt,beg,T,d,N,Delta,description_file,data_folder,save_folder,dyncode_filename,save_data_arr):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['omega', 'alpha0','sigma','tau','D','number'])
    plot_consecutiveness_activity__ = partial(plot_consecutiveness_activity_ou_D_,dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr)
    pool.map(plot_consecutiveness_activity__,tuple_)
    
    plot_time_series_square_dataset_ = partial(plot_time_series_square_dataset_ou_D,dt,beg,T,d,N,Delta,data_folder,save_folder)
    pool.map(plot_time_series_square_dataset_,tuple_)
    
    pool.close()
    pool.join()
    return (2)

def plot_consecutiveness_activity_ou_D_(dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr,tuple_):
    
    violet = sns.color_palette("deep",6)[4]
    green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]

    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=6, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    (omega,alpha0,sigma,tau,D,number),dataset = tuple_[0],tuple_[1]
    delta0 = np.round(alpha0/omega,4)  
# =============================================================================
#     quantifiers hist plot
# =============================================================================
    gs_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=3, subplot_spec=gs_main[0:2])
    DT,IPI,joint_duration,dm,pulse_rate = download_quantifiers(dataset,data_folder,T,dt,d,False)
    dyncode_df = get_conc_data(dyncode_filename)['an_WT_ESL']
    
    
    ax1 = plt.subplot(gs_row_1[1,0])
    ax1_dc = plt.subplot(gs_row_1[0,0])
    bins_dc = ax1_dc.hist(dyncode_df.dt_peaks.dropna().values/3,bins=np.linspace(0,20,21),density=True,color=green,alpha=1,linewidth=0); 
    compute_st_values(ax1_dc,dyncode_df.dt_peaks.dropna().values/3,bins_dc,1,10)   
    
    if len(DT) > 0:
        
        bins = ax1.hist(DT,bins=np.linspace(0,20,21),density=True,color = violet,alpha=1,linewidth=0); 
        ax1.axvspan(6, 8.33, color=green, alpha=0.3, lw=0)
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        compute_st_values(ax1,DT,bins,1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax1,ax1_dc]:
        ax.set_ylim([0,0.2]);
        ax.set_xlim([0,20])
        set_scale(ax,[0,5,10,15,20], [0,0.2])
        ax.set_yticklabels([0,0.2])
        ax.tick_params(labelsize=10)
    ax1.set_xticklabels([0,5,10,15,20])
    ax1_dc.set_xticklabels([])
    ax1.set_xlabel('duración' ,fontsize=10); 



    ax2 = plt.subplot(gs_row_1[1,1])
    ax2_dc = plt.subplot(gs_row_1[0,1])
    bins_dc = ax2_dc.hist(dyncode_df.IPI.dropna().values/3,np.linspace(0,40,21),density=True,color = green,alpha=1,linewidth=0); 
    compute_st_values(ax2_dc,dyncode_df.IPI.dropna().values/3,bins_dc,1,10)   

    if len(DT) > 0:
        ax2.axvspan(8, 18.67, color=green, alpha=0.3, lw=0)
        bins = ax2.hist(IPI,bins=np.linspace(0,40,21),density=True,color = violet,alpha=1,linewidth=0); 
        #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
        compute_st_values(ax2,IPI,bins,1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax2,ax2_dc]:
        ax.set_ylim([0,0.1]);
        ax.set_xlim([0,40])
        set_scale(ax,[0,10,20,30,40], [0,0.1])
        ax.set_yticklabels([0,0.1])
        ax.tick_params(labelsize=10)
    ax2.set_xticklabels([0,10,20,30,40])
    ax2_dc.set_xticklabels([])
    ax2.set_xlabel('IPI' ,fontsize=10); 

    ax3 = plt.subplot(gs_row_1[1,2])
    ax3_dc = plt.subplot(gs_row_1[0,2])
    dyncode_pr = dyncode_df.groupby(level="cell").amp_peaks.count()/(dyncode_df.groupby(level="cell").FRAME.count()/3)
    bins_dc = ax3_dc.hist(dyncode_pr,bins=np.linspace(0,0.08,10),density=True,color = green,alpha=1,linewidth=0); 
    compute_st_values(ax3_dc,dyncode_pr,bins_dc,1,10)   

    if len(DT) > 0:
        ax3.axvspan( 0.02, 0.06, color=green, alpha=0.3, lw=0)
        bins = ax3.hist(pulse_rate,bins=np.linspace(0,0.08,10),density=True,color = violet,alpha=1,linewidth=0); 
        #ax3.axvspan(0.0067, 0.02, color=green, alpha=0.3, lw=0) #old,creo que es en frames
        compute_st_values(ax3,pulse_rate,bins,1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax3,ax3_dc]:
        ax.set_ylim([0,40]);
        ax.set_xlim([0,0.08])
        set_scale(ax,[0,0.08], [0,40])       
        ax.set_yticklabels([0,40])
        ax.tick_params(labelsize=10)    
    ax3.set_xticklabels([0,0.08])
    ax3_dc.set_xticklabels([])
    ax3.set_xlabel('pulse rate' ,fontsize=10); 

    
    
# =============================================================================
#     consecutiveness plot
# =============================================================================
    gs_row_3 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[3])

    (mean_trains_cons,std_trains_cons),total_pulses_median,isolated_pulses_median,consecutive_pulses_median = load_consecutive_statistics_realizations(dataset,save_data_arr,T)
    #print('total,is,con',total_pulses_median,isolated_pulses_median,consecutive_pulses_median )

    ax5 = plt.subplot(gs_row_3[0])
     
    ax5.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , color = violet,markersize=7, alpha=1,label = "model")
    ax5.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,color = violet,alpha = 0.2,linewidth=0)
    x,y = get_consecutive_data_dyncode(dyncode_filename)
    ax5.plot(x,y,linewidth=0.5, marker = "." , markersize=7, alpha=1,color = green,label = "exp")

    ax5.set_xlim([0,10.5]);
    ax5.set_yscale('log')
    ax5.set_ylim([10**(-2),10**1])
    ax5.set_ylabel('counts x trace',fontsize=10); ax5.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax5.xaxis.set_label_coords(0.5, -0.08);ax5.yaxis.set_label_coords(-0.2,0.5);
    ax5.set_xticks([1,4,7,10])  
    ax5.legend()       

# =============================================================================
#     consecutiveness boxplot
# =============================================================================
    
    ax6 = plt.subplot(gs_row_3[1])
    ax6.axhline(y = 1,color = green,linewidth=0.5,alpha = 1,linestyle = 'dashed')

    #hay que hacer varios trials para tener este plot
    
    total_median,isolated_median,consecutive_median = get_exp_N_total_isolated_consecutive(dyncode_filename) 
    #print( 'dyncode',total_median,isolated_median,consecutive_median)

    total_pulses_normed = [i/total_median for i in total_pulses_median]
    isolated_pulses_normed = [i/isolated_median for i in isolated_pulses_median]
    consecutive_pulses_normed = [i/consecutive_median for i in consecutive_pulses_median]
    #print( total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed)
    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
    #print('len total pulses normed' , len(total_pulses_normed))
    
    X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
    bp1 = ax6.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp1['boxes']):
         box_.set( color=violet, linewidth=0.0,facecolor=violet,alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp1['whiskers']):
        whisker.set(color=violet,linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp1['caps']):
        cap.set(color=violet,linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp1['medians']):
        median.set(color=violet,linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp1['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians
    
    for i in range(len(X1)):
        xA = np.random.normal(0, 0.1, len(arr[i])), 
        ax6.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

#    total_N,isolated_N,consecutive_N = get_exp_N_total_isolated_consecutive(dyncode_filename) 
#    total_pulses_normed = sum([i/total_N for i in total_pulses])
#    isolated_pulses_normed = sum([i/isolated_N for i in isolated_pulses])
#    consecutive_pulses = [t-i for t,i in zip(total_pulses,isolated_pulses)]
#    consecutive_pulses_normed = sum([i/consecutive_N for i in consecutive_pulses])
#    
#    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
#    print(total_pulses_normed,total_pulses,total_N)
    
#    ax6.plot([1,2,3],arr,' o ')
    ax6.tick_params(axis='x', labelsize=8,length=2); 
    ax6.tick_params(axis='y', labelsize=8,length=2)
    ax6.set_xlabel(' ',fontsize=8)
    ax6.set_ylabel('cuentas normalizadas',fontsize=8)
    ax6.set_ylim([0.0,2.5])
    ax6.xaxis.set_label_coords(0.5, -0.12);ax6.yaxis.set_label_coords(-0.05,0.5)
    ax6.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax6.set_xticklabels([' total' ,'isolated','consecutive'],rotation = 0)


# =============================================================================
#     activity population and mean plot
# =============================================================================
    gs_row_2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_main[2])

    ax3 = plt.subplot(gs_row_2[0:2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    activity,silent,n_cell = load_activity(dataset,data_folder,dt,T,d)
    
    #population activity
    if len(activity) > 0:
        p1 = ax3.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0)
        p2 = ax3.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=0.9,alpha=0.8,color = violet,linewidth=0.0)
        x , y , silent_experiment = get_activity_data_dyncode(dyncode_filename)
        p3 = ax3.bar(x,y,bottom=silent_experiment,width=0.9,alpha=0.3,linewidth=0.0,color = green)
        
    ax3.set_xlim([0,n_cell]);ax3.set_ylim([0,100])
    ax3.set_xlabel( ' trazas ',fontsize=8); 
    ax3.set_xticks([1,n_cell + 1])
    ax3.set_yticks([0,50,100])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax3.xaxis.set_label_coords(0.5,-0.06)
    
    #mean activity
    ax4 = plt.subplot(gs_row_2[2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    
    if len(activity) > 0:
        p1 = ax4.barh(1,width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
        p2 = ax4.barh(1,width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),color = violet,alpha=0.8,linewidth=0.0,height=0.6)
        p3 = ax4.barh(1,width = np.mean(x),left=np.mean(silent_experiment),alpha=0.3,linewidth=0.0,height=0.6,color = green)

    ax4.set_xticks([0,50,100])
    ax4.set_xlim([0,100])
    #ax4.set_yticks([1,n_cell + 1])
    ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax4.invert_yaxis()
    
    
    ax3.set_ylabel('fraction of cell track' ,fontsize=8); 
    ax4.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_folder+ 'consecutiveness_activity_'+str(delta0)+'_'+str(sigma)+'_'+str(tau)+str(D)+'.pdf', format='pdf')
    


def plot_time_series_square_dataset_ou_D(dt,beg,T,d,N,Delta,data_folder,save_path_name,tuple_):
    '''
    
    '''
    (omega,alpha0,sigma,tau,D,number),dataset = tuple_[0],tuple_[1]
    delta0 = np.round(alpha0/omega,4)  
    violet = sns.color_palette("deep",6)[4]
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.2,2.2] ;         

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(10, 7, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    axs = axs.ravel(); 
        
    for ix, (order,data)  in  enumerate(dataset.groupby(['order'])):
        if ix < len(axs):
            ax = axs[ix]
          #  print(delta,D)
            
            ################################################
            #### download data
            ################################################
            file_name =  str(number)+'_'+str(order)+'.pkl'
            if check_file(file_name,data_folder):            
                    
                theta = download_data(data_folder + file_name) 
                t = time(dt,T,d)
                end = len(t)
                beg_ = int(beg/(dt*d))
                assert len(t) == len(theta), (len(t),len(theta))
                
                ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color = violet)
            
            if check_file('max_'+file_name,data_folder):            
                    
                MAX          = mask_arr(beg_,end, download_data(data_folder + 'max_'+file_name))
                left_minima  = mask_arr(beg_,end, download_data(data_folder + 'left_minima_'+ file_name) )
                right_minima = mask_arr(beg_,end, download_data(data_folder + 'right_minima_'+ file_name) )
                
                if len(MAX) > 0:
                    ax.plot(t[beg_:end][MAX],(1+ np.sin(theta))[beg_:end][MAX],'o',color = 'blue',markersize = 8)
                    ax.plot(t[beg_:end][left_minima],(1+ np.sin(theta))[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                    ax.plot(t[beg_:end][right_minima],(1+ np.sin(theta))[beg_:end][right_minima],'>',color='black',markersize = 8)


            ax.set_ylim(ylim);
            ax.set_xlim(xlim)
            
            ###############################################
            #### Plotting
            ################################################
            text = str(ix)
            ax.text(0.9, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            
            if (ix == len(axs)-6):
                ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                ax.set_xlabel('time (min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
            
            set_scale(ax,[beg,T], [0,2])
            ax.set_xticklabels([beg,T])
            ax.set_yticklabels([0,2])
            ax.tick_params(labelsize=20)
    

    
    plt.savefig(save_path_name + 'time_series'+str(delta0)+'_'+str(sigma)+'_'+str(tau)+'_'+str(D)+'.pdf', format='pdf')
    return(0)