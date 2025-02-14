import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from math import ceil
import multiprocessing as mp
from functools import partial 
import seaborn as sns
import os

from adler.data_managing_functions import download_data,check_file,time

from adler.plotting.plotting_main import set_scale,mask_arr,load_activity_realizations,compute_st_values,print_st_values, download_quantifiers_realizations

from adler.plotting.dyncode_main import get_consecutive_data_dyncode,get_exp_N_total_isolated_consecutive_mean,get_activity_data_dyncode,get_conc_data
from adler.plotting.plotting_consecutive import get_mean_value_place
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
#%%
def get_realizations_parameters(main_path,n):
    save_description_arr = [];save_data_arr = []
    for i in range(n):
        save_description_arr.append(main_path + '/descriptions/description_'+str(i)+'.xlsx')
        save_data_arr.append( main_path + '/realizations/data_'+str(i)+'/')
    return save_description_arr,save_data_arr


def load_consecutive_statistics_realizations_mean(dataset,save_data_arr,T):
    mean_trains_cons_trials,total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials = [],[],[],[]
    
    for data_folder in save_data_arr:
        #para cada trial
        mean_trains_cons,total_pulses,isolated_pulses,consecutive_pulses = load_consecutive_statistics_mean(dataset,data_folder,T)
        mean_trains_cons_trials.append(mean_trains_cons)
        total_pulses_trials.append(total_pulses/T)
        isolated_pulses_trials.append(isolated_pulses/T)
        consecutive_pulses_trials.append(consecutive_pulses/T)
        
    return get_mean_value_place(mean_trains_cons_trials,False),total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials


def load_consecutive_statistics_mean(dataset,data_folder,T):
        ''' le pasas un experimento  y te devuelve la estadistica de pulsos cons
        
        FALTAAA VER CONSECUTIVE TRIAL

        para total, consec, isolated te devuelve el promedio del pulse rate (para un solo trial que le pasaste)
        
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
        return get_mean_value_place(consecutive_trains_dataset,True),np.mean(total_pulses_dataset),np.mean(isolated_pulses_dataset),np.mean(consecutive_pulses_dataset)

#%%


def load_consecutive_statistics_realizations_mean_dist(dataset,save_data_arr,T):
    mean_trains_cons_trials,total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials = [],[],[],[]
    
    for data_folder in save_data_arr:
        #para cada trial
        mean_trains_cons,total_pulses,isolated_pulses,consecutive_pulses = load_consecutive_statistics_mean_dist(dataset,data_folder,T)
        mean_trains_cons_trials.append(mean_trains_cons)
        total_pulses_trials.append(total_pulses/T)
        isolated_pulses_trials.append(isolated_pulses/T)
        consecutive_pulses_trials.append(consecutive_pulses/T)
        
    return get_mean_value_place(mean_trains_cons_trials,False),total_pulses_trials,isolated_pulses_trials,consecutive_pulses_trials


def load_consecutive_statistics_mean_dist(dataset,data_folder,T):
        ''' le pasas un experimento  y te devuelve la estadistica de pulsos cons
        
        FALTAAA VER CONSECUTIVE TRIAL

        para total, consec, isolated te devuelve el promedio del pulse rate (para un solo trial que le pasaste)
        
        '''
        isolated_pulses_dataset = []
        total_pulses_dataset = []
        consecutive_pulses_dataset = []
        consecutive_trains_dataset = []
        n_cell = 0
        
        for (order, alpha), row in dataset.groupby(['order', 'alpha']):
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
        return get_mean_value_place(consecutive_trains_dataset,True),np.mean(total_pulses_dataset),np.mean(isolated_pulses_dataset),np.mean(consecutive_pulses_dataset)



#%%

def figure3_m3a2(dt,beg,T,d,N,Delta,description_file,data_folder,save_folder,dyncode_filename,save_data_arr):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['omega', 'alpha0','sigma','tau','D','number'])
    plot_consecutiveness_activity__ = partial(figure3_m3a2_,dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr)
    pool.map(plot_consecutiveness_activity__,tuple_)
    
    plot_time_series_square_dataset_ = partial(ts_figure3_m3a2_,dt,beg,T,d,N,Delta,data_folder,save_folder)
    pool.map(plot_time_series_square_dataset_,tuple_)
    
    pool.close()
    pool.join()
    return (2)

def figure3_m3a2_(dt,T,d,data_folder,save_folder,dyncode_filename,save_data_arr,tuple_):
    
    red = "r"#sns.color_palette("deep",6)[3]
    green =  sns.color_palette(sns.dark_palette("#2ecc71",30,reverse=False))[15]

    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=8, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    
    (omega,alpha0,sigma,tau,D,number),dataset = tuple_[0],tuple_[1]
    delta0 = np.round(alpha0/omega,4)  
    
# =============================================================================
#     quantifiers hist plot
# =============================================================================
    gs_row_1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main[0])
    DT,IPI,joint_duration,dm,pulse_rate =  download_quantifiers_realizations(dataset,save_data_arr,T,dt,d) 
    dyncode_df = get_conc_data(dyncode_filename)['an_WT_ESL']
    
    #dynocde data
    ax1 = plt.subplot(gs_row_1[0,1]) #second column
    #bins_dc = ax1.hist(dyncode_df.dt_peaks.dropna().values/3,bins=np.linspace(0,20,21),density=True,color=green,histtype='step',alpha=1,linewidth=1,hatch='\\\\'); 
    hist_dc, edges_dc = np.histogram(dyncode_df.dt_peaks.dropna().values/3,bins=np.linspace(0,20,21),density=True); 
    bin_centers_dc = (edges_dc[:-1] + edges_dc[1:]) / 2
    ax1.plot(bin_centers_dc, hist_dc, linewidth=0.5, marker = "." , color = green, markersize=7, alpha=1)
    print("duration \n dynode:\n")
    #print_st_values(ax1,dyncode_df.dt_peaks.dropna().values/3,[bin_centers_dc,hist_dc],1,10)   
    
    #synthetic data
    if len(DT) > 0:
        #bins = ax1.hist(DT,bins=np.linspace(0,20,21),density=True,color = red,histtype='step',alpha=1,linewidth=1)#,hatch='///'); 
        #bins = ax1.hist(DT,bins=np.linspace(0,20,21),density=True,color = red,histtype='stepfilled',alpha=0.2,linewidth=0); 
        hist, edges = np.histogram(DT,bins=np.linspace(0,20,21),density=True); 
        bin_centers = (edges[:-1] + edges[1:]) / 2
        ax1.plot(bin_centers, hist, linewidth=0.5, marker = "." , color = red, markersize=7, alpha=1)
        print("duration \n model:\n")
       # print_st_values(ax1,DT,[bin_centers,hist],1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax1]:
        ax.set_ylim([0,0.3]);
        ax.set_xlim([0,20])
        set_scale(ax,[0,5,10,15,20], [0,0.3])
        ax.set_yticklabels([0,0.3])
        ax.tick_params(labelsize=10)
    ax1.set_xticklabels([0,5,10,15,20])
    ax1.set_xlabel('duration (min)' ,fontsize=10); 



    ax2 = plt.subplot(gs_row_1[0,0])
    #bins_dc = ax2.hist(dyncode_df.IPI.dropna().values/3,np.linspace(0,40,21),density=True,color = green,histtype='step',alpha=1,linewidth=1,hatch='\\\\');  
    hist_dc, edges_dc = np.histogram(dyncode_df.IPI.dropna().values/3,np.linspace(0,20,11),density=True); 
    bin_centers_dc = (edges_dc[:-1] + edges_dc[1:]) / 2
    ax2.plot(bin_centers_dc, hist_dc, linewidth=0.5, marker = "." , color = green, markersize=7, alpha=1)
    print("IPI \n dyncode:\n")
    #print_st_values(ax2,dyncode_df.IPI.dropna().values/3,[bin_centers_dc,hist_dc],1,10)   

    if len(DT) > 0:
        #bins = ax2.hist(IPI,bins=np.linspace(0,20,21),density=True,color = red,histtype='step',alpha=1,linewidth=1)#,hatch='///');  
       # bins = ax2.hist(IPI,bins=np.linspace(0,20,21),density=True,color = red,histtype='stepfilled',alpha=0.2,linewidth=0); 
        hist, edges = np.histogram(IPI,bins=np.linspace(0,20,11),density=True); 
        bin_centers = (edges[:-1] + edges[1:]) / 2
        ax2.plot(bin_centers, hist, linewidth=0.5, marker = "." , color = red, markersize=7, alpha=1)
        print("IPI \n model:\n")
        #print_st_values(ax2,IPI,[bin_centers,hist],1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax2]:
        ax.set_ylim([0,0.15]);
        ax.set_xlim([0,20])
        set_scale(ax,[0,5,10,15,20], [0,0.15])
        ax.set_yticklabels([0,0.15])
        ax.tick_params(labelsize=10)
    ax2.set_xticklabels([0,5,10,15,20])
    ax2.set_xlabel('Interpulse interval (min)' ,fontsize=10); 



    ax3 = plt.subplot(gs_row_1[0,2])
    dyncode_pr = dyncode_df.groupby(level="cell").amp_peaks.count()/(dyncode_df.groupby(level="cell").FRAME.count()/3)
   # bins_dc = ax3.hist(dyncode_pr,bins=np.linspace(0,0.08,10),density=True,color = green,histtype='step',alpha=1,linewidth=1,hatch='\\\\'); 
    hist_dc, edges_dc = np.histogram(dyncode_pr,bins=np.linspace(0,0.08,10),density=True); 
    bin_centers_dc = (edges_dc[:-1] + edges_dc[1:]) / 2
    ax3.plot(bin_centers_dc, hist_dc, linewidth=0.5, marker = "." , color = green, markersize=7, alpha=1)
    print("pulse_rate \n dyncode:\n")
    #print_st_values(ax3,dyncode_pr,[bin_centers,hist],1,10)   

    if len(DT) > 0:
        #bins = ax3.hist(pulse_rate,bins=np.linspace(0,0.08,10),density=True,color = red,histtype='step',alpha=1,linewidth=1)#,hatch='///'); 
       # bins = ax3.hist(pulse_rate,bins=np.linspace(0,0.08,10),density=True,color = red,histtype='stepfilled',alpha=0.2,linewidth=0); 
        hist, edges = np.histogram(pulse_rate,bins=np.linspace(0,0.08,10),density=True); 
        bin_centers = (edges[:-1] + edges[1:]) / 2
        ax3.plot(bin_centers, hist, linewidth=0.5, marker = "." , color = red, markersize=7, alpha=1)
        print("pulse_rate \n model:\n")
       # print_st_values(ax3,pulse_rate,[bin_centers,hist],1,10)   
    else:
        print(delta0,"no data")
    
    for ax in [ax3]:
        ax.set_ylim([0,30]);
        ax.set_xlim([0,0.08])
        set_scale(ax,[0,0.08], [0,30])       
        ax.set_yticklabels([0,30])
        ax.tick_params(labelsize=10)    
    ax3.set_xticklabels([0,0.08])
    ax3.set_xlabel('pulse rate (1/min)' ,fontsize=10); 

    
    
# =============================================================================
#     consecutiveness plot
# =============================================================================
    (mean_trains_cons,std_trains_cons),total_pulses_mean,isolated_pulses_mean,consecutive_pulses_mean = load_consecutive_statistics_realizations_mean(dataset,save_data_arr,T)

    ax5 = plt.subplot(gs_row_1[0,3])
     
    ax5.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , color = red,markersize=7, alpha=1,label = "model")
    ax5.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,color = red,alpha = 0.2,linewidth=0)
    x,y = get_consecutive_data_dyncode(dyncode_filename)
    ax5.plot(x,y,linewidth=0.5, marker = "." , markersize=7, alpha=1,color = green,label = "exp")

    ax5.set_xlim([0,8.5]);
    ax5.set_yscale('log')
    ax5.set_ylim([10**(-2),10**1])
    ax5.set_ylabel('counts x trace',fontsize=10); ax5.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax5.xaxis.set_label_coords(0.5, -0.08);ax5.yaxis.set_label_coords(-0.2,0.5);
    ax5.set_xticks([1,4,8])  
    ax5.legend()       


# =============================================================================
#     consecutiveness boxplot
# =============================================================================
    gs_row_2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=4, subplot_spec=gs_main[1])
    
    ax6 = plt.subplot(gs_row_2[0,3])
    ax6.axhline(y = 1,color = green,linewidth=0.5,alpha = 1,linestyle = 'dashed')
    
    total_mean, isolated_mean, consecutive_mean = get_exp_N_total_isolated_consecutive_mean(dyncode_filename) 
    (activity,activity_err),(silent,_),n_cell,activity_mean,silent_mean = load_activity_realizations(dataset,save_data_arr,dt,T,d)
    x , y , silent_experiment = get_activity_data_dyncode(dyncode_filename)


    total_pulses_normed = [i/total_mean for i in total_pulses_mean]
    isolated_pulses_normed = [i/isolated_mean for i in isolated_pulses_mean]
    consecutive_pulses_normed = [i/consecutive_mean for i in consecutive_pulses_mean]
    
    #activity_mean_normed = [i/np.mean(x) for i in activity_mean]

    
    #arr = [activity_mean_normed,total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
    arr = [total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]


    
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


    ax6.tick_params(axis='x', labelsize=8,length=2); 
    ax6.tick_params(axis='y', labelsize=8,length=2)
    ax6.set_xlabel(' ',fontsize=8)
    ax6.set_ylabel('normalized counts',fontsize=8)
    ax6.set_ylim([0.0,1.5])
    ax6.xaxis.set_label_coords(0.5, -0.12);ax6.yaxis.set_label_coords(-0.05,0.5)
    ax6.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax6.set_xticklabels([' total' ,'isolated','consecutive'],rotation = 0)


# =============================================================================
#     activity population and mean plot
# =============================================================================
   # gs_row_2_in = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_row_2[0,1:3])

    ax3 = plt.subplot(gs_row_2[0,1]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
   # (activity,activity_err),(silent,_),n_cell,activity_mean,silent_mean = load_activity_realizations(dataset,save_data_arr,dt,T,d)
    
    #population activity
    if len(activity) > 0:
        #p1 = ax3.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0,yerr=activity_err)
        #p2 = ax3.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=0.9,alpha=0.8,color = red,linewidth=0.0)
        ax3.plot(np.arange(1 ,n_cell + 1),activity[::-1],linewidth=0.5, marker = "." , color = red, markersize=0, alpha=1)
        ax3.fill_between(np.arange(1, n_cell + 1), (activity - activity_err)[::-1], (activity + activity_err)[::-1], color=red, alpha=0.2,linewidth=0)
       # x , y , silent_experiment = get_activity_data_dyncode(dyncode_filename)
        #p3 = ax3.bar(x,y,bottom=silent_experiment,width=0.9,alpha=0.3,linewidth=0.0,color = green)
        ax3.plot(x,y[::-1],linewidth=0.5, marker = "." , color = green, markersize=0, alpha=1)
        
    ax3.set_xlim([0,n_cell]);ax3.set_ylim([0,100])
    ax3.set_xlabel( ' trazas ',fontsize=8); 
    ax3.set_xticks([1,n_cell + 1])
    ax3.set_yticks([0,50,100])
    ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax3.xaxis.set_label_coords(0.5,-0.06)
    

    
    
    
    #mean activity
    gs_row_2_in = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs_row_2[0,0])
    ax4 = plt.subplot(gs_row_2_in[2]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 

    
    if len(activity) > 0:
        #arr = [activity_mean_normed,total_pulses_normed,isolated_pulses_normed,consecutive_pulses_normed]
        arr = [activity_mean]


        
        X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
        bp1 = ax4.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

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
            ax4.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)


        ax4.tick_params(axis='x', labelsize=8,length=2); 
        ax4.tick_params(axis='y', labelsize=8,length=2)
        ax4.set_xlabel(' ',fontsize=8)
        ax4.set_ylabel('proportion of\ncell track',fontsize=8)
        ax4.set_ylim([0,100])
        ax4.xaxis.set_label_coords(0.5, -0.12);ax6.yaxis.set_label_coords(-0.05,0.5)
        ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)
        ax4.set_xticklabels(['mean \n activity'],rotation = 0)    
        
        
    #     p1 = ax4.barh(1,width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
    #     p2 = ax4.barh(1,width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),color = red,alpha=0.8,linewidth=0.0,height=0.6)
    #     p3 = ax4.barh(1,width = np.mean(x),left=np.mean(silent_experiment),alpha=0.3,linewidth=0.0,height=0.6,color = green)

    # ax4.set_xticks([0,50,100])
    # ax4.set_xlim([0,100])
    # #ax4.set_yticks([1,n_cell + 1])
    # ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)
    # ax4.invert_yaxis()
    
    
    # ax3.set_ylabel('fraction of cell track' ,fontsize=8); 
    # ax4.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_folder+ 'figure3_m3a2_'+str(delta0)+'_'+str(sigma)+'_'+str(tau)+str(D)+'.pdf', format='pdf')
    


def ts_figure3_m3a2_(dt,beg,T,d,N,Delta,data_folder,save_path_name,tuple_):
    '''
    
    '''
    (omega,alpha0,sigma,tau,D,number),dataset = tuple_[0],tuple_[1]
    delta0 = np.round(alpha0/omega,4)  
    red = "r"#sns.color_palette("deep",6)[3]

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
                
                ax.plot(t[beg_:end:1],(1+ np.sin(theta))[beg_:end:1],linewidth=2,color = red)

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
            #text = str(ix)
           # ax.text(0.9, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            
            if (ix == len(axs)-6):
                ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                ax.set_xlabel('time (min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
            
            set_scale(ax,[beg,T], [0,1]) #set_scale(ax,[beg,T], [0,2])
            ax.set_xticklabels([beg,T])
            ax.set_yticklabels([0,1])#ax.set_yticklabels([0,2])
            ax.tick_params(labelsize=20)

    
    plt.savefig(save_path_name + 'ts_figure3_m3a2_'+str(delta0)+'_'+str(sigma)+'_'+str(tau)+'_'+str(D)+'.pdf', format='pdf')
    return(0)


#%%

def figure4_m3a2(dt,T,d,mother_path,root_path,params,n):

    params_keys = list(params.keys())
    params_grouped_keys = [[key, f'{key}_m', f'{key}_p'] for key in params_keys]
    colors = ["r","c","b"] 
    labels = ["0","-","+"]
    
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=10, ncols=6, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    
# =============================================================================
#     consecutiveness plot
# =============================================================================
    
    for j,params_group_keys in enumerate(params_grouped_keys): #por cada parametro
        ax5 = plt.subplot(gs_main[j,1])
        ax3 = plt.subplot(gs_main[j,0]); #plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
        
        ax1 = plt.subplot(gs_main[j,2])
        ax2 = plt.subplot(gs_main[j,3])
        ax4 = plt.subplot(gs_main[j,4])
        
        for i,p in enumerate(params_group_keys): #por cada variacion de parametro
            path = mother_path+p+"/" if i != 0 else root_path #esto es donde encontramos realizations  
            
            if os.path.exists(path+"descriptions/"):
                if len(os.listdir(path+"descriptions/"))==100:
                    # para definir dataset
                    ref = pd.read_excel(path+"descriptions/description_0.xlsx",sheet_name= 'File_references');ref.set_index('Unnamed: 0',inplace=True);
                    tuple_ = list(ref.groupby(['omega', 'alpha0','sigma','tau','D','number'])); assert len(tuple_) == 1
                    _ ,dataset = tuple_[0][0],tuple_[0][1]
                 
                    #para definir save_data_arr
                    _ , save_data_arr = get_realizations_parameters(path,n)
                    
                    (mean_trains_cons,std_trains_cons),total_pulses_mean,isolated_pulses_mean,consecutive_pulses_mean = load_consecutive_statistics_realizations_mean(dataset,save_data_arr,T)
                    (activity,activity_err),(silent,_),n_cell,activity_mean,silent_mean = load_activity_realizations(dataset,save_data_arr,dt,T,d)
                    
                    #consecutive
                    ax5.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , color = colors[i],markersize=7, alpha=1,label = labels[i])
                    ax5.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,color = colors[i],alpha = 0.2,linewidth=0)
                
                    #activity
                    if len(activity) > 0:
                        ax3.plot(np.arange(1 ,n_cell + 1),activity[::-1],linewidth=0.5, marker = "." , color = colors[i], markersize=0, alpha=1)
                        ax3.fill_between(np.arange(1, n_cell + 1), (activity - activity_err)[::-1], (activity + activity_err)[::-1], color = colors[i], alpha=0.2,linewidth=0)
                
                
                    # histograms
                    DT,IPI,joint_duration,dm,pulse_rate =  download_quantifiers_realizations(dataset,save_data_arr,T,dt,d) 
                    
                    
                    hist, edges = np.histogram(DT,bins=np.linspace(0,20,21),density=True); 
                    bin_centers = (edges[:-1] + edges[1:]) / 2
                    ax1.plot(bin_centers, hist, linewidth=0.5, marker = "." , color = colors[i], markersize=7, alpha=1)

                    
                    hist, edges = np.histogram(IPI,bins=np.linspace(0,20,11),density=True); 
                    bin_centers = (edges[:-1] + edges[1:]) / 2
                    ax2.plot(bin_centers, hist, linewidth=0.5, marker = "." , color = colors[i], markersize=7, alpha=1)

                
                    hist, edges = np.histogram(pulse_rate,bins=np.linspace(0,0.08,10),density=True); 
                    bin_centers = (edges[:-1] + edges[1:]) / 2
                    ax4.plot(bin_centers, hist, linewidth=0.5, marker = "." , color = colors[i], markersize=7, alpha=1)

                
                else:
                    n_cell = 69
            else:
                n_cell = 69
            
        ax3.set_xlim([0,n_cell]);ax3.set_ylim([0,100])
        ax3.set_xticks([1,n_cell + 1])
        ax3.set_yticks([0,50,100])
        ax3.tick_params(labelsize=6,direction='out', pad=1,length=2)
        ax3.set_ylabel(params_group_keys[0],fontsize=8); 


        ax5.set_xlim([0,8.5]);
        ax5.set_yscale('log')
        ax5.set_ylim([10**(-2),10**1])
        ax5.set_xticks([1,4,8])  
        ax5.set_yticks([10**(-2),10**1])

        for ax in [ax1]:
            ax.set_ylim([0,0.3]);
            ax.set_xlim([0,20])
            set_scale(ax,[0,5,10,15,20], [0,0.3])
            ax.set_yticklabels([0,0.3])
            ax.tick_params(labelsize=10)

        
        for ax in [ax2]:
            ax.set_ylim([0,0.15]);
            ax.set_xlim([0,20])
            set_scale(ax,[0,5,10,15,20], [0,0.15])
            ax.set_yticklabels([0,0.15])
            ax.tick_params(labelsize=10)

        
        for ax in [ax4]:
            ax.set_ylim([0,30]);
            ax.set_xlim([0,0.08])
            set_scale(ax,[0,0.08], [0,30])       
            ax.set_yticklabels([0,30])
            ax.tick_params(labelsize=10)    



    ax5.set_ylabel('counts x trace',fontsize=10); ax5.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax5.xaxis.set_label_coords(0.5, -0.08);ax5.yaxis.set_label_coords(-0.2,0.5);
    ax5.legend()      
    ax3.xaxis.set_label_coords(0.5,-0.06)
    ax3.set_xlabel( ' trazas ',fontsize=8); 
    ax1.set_xticklabels([0,5,10,15,20])
    ax1.set_xlabel('duration (min)' ,fontsize=10); 
    ax2.set_xticklabels([0,5,10,15,20])
    ax2.set_xlabel('Interpulse interval (min)' ,fontsize=10); 
    ax4.set_xticklabels([0,0.08])
    ax4.set_xlabel('pulse rate (1/min)' ,fontsize=10); 


    plt.savefig(root_path+'figures/figure4_m3a2.pdf', format='pdf')