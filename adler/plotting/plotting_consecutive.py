import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from math import ceil
import multiprocessing as mp
from functools import partial 

from adler.data_managing_functions import download_data,check_file,time
from adler.plotting.plotting_main import set_scale,mask_arr,load_activity



def mean_consecutive_value(trials):
    '''
    agarra una lista de trials de consecutives y te los devuelve en mean y std
    '''
    if len(trials) > 0:
        arr_aux = []
        for j in range(np.max([len(i) for i in trials])):
            arr_aux.append([i[j] for i in trials if len(i) > j])
        return (np.array([np.mean(k) for k in arr_aux]),np.array([np.std(k) for k in arr_aux]))
    else:
        return ([],[])



def load_consecutive_statistics(dataset,data_folder):
        ''' le pasas un experimento  y te devuelve la estadistica de pulsos cons'''
        isolated_pulses_dataset = []
        total_pulses_dataset = []
        consecutive_trains_dataset = []
        
        for (order,row) in dataset.groupby(['order']):
            number      = int(row.number)
            file_name   =  str(number)+'_'+str(order)+'.pkl'
            isolated_pulses_dataset.append(download_data(data_folder+'i_'+file_name))
            
            consecutive_trial = download_data(data_folder+'c_'+file_name)
            consecutive_trains_dataset.append(consecutive_trial)
            total_pulses_dataset.append(consecutive_trial[0])
        
        return (mean_consecutive_value(consecutive_trains_dataset),total_pulses_dataset,isolated_pulses_dataset)
        
#%%
        
def plot_consecutiveness_activity(dt,beg,T,d,N,Delta,description_file,data_folder,save_folder):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['omega', 'alpha','D','number'])
    plot_consecutiveness_activity__ = partial(plot_consecutiveness_activity_,dt,T,d,data_folder,save_folder)
    pool.map(plot_consecutiveness_activity__,tuple_)
    
    plot_time_series_square_dataset_ = partial(plot_time_series_square_dataset,dt,beg,T,d,N,Delta,data_folder,save_folder)
    pool.map(plot_time_series_square_dataset_,tuple_)
    
    pool.close()
    pool.join()
    return (2)

def plot_consecutiveness_activity_(dt,T,d,data_folder,save_folder,tuple_):

    
# =============================================================================
#     consecutiveness plot
# =============================================================================
    (omega,alpha,D,number),dataset = tuple_[0],tuple_[1]
    delta = np.round(alpha/omega,4)  
    (mean_trains_cons,std_trains_cons),total_pulses,isolated_pulses = load_consecutive_statistics(dataset,data_folder)

    colors = ['r','g', 'b']
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

    ax3 = plt.subplot(gs_main[1,0])
     
    ax3.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , markersize=7, alpha=1)
    ax3.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,alpha = 0.2)
   
    #X_lim = [0,50]
    #ax3.set_xlim(X_lim);
    ax3.set_yscale('log')
    #ax3.set_ylim(YC_lim)
    ax3.set_ylabel('counts',fontsize=10); ax3.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax3.xaxis.set_label_coords(0.5, -0.08);ax3.yaxis.set_label_coords(-0.2,0.5);
    #ax3.set_xticks([0,3,6,9,12,15])        

# =============================================================================
#     consecutiveness boxplot
# =============================================================================
    
    ax4 = plt.subplot(gs_main[1,1])
    arr = [total_pulses,isolated_pulses,[t-i for t,i in zip(total_pulses,isolated_pulses)]]
    
    X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
    bp1 = ax4.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp1['boxes']):
         box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp1['whiskers']):
        whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp1['caps']):
        cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp1['medians']):
        median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp1['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians
    
    for i in range(len(X1)):
        xA = np.random.normal(0, 0.1, len(arr[i])), 
        ax4.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

    ax4.tick_params(axis='x', labelsize=8,length=2); 
    ax4.tick_params(axis='y', labelsize=8,length=2)
    ax4.set_xlabel('total,isolated,consecutive',fontsize=8)
    ax4.set_ylabel('counts',fontsize=8)
    #ax4.set_ylim([-1,200])
    ax4.xaxis.set_label_coords(0.5, -0.12);ax4.yaxis.set_label_coords(-0.05,0.5)
    ax4.tick_params(labelsize=6,direction='out', pad=1,length=2)


# =============================================================================
#     activity population and mean plot
# =============================================================================
    ax2 = plt.subplot(gs_main[0,1]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    activity,silent,n_cell = load_activity(dataset,data_folder,dt,T,d)
    
    #population activity
    if len(activity) > 0:
        p1 = ax2.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0)
        p2 = ax2.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=1,alpha=0.8,linewidth=0.0)
        
    ax2.set_xlim([0,n_cell]);ax2.set_ylim([0,100])
    ax2.set_xlabel( ' trazas ',fontsize=8); 
    ax2.set_xticks([1,n_cell + 1])
    ax2.set_yticks([0,50,100])
    ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax2.xaxis.set_label_coords(0.5,-0.06)
    
    #mean activity
    ax1 = plt.subplot(gs_main[0,0]); plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    
    
    if len(activity) > 0:
        p1 = ax2.barh(np.arange(n_cell),width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
        p2 = ax2.barh(np.arange(n_cell),width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),alpha=0.8,linewidth=0.0,height=0.6)

    ax2.set_xticks([0,50,100])
    ax2.set_ylim([0,100])
    ax2.set_yticks([1,n_cell + 1])
    ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
    ax2.invert_yaxis()
    
    
    ax1.set_ylabel('fraction of cell track' ,fontsize=8); 
    ax2.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_folder+ 'consecutiveness_activity_'+str(delta)+'_'+str(D)+'.pdf', format='pdf')
    


def plot_time_series_square_dataset(dt,beg,T,d,N,Delta,data_folder,save_path_name,tuple_):
    '''
    
    '''
    (omega,alpha,D,number),dataset = tuple_[0],tuple_[1]
    delta = np.round(alpha/omega,4)  
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.05,2.05] ;         

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(10, 7, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    axs = axs.ravel(); 
        
    for ix, (order,data)  in  enumerate(dataset.groupby(['order'])):
        if ix < len(axs):
            ax = axs[ix]
            print(delta,D)
            
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
                
                ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2)
            
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
    

    
    plt.savefig(save_path_name + 'time_series'+str(delta)+'_'+str(D)+'.pdf', format='pdf')
    return(0)

#%%
    

    