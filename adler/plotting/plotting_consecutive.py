import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from math import ceil
import multiprocessing as mp
from functools import partial 

from adler.data_managing_functions import download_data


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
def plot_consecutiveness_(data_folder,save_folder,tuple_):
   
    (alpha,D,number),dataset = tuple_[0],tuple_[1]
    omega =  dataset.omega.unique()[0]
    delta = np.round(alpha/omega,4)  
    mean_trains_cons,std_trains_cons,total_pulses,isolated_pulses = load_consecutive_statistics(dataset,data_folder)

    colors = ['r','g', 'b']
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

    ax1 = plt.subplot(gs_main[1,0])
     
    ax1.plot(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons, linewidth=0.5, marker = "." , markersize=7, alpha=1)
    ax1.fill_between(np.arange(1,len(mean_trains_cons)+1),mean_trains_cons-std_trains_cons,mean_trains_cons+std_trains_cons,alpha = 0.2)
   
    #X_lim = [0,50]
    #ax1.set_xlim(X_lim);
    ax1.set_yscale('log')
    #ax1.set_ylim(YC_lim)
    ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
    ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
    #ax1.set_xticks([0,3,6,9,12,15])        
    
    
    ax2 = plt.subplot(gs_main[1,1])
    arr = [total_pulses,isolated_pulses,[t-i for t,i in zip(total_pulses,isolated_pulses)]]
    
    X1 = [np.ones(len(arr[i]))*(i+1) for i in range(0,len(arr))]
    bp1 = ax2.boxplot(arr,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

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
        ax2.scatter(xA+X1[i],arr[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

    ax2.tick_params(axis='x', labelsize=8,length=2); 
    ax2.tick_params(axis='y', labelsize=8,length=2)
    ax2.set_xlabel('total,isolated,consecutive',fontsize=8)
    ax2.set_ylabel('counts',fontsize=8)
    #ax2.set_ylim([-1,200])
    ax2.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
    ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
  

    plt.savefig(save_folder+ 'consecutiveness_'+str(delta)+'_'+str(D)+'.pdf', format='pdf')
    
def plot_consecutiveness(description_file,data_folder,save_folder):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha','D','number'])
    plot_consecutiveness__ = partial(plot_consecutiveness_,data_folder,save_folder)
    pool.map(plot_consecutiveness__,tuple_)
    pool.close()
    pool.join()
    return (2)
