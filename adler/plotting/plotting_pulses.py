import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import ceil
from functools import partial 
import matplotlib.gridspec as gridspec
import seaborn as sns

from plotting.plotting_main import set_scale,download_data,check_file
from pulse_detection.pulse_detection_main import get_fixed_points

from pulse_detection.quantifiers_main import time



#%%

#%%



#%%
###############################################################################
### Function for paralelizing the pulses plotting 
###############################################################################
def mask_arr(end,arr):
    return ([i for i in arr if i  < end])
    

def plot_pulses(description_file,data_folder_ts,data_folder_pulses,save_path_name,dt,T,d,TS,N,delta):
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    plot_pulses_alpha_ = partial(plot_pulses_alpha,data_folder_ts,data_folder_pulses,save_path_name,dt,T,d,TS,N,delta)
    pool.map(plot_pulses_alpha_,ref.groupby(['alpha']) )
    pool.close()
    pool.join()
    return (1)
# =============================================================================
# plotting pulses module
# =============================================================================

def plot_pulses_alpha(data_folder_ts,data_folder_pulses,save_path_name,dt,T,d,TS,N,delta,tuple_):
    '''
    plor pulses for a certain alpha
    data_folder_ts : carpeta donde estan las series temporales
    data_folder_pulses : donde está la data de los pulsos
    save_path_name : donde queres guardar la figura  
    plotea el primer archivo solamente de la serie temporal! (variable j)
    '''
    i,rows = tuple_[0],tuple_[1]; j = 0
###############################################################################
### Parameters
###############################################################################

    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    PFE , PFI = get_fixed_points(alpha)
    print('alpha = ', alpha)
    T_n = ceil(int(T/dt)/d)

###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5,T+5] ; ylim = [-1.1,1.1] ;         
    Cols = TS; Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    colors =  sns.color_palette(sns.color_palette("viridis",Rows*1))
    colors =  colors[::1]

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.1)
    text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' #+' ~ ' + r'$ \Delta\theta = $'+ text_dtheta
    axs[0].text(0,1.5, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=35)
    
    
    for k,(ix,row_) in  enumerate(rows.groupby(['D'])):
        
        row = row_.iloc[j]; order = int(row.order); number = int(row.number)

        file_name =  str(number)+'_'+str(order)+'.pkl'
        file_name_max = 'max_xf_'+file_name

        theta = download_data(data_folder_ts + file_name)[:T_n:delta]
        t = time(dt,T,d)[::delta]
        print(row.D,'plotting ',order,number);D = row.D
        ax = axs[k]; ax.grid(False);
        ax.plot(t,np.cos(theta),linewidth=0.8,color = colors[k])
        
        
        if (check_file(file_name_max,data_folder_pulses)):
                        
            MAX          = mask_arr(T_n, download_data(data_folder_pulses + file_name_max))
            MIN          = mask_arr(T_n, download_data(data_folder_pulses + 'min_xf_'+ file_name))
            left_minima  = mask_arr(T_n, download_data(data_folder_pulses + 'left_minima_'+ file_name) )
            right_minima = mask_arr(T_n, download_data(data_folder_pulses + 'right_minima_'+ file_name) )
            
            ax.plot(t[MAX],np.cos(theta)[MAX],'o',color = 'red',markersize = 8)
            ax.plot(t[MIN],np.cos(theta)[MIN],'o',color = 'blue',markersize =8)
            ax.plot(t[left_minima],np.cos(theta)[left_minima],'<',color = 'black',markersize = 8)
            ax.plot(t[right_minima],np.cos(theta)[right_minima],'>',color='black',markersize = 8)
            
        ################################################
        #### Plotting
        ################################################
        text = 'D = ' + str(np.round(D,5))
        ax.text(1.05, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        if k == len(axs)-Cols:
            ax.set_ylabel(r'$\cos(\theta)$', fontsize=30);
            ax.set_xlabel('time', fontsize=30)
            ax.xaxis.set_label_coords(0.5, -0.1);
            ax.yaxis.set_label_coords(-0.05, 0.5)

        set_scale(ax, [0,T], [-1,1])
        ax.set_xticklabels([0,T])
        ax.tick_params(labelsize=20)


    plt.savefig(save_path_name + 'pulses_alpha_'+str(alpha)+'.pdf', format='pdf')
    return(0)

#%%
def plot_pulses_sine(description_file,data_folder_ts,data_folder_pulses,save_path_name,dt,T,d,TS,N,delta):
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    plot_pulses_alpha_ = partial(plot_pulses_alpha_sine,data_folder_ts,data_folder_pulses,save_path_name,dt,T,d,TS,N,delta)
    pool.map(plot_pulses_alpha_,ref.groupby(['alpha']) )
    pool.close()
    pool.join()
    return (1)
# =============================================================================
# plotting pulses module
# =============================================================================

def plot_pulses_alpha_sine(data_folder_ts,data_folder_pulses,save_path_name,dt,T,d,TS,N,delta,tuple_):
    '''
    plor pulses for a certain alpha
    data_folder_ts : carpeta donde estan las series temporales
    data_folder_pulses : donde está la data de los pulsos
    save_path_name : donde queres guardar la figura  
    plotea el primer archivo solamente de la serie temporal! (variable j)
    '''
    i,rows = tuple_[0],tuple_[1]; j = 0
###############################################################################
### Parameters
###############################################################################

    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    PFE , PFI = get_fixed_points(alpha)
    print('alpha = ', alpha)
    T_n = ceil(int(T/dt)/d)

###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5,T+5] ; ylim = [-1.1,1.1] ;         
    Cols = TS; Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    colors =  sns.color_palette(sns.color_palette("viridis",Rows*1))
    colors =  colors[::1]

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.1)
    text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' #+' ~ ' + r'$ \Delta\theta = $'+ text_dtheta
    axs[0].text(0,1.5, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=35)
    
    
    for k,(ix,row_) in  enumerate(rows.groupby(['D'])):
        
        row = row_.iloc[j]; order = int(row.order); number = int(row.number)

        file_name =  str(number)+'_'+str(order)+'.pkl'
        file_name_max = 'max_xf_'+file_name

        theta = download_data(data_folder_ts + file_name)[:T_n:delta]
        t = time(dt,T,d)[::delta]
        print(row.D,'plotting ',order,number);D = row.D
        ax = axs[k]; ax.grid(False);
        ax.plot(t,np.sin(theta),linewidth=0.8,color = colors[k])
        
        
        if (check_file(file_name_max,data_folder_pulses)):
                        
            MAX          = mask_arr(T_n, download_data(data_folder_pulses + file_name_max))
            MIN          = mask_arr(T_n, download_data(data_folder_pulses + 'min_xf_'+ file_name))
            left_minima  = mask_arr(T_n, download_data(data_folder_pulses + 'left_minima_'+ file_name) )
            right_minima = mask_arr(T_n, download_data(data_folder_pulses + 'right_minima_'+ file_name) )
            
            ax.plot(t[MAX],np.sin(theta)[MAX],'o',color = 'red',markersize = 8)
            ax.plot(t[MIN],np.sin(theta)[MIN],'o',color = 'blue',markersize =8)
            ax.plot(t[left_minima],np.sin(theta)[left_minima],'<',color = 'black',markersize = 8)
            ax.plot(t[right_minima],np.sin(theta)[right_minima],'>',color='black',markersize = 8)
            
        ################################################
        #### Plotting
        ################################################
        text = 'D = ' + str(np.round(D,5))
        ax.text(1.05, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
        ax.set_ylim(ylim);
        ax.set_xlim(xlim)
        if k == len(axs)-Cols:
            ax.set_ylabel(r'$\cos(\theta)$', fontsize=30);
            ax.set_xlabel('time', fontsize=30)
            ax.xaxis.set_label_coords(0.5, -0.1);
            ax.yaxis.set_label_coords(-0.05, 0.5)

        set_scale(ax, [0,T], [-1,1])
        ax.set_xticklabels([0,T])
        ax.tick_params(labelsize=20)


    plt.savefig(save_path_name + 'pulses_sine_alpha_'+str(alpha)+'.pdf', format='pdf')
    return(0)

#%%
# =============================================================================
#     pulses quantifiers plotting block
# =============================================================================

def tune_plot(ax,x_label,y_label,xlim,xscale,ylim,yscale,axlabel_fs = 10,ticklabel_fs = 8):
    ''' Esto es una funcion auxiliar para plotear. La idea es escribir una vez
    el tamano de las labels, las escalas y esas cosas'''
    
    ax.set_xlabel(x_label, fontsize=axlabel_fs);
    ax.set_ylabel(y_label, fontsize=axlabel_fs);
    ax.set_xlim(xlim)
    ax.set_xticklabels([np.round(item/xscale,3) for item in ax.get_xticks()],fontsize=ticklabel_fs)
    ax.set_ylim(ylim)
    ax.set_yticklabels([np.round(item*yscale,3) for item in ax.get_yticks()],fontsize=ticklabel_fs)

def compute_st_values(ax,samples,bins,scale,fs = 6):
    ''' Esto esuna función auxiliar paraplotear.
    computa la moda de samples usando bins, y los cuartiles de samples usando samples. 
    samples es el dataset que queres estudiar, y bins es el histograma que estas graficando.
    scale es por loque tengo que dividir samples para llegar a minutos. fs es el fontsize del texto. '''           
    
    mode = (bins[1][np.argmax(bins[0])] + bins[1][np.argmax(bins[0])+1])/2 ; 
    mode = 'mode: '+str(np.round(mode/scale,2))+' min \n'
    ax.text(1, 0.85, mode, ha='right', va='center', transform=ax.transAxes, fontsize=fs) 
    ax.text(1, 0.9,r'$Q$: '+str(np.round(np.quantile(samples,0.25)/scale,2))+' ; '+str(np.round(np.quantile(samples,0.5)/scale,2))+' ; '
            +str(np.round(np.quantile(samples,0.75)/scale,2)) , ha='right', va='center', transform=ax.transAxes, fontsize=fs)
    ax.text(1, 0.8, 'total data: ' + str(len(samples)), ha='right', va='center', transform=ax.transAxes, fontsize=fs) 

#chequeamos que no haya numeros raros
    if True:
        neg = sum([1 for i in samples if i < bins[1][0]]) 
        pos =  sum([1 for i in samples if i > bins[1][-1]])
        text_borders = '< lower bound :'+str(neg) + '\n > upper bound :'+str(pos)
        print(text_borders)
        #ax.text(1, 0.75, text_borders, ha='right', va='center', transform=ax.transAxes, fontsize=fs) 

    
    
def plot_quantifiers_histograms(description_file,data_folder,save_path_name):
    '''
    data folder: donde están los pulsos
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha'])
    plot_quantifiers_histograms_D_ = partial(plot_quantifiers_histograms_D,data_folder,save_path_name)
    pool.map(plot_quantifiers_histograms_D_,tuple_)
    pool.close()
    pool.join()
    return (2)


def plot_quantifiers_histograms_D(data_folder,save_path_name,tuple_):
    plt.close('all')
    i,rows = tuple_[0],tuple_[1]
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    colors =  sns.color_palette(sns.color_palette("viridis",len(rows.groupby(['D'])) ))
    print(len(colors))
    colors =  colors[::1]
    for k,(D,row_) in  enumerate(rows.groupby(['D'])):
        dt = []; IPI = []; joint_duration = []; dm = []
        print(alpha,D,'plotting ');

        for (order,row) in row_.groupby(['order']):
        #para cada ensayo
            number = int(row.number)
            file_name =  str(number)+'_'+str(order)+'.pkl'

            if (check_file('dt_xf_'+file_name,data_folder)):        
            
                dt = dt + download_data(data_folder+'dt_xf_'+file_name)
                IPI = IPI + download_data(data_folder+'IPI_xf_'+file_name)
                dm = dm + download_data(data_folder+'dm_xf_'+file_name)
                joint_duration = joint_duration + download_data(data_folder+'joint_duration_xf_'+file_name)
                
            else:
                pass
                
           
        if len(dt) == 0:
            print('no quantifiers :( ')
        else:
            fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
            gs_row = gridspec.GridSpec(nrows=3, ncols=3, figure=fig, wspace=0.4, hspace=0.5)
            ax = plt.subplot(gs_row[0,0])            
            text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' +' ~ ' + r'$ D = $'+ str(D)
            ax.text(0.5,1.2, text, ha='center', va='center', transform=ax.transAxes, fontsize=10)
            
            #(0,15,76)
            #(0,30,151)
            #(0,45,226)
            #(0,90,451)
                
            if True:
                print(len(dt))
                bins = ax.hist(dt,bins=np.linspace(0,90,451)*1000, density=1,alpha=1,linewidth=1,color = colors[k]); 
                ax.plot(0.5*(bins[1][1:]+bins[1][:-1]),bins[0],alpha = 1,linewidth=1,color = colors[k] )                
                x_label = 'dt (min)'
                tune_plot(ax,x_label,r'probability density $(\times 1000)$',[0,30*1000],1000,[0,0.0003],1000)
                compute_st_values(ax,dt,bins,1000)   
            
                ax = plt.subplot(gs_row[0,1])            
                bins = ax.hist(IPI,bins=np.linspace(0,90,451)*1000, density=1, alpha=1,linewidth=1,color = colors[k]); 
                ax.plot(0.5*(bins[1][1:]+bins[1][:-1]),bins[0],alpha = 1,linewidth=1,color = colors[k] )                
                x_label = 'IPI (min)'
                tune_plot(ax,x_label,'',[0,40*1000],1000,[0,0.0003],1000)
                compute_st_values(ax,IPI,bins,1000)            
     
                ax = plt.subplot(gs_row[0,2])            
                bins = ax.hist(IPI,bins=np.linspace(0,90,451)*1000, density=1, alpha=1,linewidth=1,color = colors[k]); 
                ax.plot(0.5*(bins[1][1:]+bins[1][:-1]),bins[0],alpha = 1,linewidth=1,color = colors[k] )                
                ax.set_yscale('log')
                x_label = 'IPI (min)'
                tune_plot(ax,x_label,'',[0,40*1000],1000,[0,0.0003],1000)
                
                ax = plt.subplot(gs_row[1,1])      
                bins = ax.hist(dm,bins=np.linspace(0,90,451)*1000, density=1, alpha=1,linewidth=1,color = colors[k]); 
                ax.plot(0.5*(bins[1][1:]+bins[1][:-1]),bins[0],alpha = 1,linewidth=1,color = colors[k] )                
                x_label = 'silent intervals (min)'
                tune_plot(ax,x_label,'',[0,40*1000],1000,[0,0.0003],1000)
                compute_st_values(ax,dm,bins,1000) 
    
                ax = plt.subplot(gs_row[1,2])      
                bins = ax.hist(dm,bins=np.linspace(0,90,451)*1000, density=1, alpha=1,linewidth=1,color = colors[k]); 
                ax.plot(0.5*(bins[1][1:]+bins[1][:-1]),bins[0],alpha = 1,linewidth=1,color = colors[k] )   
                ax.set_yscale('log')
                x_label = 'silent intervals (min)'
                tune_plot(ax,x_label,'',[0,40*1000],1000,[0,0.0003],1000)
    
    
                ax = plt.subplot(gs_row[1,0])
                bins = ax.hist(joint_duration,bins=np.linspace(0,90,451)*1000, density=1, alpha=1,linewidth=1,color = colors[k]); 
                ax.plot(0.5*(bins[1][1:]+bins[1][:-1]),bins[0],alpha = 1,linewidth=1,color = colors[k] )                
                x_label = 'joint duration (min)'
                tune_plot(ax,x_label,r'probability density $(\times 1000)$',[0,30*1000],1000,[0,0.0003],1000)
                compute_st_values(ax,joint_duration,bins,1000)            
                
            else:
                bins = ax.hist(dt,bins=300, density=0, alpha=1,linewidth=1,color = colors[k]); 
                x_label = 'dt (min)'
                tune_plot(ax,x_label,r'counts $(\times 1000)$',[0,30*1000],1000,[0,5*1000],1000)
                compute_st_values(ax,dt,bins,1000)   
            
                ax = plt.subplot(gs_row[0,1])            
                bins = ax.hist(IPI,bins=300, density=0, alpha=1,linewidth=1,color = colors[k]); 
                x_label = 'IPI (min)'
                tune_plot(ax,x_label,r'counts $(\times 1000)$',[0,40*1000],1000,[0,5*1000],1000)
                compute_st_values(ax,IPI,bins,1000)            
 
                ax = plt.subplot(gs_row[1,0])      
                bins = ax.hist(dm,bins=300, density=0, alpha=1,linewidth=1,color = colors[k]); 
                x_label = 'silent intervals (min)'
                tune_plot(ax,x_label,r'counts $(\times 1000)$',[0,25*1000],1000,[0,20*1000],1000)
                compute_st_values(ax,dm,bins,1000) 

             
                ax = plt.subplot(gs_row[1,1])          
                bins = ax.hist(joint_duration,bins=300, density=0, alpha=1,linewidth=1,color = colors[k]); 
                x_label = 'joint duration (min)'
                tune_plot(ax,x_label,r'counts $(\times 1000)$',[0,90*1000],1000,[0,45*1000],1000)
                compute_st_values(ax,joint_duration,bins,1000)    
        
#                ax = plt.subplot(gs_row[2,0])            
#                ax.plot(dm,joint_duration,'o',alpha = 0.8, linewidth = 0)
#                ax.plot(np.arange(np.max(joint_duration)),np.arange(np.max(joint_duration))*2,linestyle='-',linewidth=0.5,color='black')
#                ax.set_xlim([0,300*1000]);ax.set_ylim([0,40*1000])
#                labels = [item/1000 for item in ax.get_xticks()]
#                ax.set_xticklabels(labels)
#                labels = [item/1000 for item in ax.get_yticks()]
#                ax.set_yticklabels(labels)
#                ax.set_xlabel('silence (min)', fontsize=8);
#                ax.set_ylabel('joint duration (min)', fontsize=8)
#            
        
            fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
            fig.savefig(save_path_name + 'histograms_alpha_'+str(alpha)+'D_'+str(D)+'.pdf', format='pdf')
            plt.close()
#%%

# =============================================================================
#     activity plot blocking
# =============================================================================     

def load_activity(row_,data_folder,dt,T,d,delta):
    '''
    load_activity(row_,dt,T,d,delta). Calcula la activity para cada par D, alpha.
    
    row_: DataFrame
        variable que viene del excel. Es el grupo de filas del excel de descripción,
        con un sólo D y alpha, pero con todos los trials. 
    '''    

    activity = []; n_cell = 0    
    for (order,row) in row_.groupby(['order']):
        
        number      = int(row.number)
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        
        if (check_file('dt_xf_'+file_name,data_folder)):        
            
            duration_cell   = download_data(data_folder+'dt_xf_'+file_name)     
            n_cell = n_cell + 1
            activity = activity + [sum(duration_cell) / len(time(dt,T,d)) *  100]
            
    activity = np.sort(activity)[::-1] #orden descendente
    silent = np.ones(len(activity)) * 100 - activity

    return activity,silent,n_cell


def plot_activity(dt,T,d,delta,description_file,data_folder,save_path_name):
    '''
    data folder: donde están los dt
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha'])
    plot_activity_D_ = partial(plot_activity_D,dt,T,d,delta,data_folder,save_path_name)
    pool.map(plot_activity_D_,tuple_)
    pool.close()
    pool.join()
    return (2)

def plot_activity_D(dt,T,d,delta,data_folder,save_path_name,tuple_):

    i,rows = tuple_[0],tuple_[1]
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    colors =  sns.color_palette(sns.color_palette("viridis",len(rows.groupby(['D'])) ))
    colors =  colors[::1]
    
    plt.close('all')
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=1, ncols=2, figure=fig,width_ratios= [0.75,0.25]); 
    gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.0,wspace= 0.5)
    gs_col0 = gridspec.GridSpecFromSubplotSpec(nrows=len(rows.groupby(['D'])), ncols=1, subplot_spec=gs_main[0], wspace=0.0, hspace=0.3)
    gs_col1 = gridspec.GridSpecFromSubplotSpec(nrows=len(rows.groupby(['D'])), ncols=1, subplot_spec=gs_main[1], wspace=0.0, hspace=0.3)
    
    population_activity,population_silent,D_labels = [],[],[]
 
    for k,(D,row_) in  enumerate(rows.groupby(['D'])):
        plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
        activity,silent,n_cell = load_activity(row_,data_folder,dt,T,d,delta)
        ax1 = plt.subplot(gs_col0[k])
        
        if len(activity) > 0:
            p1 = ax1.bar(np.arange(1 ,n_cell + 1),silent,width=1,color='darkgray',alpha=0.5,linewidth=0.0)
            p2 = ax1.bar(np.arange(1 ,n_cell + 1),activity,bottom=silent,width=1,color=colors[k],alpha=0.8,linewidth=0.0)
            population_activity.append(activity); population_silent.append(silent); D_labels.append(D)
            
        ax1.set_xlim([0,n_cell + +2 ]);ax1.set_ylim([0,100])
        ax1.set_xlabel( 'D : ' + str(D),fontsize=8); 
        ax1.set_xticks([1,n_cell + 1])
        ax1.set_yticks([0,50,100])
        ax1.tick_params(labelsize=6,direction='out', pad=1,length=2)
        ax1.xaxis.set_label_coords(0.5,-0.06)
        
        #mean activity
        plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
        ax2 = plt.subplot(gs_col1[k]);
        
        if len(activity) > 0:
            p1 = ax2.barh(np.arange(n_cell),width = np.mean(silent),xerr=np.std(silent),left =0,color='darkgray',alpha=0.5,linewidth=0.0,height=0.6)
            p2 = ax2.barh(np.arange(n_cell),width = np.mean(activity),left=np.mean(silent),xerr = np.std(activity),color=colors[k],alpha=0.8,linewidth=0.0,height=0.6)

        ax2.set_xticks([0,50,100])
        ax2.set_ylim([0,100])
        ax2.set_yticks([1,n_cell + 1])
        ax2.tick_params(labelsize=6,direction='out', pad=1,length=2)
        ax2.invert_yaxis()
    
    
    
    ax1.set_ylabel('fraction of cell track' ,fontsize=8); 
    ax2.set_xlabel('fraction of cell track' ,fontsize=8); 

    plt.savefig(save_path_name+'activity_'+str(alpha)+'.pdf', format='pdf')
    
    
    #boxplots
    plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); 
    ax_list = [plt.subplot(gs_main[0,0]),plt.subplot(gs_main[0,1])]
    ax_labels = ['activity', 'silent']
    #population activity/silent boxplot
    for n,j in enumerate([population_activity,population_silent]):
        ax = ax_list[n]
        X = [np.ones(len(j[i]))*(i+1) for i in range(0,len(j))]
        bp = ax.boxplot(j,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )
    
        for i,box_ in enumerate(bp['boxes']):
             box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
        for i,whisker in enumerate(bp['whiskers']):
            whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
        for i,cap in enumerate(bp['caps']):
            cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
        for i,median in enumerate(bp['medians']):
            median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
        for i,flyer in enumerate(bp['fliers']):
            flyer.set(markeredgecolor='black')## change color and linewidth of the medians
    
        for i in range(len(X)):
            xA = np.random.normal(0, 0.1,len(j[i])), 
            ax.scatter(xA+X[i],j[i], alpha=1,s = 0.2,color='black',edgecolors='black',linewidths=0.0)
           
        ax.tick_params(axis='x', labelsize=8,length=2); 
        ax.tick_params(axis='y', labelsize=8,length=2)
        ax.set_xticklabels([str(D) for D in D_labels],rotation = 0)
    
        ax.set_ylabel('fraction of cell track',fontsize=8)
        ax.set_xlabel('noise \n' + ax_labels[n],fontsize=8)
        ax.xaxis.set_label_coords(0.5, -0.05);ax.yaxis.set_label_coords(-0.05,0.5)
        ax.tick_params(labelsize=6,direction='out', pad=1,length=2)
        ax.set_ylim([0,100])
        ax.set_yticks([0,50,100])
        ax.tick_params(labelsize=6,direction='out', pad=1,length=2)


    plt.savefig(save_path_name+'bxp_activity_'+str(alpha)+'.pdf', format='pdf')

#%%

# =============================================================================
#     consecutiveness plot 
# =============================================================================     

def mean_consecutive_value(trials):
    if len(trials) > 0:
        arr_aux = []
        for j in range(np.max([len(i) for i in trials])):
            arr_aux.append([i[j] for i in trials if len(i) > j])
        return (np.array([np.mean(k) for k in arr_aux]),np.array([np.std(k) for k in arr_aux]))
    else:
        return ([],[])

def compute_total_consecutive(box_consecutiveness_cumulative):
    total_consecutive = []
    for i in box_consecutiveness_cumulative:
        if len(i) > 1:
            total_consecutive.append(i[1])
        else:
            total_consecutive.append(0)
    return total_consecutive

def load_consecutiveness_st(row_,data_folder):
    '''
    load_consecutiveness_st(row_,dt,T,d,delta). Calcula la estadistica de la consecutividad para cada par D, alpha.
    
    row_: DataFrame
        variable que viene del excel. Es el grupo de filas del excel de descripción,
        con un sólo D y alpha, pero con todos los trials. 
    '''    

    box_consecutiveness_cumulative , box_consecutiveness_non_cumulative= [],[]
    for (order,row) in row_.groupby(['order']):
        
        number      = int(row.number)
        file_name   =  str(number)+'_'+str(order)+'.pkl'
        if (check_file('c_nc_'+file_name,data_folder)):        
            box_consecutiveness_cumulative.append(download_data(data_folder+'c_c_'+file_name) )
            box_consecutiveness_non_cumulative.append(download_data(data_folder+'c_nc_'+file_name) )
    
    total,total_consecutive,total_isolated= [0],[0],[0]
    if len(box_consecutiveness_cumulative) > 0:
        total = [i[0] for i in box_consecutiveness_cumulative] 
        total_isolated = [i[0] for i in box_consecutiveness_non_cumulative] 
        total_consecutive = compute_total_consecutive(box_consecutiveness_cumulative)
        return mean_consecutive_value(box_consecutiveness_cumulative),mean_consecutive_value(box_consecutiveness_non_cumulative),(total, total_consecutive,total_isolated),([0,len(box_consecutiveness_cumulative)+1],[0.5,np.max(box_consecutiveness_cumulative)],[0.5,np.max(box_consecutiveness_non_cumulative)])
    else:
        return( mean_consecutive_value(box_consecutiveness_cumulative),mean_consecutive_value(box_consecutiveness_non_cumulative),(total, total_consecutive,total_isolated),([0,len(box_consecutiveness_cumulative)+1],[0.5,1],[0.5,1]))

def plot_consecutiveness(description_file,data_folder,save_folder):
    '''
    data folder: donde están los cc
    '''

    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha'])
    plot_consecutiveness_D_ = partial(plot_consecutiveness_D,data_folder,save_folder)
    pool.map(plot_consecutiveness_D_,tuple_)
    pool.close()
    pool.join()
    return (2)

def plot_consecutiveness_D(data_folder,save_folder,tuple_):
    i,rows = tuple_[0],tuple_[1]
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    colors =  sns.color_palette(sns.color_palette("viridis",len(rows.groupby(['D'])) ))
    colors =  colors[::1]
    plt.close('all')
    plt.rc('axes.spines', top=False, bottom=True, left=True, right=False); 
    fig = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    #gs_main = gridspec.GridSpec(nrows=2, ncols=2, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)
    gs_main = gridspec.GridSpec(nrows=2, ncols=1, figure=fig); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

    ax1 = plt.subplot(gs_main[0])
    ax2 = plt.subplot(gs_main[1])
    #ax3 = plt.subplot(gs_main[1,0])
    #ax4 = plt.subplot(gs_main[1,1])
    box_total,box_consecutive,box_isolated,D_labels = [],[],[],[]
    for k,(D,row_) in  enumerate(rows.groupby(['D'])):
        (mean_c_c,std_c_c),(mean_c_nc,std_c_nc),(total,total_consecutive,total_isolated),(X_lim,YC_lim,YNC_lim) = load_consecutiveness_st(row_,data_folder)
        box_total.append(total);box_consecutive.append(total_consecutive);box_isolated.append(total_isolated)

        if len(mean_c_c) > 0:
            ax1.plot(np.arange(1,len(mean_c_c)+1),mean_c_c, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[k],label=str(D))
            ax1.fill_between(np.arange(1,len(mean_c_c)+1),mean_c_c-std_c_c,mean_c_c+std_c_c,color =colors[k],alpha = 0.2)
            X_lim = [0,50]
            ax1.set_xlim(X_lim);
            ax1.set_yscale('log')
            #ax1.set_ylim(YC_lim)
            ax1.set_ylabel('counts',fontsize=10); ax1.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
            ax1.xaxis.set_label_coords(0.5, -0.08);ax1.yaxis.set_label_coords(-0.2,0.5);
            #ax1.set_xticks([0,3,6,9,12,15])
            
            ax2.plot(np.arange(1,len(mean_c_nc)+1),mean_c_nc, color=colors[k], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=str(D))
            ax2.fill_between(np.arange(1,len(mean_c_nc)+1),mean_c_nc-std_c_nc,mean_c_nc+std_c_nc,color =colors[k],alpha = 0.2)        
            ax2.set_yscale('log');    X_lim = [0,50]
            ax2.set_xlim(X_lim); #ax2.set_ylim(YNC_lim)
            ax2.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses ',fontsize=10)
            ax2.xaxis.set_label_coords(0.5, -0.08);
            ax2.yaxis.set_label_coords(-0.2,0.5);
            #ax2.set_xticks([0,3,6,9,12,15])
            
#            mean_c_c = np.array([i/mean_c_c[0] for i in mean_c_c])
#            std_c_c = np.array([i/mean_c_c[0] for i in std_c_c])
#            ax3.plot(np.arange(1,len(mean_c_c)+1),mean_c_c, linewidth=0.5, marker = "." , markersize=7, alpha=1,color =colors[k],label=str(D))
#            ax3.fill_between(np.arange(1,len(mean_c_c)+1),mean_c_c-std_c_c,mean_c_c+std_c_c,color =colors[k],alpha = 0.2)
#            ax3.set_xlim(X_lim);
#            ax3.set_yscale('log')
#            #ax3.set_ylim([0.5,10])
#            ax3.set_ylabel('counts',fontsize=10); ax3.set_xlabel('length of sequence of \n consecutive pulses',fontsize=10)
#            ax3.xaxis.set_label_coords(0.5, -0.08);ax3.yaxis.set_label_coords(-0.2,0.5);
#            #ax3.set_xticks([0,3,6,9,12,15])
#
#            mean_c_nc = np.array([i/mean_c_nc[0] for i in mean_c_nc])
#            std_c_nc = np.array([i/mean_c_nc[0] for i in std_c_nc])
#            ax4.plot(np.arange(1,len(mean_c_nc)+1),mean_c_nc, color=colors[k], linewidth=0.5, marker = "." , markersize=7, alpha=1,label=str(D))
#            ax4.fill_between(np.arange(1,len(mean_c_nc)+1),mean_c_nc-std_c_nc,mean_c_nc+std_c_nc,color =colors[k],alpha = 0.2)        
#            ax4.set_yscale('log')
#            ax4.set_xlim(X_lim); #ax4.set_ylim([0.5,10])
#            ax4.set_ylabel('counts',fontsize=10); ax2.set_xlabel('length of sequence of \n consecutive pulses ',fontsize=10)
#            ax4.xaxis.set_label_coords(0.5, -0.08);
#            ax4.yaxis.set_label_coords(-0.2,0.5);
#            #ax4.set_xticks([0,3,6,9,12,15])


    ax1.legend(fontsize=6, ncol=1, framealpha=0, fancybox=True)
    plt.savefig(save_folder+'consecutiveness_'+str(alpha)+'.pdf', format='pdf')
    
    figB = plt.figure(constrained_layout=False, figsize=(8.27, 11.692))
    gs_mainB = gridspec.GridSpec(nrows=3, ncols=1, figure=figB); gs_main.update(left=0.1, right=0.9, bottom=0.1, top=0.90, hspace=0.3,wspace=0.3)

    ax1B = plt.subplot(gs_mainB[0,0])
    ax2B = plt.subplot(gs_mainB[1,0])
    ax3B = plt.subplot(gs_mainB[2,0])
    X1 = [np.ones(len(box_total[i]))*(i+1) for i in range(0,len(box_total))]
    bp1 = ax1B.boxplot(box_total,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

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
        xA = np.random.normal(0, 0.1, len(box_total[i])), 
        ax1B.scatter(xA+X1[i],box_total[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)

    ax1B.tick_params(axis='x', labelsize=8,length=2); 
    ax1B.tick_params(axis='y', labelsize=8,length=2)
    ax1B.set_xticklabels(D_labels,rotation = 0)
    ax1B.set_xlabel('total pulses',fontsize=8)
    ax1B.set_ylabel('counts',fontsize=8)
    ax1B.set_ylim([-1,200])
    ax1B.xaxis.set_label_coords(0.5, -0.12);ax1B.yaxis.set_label_coords(-0.05,0.5)
    ax1B.tick_params(labelsize=6,direction='out', pad=1,length=2)

    X2 = [np.ones(len(box_isolated[i]))*(i+1) for i in range(0,len(box_isolated))]
    bp2 = ax2B.boxplot(box_isolated,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp2['boxes']):
         box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp2['whiskers']):
        whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp2['caps']):
        cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp2['medians']):
        median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp2['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians

    for i in range(len(X2)):
        xA = np.random.normal(0, 0.1, len(box_isolated[i])), 
        ax2B.scatter(xA+X2[i],box_isolated[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
           
        
    ax2B.tick_params(axis='x', labelsize=8,length=2); 
    ax2B.tick_params(axis='y', labelsize=8,length=2)
    ax2B.set_xticklabels(D_labels,rotation = 0)
    ax2B.set_xlabel('isolated pulses',fontsize=8)
    ax2B.set_ylabel('counts',fontsize=8)
    ax2B.xaxis.set_label_coords(0.5, -0.12);ax2B.yaxis.set_label_coords(-0.05,0.5)
    ax2B.set_ylim([-1,200])
    ax2B.tick_params(labelsize=6,direction='out', pad=1,length=2)
    

    X3 = [np.ones(len(box_consecutive[i]))*(i+1) for i in range(0,len(box_consecutive))]
    bp3 = ax3B.boxplot(box_consecutive,vert=True,whis=[5, 95],patch_artist=True,showmeans=False,meanline=True,showfliers=False )

    for i,box_ in enumerate(bp3['boxes']):
         box_.set( color=colors[i], linewidth=0.0,facecolor=colors[i],alpha = 0.1)# change outline color
    for i,whisker in enumerate(bp3['whiskers']):
        whisker.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)
    for i,cap in enumerate(bp3['caps']):
        cap.set(color=colors[i//2],linestyle = '-', linewidth=1,alpha=0.3)## change color and linewidth of the caps
    for i,median in enumerate(bp3['medians']):
        median.set(color=colors[i],linestyle = '-', linewidth=1.5)## change color and linewidth of the medians
    for i,flyer in enumerate(bp3['fliers']):
        flyer.set(markeredgecolor='black')## change color and linewidth of the medians

    for i in range(len(X3)):
        xA = np.random.normal(0, 0.1, len(box_consecutive[i])), 
        ax3B.scatter(xA+X2[i],box_consecutive[i], alpha=1,s = 1.5,color='black',edgecolors='black',linewidths=0.0)
           
        
    ax3B.tick_params(axis='x', labelsize=8,length=2); 
    ax3B.tick_params(axis='y', labelsize=8,length=2)
    ax3B.set_xticklabels(D_labels,rotation = 0)
    ax3B.set_xlabel('consecutive pulses',fontsize=8)
    ax3B.set_ylabel('counts',fontsize=8)
    ax3B.xaxis.set_label_coords(0.5, -0.12);ax2.yaxis.set_label_coords(-0.05,0.5)
    ax3B.set_ylim([-1,200])
    ax3B.tick_params(labelsize=6,direction='out', pad=1,length=2)
    

    plt.savefig(save_folder+ 'consecutiveness_bx_'+str(alpha)+'.pdf', format='pdf')
