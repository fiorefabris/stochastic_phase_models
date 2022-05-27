import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import ceil
from functools import partial 
import seaborn as sns

from adler.data_managing_functions import download_data,check_file



#%%    
    

###############################################################################
### Function for paralelizing the FPT plotting 
###############################################################################
    
def plot_FPT(description_file,save_path_name,dt,T,d,data_folder):
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    plot_FPT_alpha_ = partial(plot_FPT_alpha,save_path_name,dt,T,d,data_folder)
    pool.map(plot_FPT_alpha_,ref.groupby(['alpha']) )
    pool.close()
    pool.join()
    return (1)


###############################################################################
### Function for plotting the FPT distribution 
###############################################################################


def plot_FPT_alpha(save_path_name,dt,T,d,data_folder,tuple_):
    '''
    FPT plotting function
    '''
    i,rows = tuple_[0],tuple_[1]


###############################################################################
### Parameters
###############################################################################
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    

###############################################################################
### Plotting parameters
###############################################################################    
    plt.rcdefaults();       
    Cols = 1; Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    ylim = [0,0.15] ; #xlim = [0,5] 
    xlim = [0,100] ; 

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel();
    text = r'$\omega = \frac{2\pi}{7 min}$' + ' ~ ' + r'$ \alpha = $'+str(alpha) + r'$ \frac{2\pi}{7 min}$' 
    axs[0].text(0.2,1.3, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=12)

    
    for k,(D,row) in enumerate(rows.groupby(['D'])):
        ax = axs[k]; ax.grid(False);
        
        if check_file('fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl',save_path_name):
            print(alpha,np.round(D,5),'available')        
            
            FPT = download_data(save_path_name+'FPT_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
            bins = ax.hist(FPT,bins=100,range=(0,200) , density=1, alpha=0.8,linewidth=0,color = 'green'); 
        
            if len(FPT) > 0:
                    mode_FPT = (bins[1][np.argmax(bins[0])] + bins[1][np.argmax(bins[0])+1])/2 ; mode_FPT = 'mode: '+str(np.round(mode_FPT,2))+' min \n'
                    ax.text(1, 0.50, mode_FPT, ha='right', va='center', transform=ax.transAxes, fontsize=10) 
                    ax.text(1, 0.8,r'$Q$: '+str(np.round(np.quantile(FPT,0.25),2))+' ; '+str(np.round(np.quantile(FPT,0.5),2))+' ; '
                                  +str(np.round(np.quantile(FPT,0.75),2)) , ha='right', va='center', transform=ax.transAxes, fontsize=10)
                    ax.text(1, 0.3, 'total data: ' + str(len(FPT)), ha='right', va='center', transform=ax.transAxes, fontsize=10) 
        else:
            print(alpha,D,'empty')
            
        ax.set_xlim(xlim);
        ax.set_ylim(ylim);
        #ax.set_xticks(xticks);
        ax.set_yticks(ylim)
        ax.tick_params(labelsize=10)
        
        text = 'D = ' + str(np.round(D,5))
        ax.text(1.10, 0.7, text , ha='center', va='center', transform=ax.transAxes, fontsize=12)

#    xticks = np.arange(0,200,14); ax.set_xticks(xticks)
#    ax.set_xticklabels(xticks, fontsize=15);
    
    ax.set_ylabel('density', fontsize=15);
    ax.set_xlabel('first passage time (mins)', fontsize=15)
    ax.xaxis.set_label_coords(0.5, -0.35);
    ax.yaxis.set_label_coords(-0.10, 1)
    
    #for m in range(Cols*Rows - k):
        #fig.delaxes(axs[-m-1])
    
    
    plt.savefig(save_path_name + 'FPT_alpha_'+str(alpha)+'.pdf', format='pdf')

    return(0)
#%%

def plot_FPT_all(description_file,save_path_name,dt,T,d,data_folder):
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    pool = mp.Pool(processes= ceil(mp.cpu_count()/4))
    plot_FPT_alpha_ = partial(plot_FPT_alpha_all,save_path_name,dt,T,d,data_folder)
    pool.map(plot_FPT_alpha_,ref.groupby(['alpha']) )
    pool.close()
    pool.join()
    return (1)

def plot_FPT_alpha_all(save_path_name,dt,T,d,data_folder,tuple_):
    '''
    FPT plotting function
    '''
    i,rows = tuple_[0],tuple_[1]


###############################################################################
### Parameters
###############################################################################
    omega =  rows.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    

###############################################################################
### Plotting parameters
###############################################################################    
    plt.rcdefaults();       
    Cols = 1; Tot = len(rows.groupby(['D'])) ;
    Rows = ceil(Tot/ Cols)
    colors =  sns.color_palette(sns.color_palette("viridis",Rows*1))
    colors =  colors[::1]
    ylim = [0,0.3] ; #xlim = [0,5] 
    xlim = [0,100] ; 

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel();
    text = r'$\omega = \frac{2\pi}{7 min}$' + ' ~ ' + r'$ \alpha = $'+str(alpha) + r'$ \frac{2\pi}{7 min}$' 
    axs[0].text(0.2,1.05, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=12)

    
    for k,(D,row) in enumerate(rows.groupby(['D'])):
        ax = axs[0]; ax.grid(False);
        
        assert check_file('fft_yf_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl',data_folder)
            
        FPT = download_data(data_folder+'FPT_'+str(omega)+'_'+str(alpha)+'_'+str(D)+'.pkl')
        n,bins,_ = ax.hist(FPT,bins=100,range=(0,200) , density=1, alpha=0.0,linewidth=0,color = colors[k]); 
        ax.plot(0.5*(bins[1:]+bins[:-1]),n,alpha = 1,linewidth=1,color = colors[k],label = str(D) )
            
        ax.set_xlim(xlim);
        ax.set_ylim(ylim);
        #ax.set_xticks(xticks);
        ax.set_yticks(ylim)
        ax.tick_params(labelsize=10)
        
        #text = 'D = ' + str(np.round(D,5))
        #ax.text(1.10, 0.7, text , ha='center', va='center', transform=ax.transAxes, fontsize=12)

#    xticks = np.arange(0,200,14); ax.set_xticks(xticks)
#    ax.set_xticklabels(xticks, fontsize=15);
    
    ax.set_ylabel('density', fontsize=15);
    ax.set_xlabel('first passage time (mins)', fontsize=15)
    ax.xaxis.set_label_coords(0.5, -0.1);
    ax.yaxis.set_label_coords(-0.10,0.5)
    ax.legend(fontsize=7, ncol=1,framealpha=0, fancybox=True)

    
    #for m in range(Cols*Rows - k):
        #fig.delaxes(axs[-m-1])
    
    
    plt.savefig(save_path_name + 'all_FPT_alpha_'+str(alpha)+'.pdf', format='pdf')

    return(0)

#%%
def plot_FPT_square(dt,T,d,description_file,data_folder,save_path_name):
    '''
    FPT plotting function
    '''
######## Getting data information
    plt.rcdefaults();
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    Cols = len(ref.groupby(['D'])) ;
    Rows = len(ref.groupby(['alpha'])) ; 
    colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
    colors =  colors[::1]

    ylim = [0,0.15] ; #xlim = [0,5] 
    xlim = [0,100] ; 

###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.3, hspace=0.3)

    
    for col,(D,col_) in  enumerate(ref.groupby(['D'])):
        if 'alpha' in ref.keys() : iterator = col_.groupby(['alpha'])
        if 'delta' in ref.keys() : iterator = col_.groupby(['delta'])
        for row, (alpha,row_)  in  enumerate(iterator):
            if 'omega' in row.keys() : omega = col_.omega.unique()[0]
            if 'T0' in row.keys() : omega = col_.T0.unique()[0]

            if 'alpha' in row.keys() :delta= np.round(alpha/omega,4)  
            if 'delta' in row.keys() :delta= alpha


            file_name = 'FPT_'+str(omega)+'_'+str(delta)+'_'+str(D)+'.pkl'
            ax = axs[row,col]; ax.grid(False);
            
            ################################################
            #### download data
            ################################################
            if check_file(file_name,data_folder):        
            
                FPT = download_data(data_folder+file_name)
                bins = ax.hist(FPT,bins=100,range=(0,200) , density=1, alpha=0.8,linewidth=0,color=colors[col]); 
        
                if len(FPT) > 0:
                        mode_FPT = (bins[1][np.argmax(bins[0])] + bins[1][np.argmax(bins[0])+1])/2 ; mode_FPT = 'mode: '+str(np.round(mode_FPT,2))+' min \n'
                        ax.text(1, 0.75, mode_FPT, ha='right', va='center', transform=ax.transAxes, fontsize=7) 
                        ax.text(1, 0.9,r'$Q$: '+str(np.round(np.quantile(FPT,0.25),2))+' ; '+str(np.round(np.quantile(FPT,0.5),2))+' ; '
                                      +str(np.round(np.quantile(FPT,0.75),2)) , ha='right', va='center', transform=ax.transAxes, fontsize=7)
                        ax.text(1, 0.7, 'total data: ' + str(len(FPT)), ha='right', va='center', transform=ax.transAxes, fontsize=7) 
                else:
                    print(delta,D,'empty')
            
            ax.set_xlim(xlim);
            ax.set_ylim(ylim);
            #ax.set_xticks(xticks);
            ax.set_yticks(ylim)
            ax.tick_params(labelsize=7)
            
            if row == 0:
                text = 'D = ' + str(np.round(D,5))
                ax.text(0.8, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=7)
            if col == 0:
                text = 'delta = ' + str(delta)
                ax.text(-0.25, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=7)

#    xticks = np.arange(0,200,14); ax.set_xticks(xticks)
#    ax.set_xticklabels(xticks, fontsize=15);
            if (row == Rows-1) and (col == 0): 
                ax.set_ylabel('density', fontsize=7)
                ax.set_xlabel('first passage time (mins)', fontsize=7)
                ax.xaxis.set_label_coords(0.5,  -0.2);
                ax.yaxis.set_label_coords(-0.05, 0.5)

    
    #for m in range(Cols*Rows - k):
        #fig.delaxes(axs[-m-1])
    
    
    plt.savefig(save_path_name + 'FPT_square.pdf', format='pdf')

    return(0)