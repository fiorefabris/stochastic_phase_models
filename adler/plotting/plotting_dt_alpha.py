import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.plotting.plotting_main import download_data,check_file

#%%
def download_dt_alpha(row,data_folder):
    '''
    te arma una lista con todos los dt, en donde cada lugar de la lista (que es una lista) es un alpha distinto
    si un dt no esta en la carpeta, devuelve un cero
    '''
    dt = []; 
    for (number,row_) in row.groupby(['number']):
        order = 0
        print('order',row_.order.values)
        file_name =  str(number)+'_'+str(order)+'.pkl'
        if (check_file('dt_xf_'+file_name,data_folder)):        
            dt.append(download_data(data_folder+'dt_xf_'+file_name))
        else:
            dt.append([0])        
    return(dt)



#%%
#tiene ax[1]
# el numero de orden order donde levantas el archivo es cero. 

def plot_dt_alpha(description_file,data_folder,save_path_name):
    
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);

    colors =  sns.color_palette(sns.color_palette("viridis",len(ref.groupby(['D']))))

###############################################################################
### Figure
###############################################################################  
    
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel();        
    

    for k,(D,row) in enumerate(ref.groupby(['D'])):
        
###############################################################################
### Parameters and data -- plotting
###############################################################################

        omega =  row.omega.unique()[0]
        alphas = row.alpha.values/omega
        
        dt =download_dt_alpha(row,data_folder)
        axs[0].plot(alphas,[np.mean(i) for i in dt],'o',color=colors[k],label = D)
        axs[1].plot(alphas,[np.mean(i) for i in dt],'o',color=colors[k],label = D)



#        ax.set_ylim(ylim);
#        ax.set_xlim(xlim)
        
    axs[1].set_yscale('log'); axs[0].set_yscale('linear')
    axs[1].set_ylabel('mean duration', fontsize=10);
    axs[1].set_xlabel('alpha/omega', fontsize=10)
    
    axs[1].xaxis.set_label_coords(0.5, -0.1);
    axs[1].yaxis.set_label_coords(-0.05, 0.5)
    
    axs[0].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    axs[1].legend(fontsize=8, ncol=1, framealpha=0, fancybox=True)
    
    axs[1].set_xticklabels([np.round(i,2) for i in alphas],fontsizte = 10)
    axs[1].tick_params(labelsize=10)

    plt.savefig(save_path_name + 'dt_alpha.pdf', format='pdf')
    return(0)

