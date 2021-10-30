import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.plotting.plotting_main import download_data,check_file


def plot_dt_alpha(description_file,data_folder,save_path_name):
    
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);

    colors =  sns.color_palette(sns.color_palette("viridis",len(ref.groupby(['D']))))
    colors =  colors[::1]

    for k,(D,row) in enumerate(ref.groupby(['D'])):

        
###############################################################################
### Parameters
###############################################################################

        omega =  row.omega.unique()[0]
        alphas = row.alpha.values/omega
    

###############################################################################
### Figure
###############################################################################    

        fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(8.27, 11.69))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
        axs = axs.ravel();        
        ax = axs[0]; ax.grid(False);

###############################################################################
### Download data
###############################################################################    
        dt = []; 
        for (number,row_) in row.groupby(['number']):
            order = row_.order
            file_name =  str(number)+'_'+str(order)+'.pkl'
            print(file_name,data_folder)
            if (check_file('dt_xf_'+file_name,data_folder)):        
                dt.append(download_data(data_folder+'dt_xf_'+file_name))
            else:
                dt.append([0])        
################################################
#### Plotting
################################################
        print(dt)
        ax.plot(alphas,[np.mean(i) for i in dt],linewidth=2,color=colors[k],label = D)

#        ax.set_ylim(ylim);
#        ax.set_xlim(xlim)
    ax.set_ylabel('mean duration', fontsize=30);
    ax.set_xlabel('alpha/omega', fontsize=30)
    ax.xaxis.set_label_coords(0.5, -0.1);
    ax.yaxis.set_label_coords(-0.05, 0.5)

    ax.set_xticklabels(alphas)
    ax.tick_params(labelsize=20)

    plt.savefig(save_path_name + 'dt_alpha.pdf', format='pdf')
    return(0)

   