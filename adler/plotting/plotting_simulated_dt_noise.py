#import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adler.data_managing_functions import download_data#check_file
from adler.get_time_series.main import all_combinations
from adler.plotting.plotting_main import set_scale,silent_ax

#%%
    
###############################################################################
###  plotting
###############################################################################  
    
def plot_theta_alpha(data_folder,save_path_name,params):

    
    #colors =  sns.color_palette(sns.color_palette("viridis",len(params['alpha'])))
    
    D_ = len(params['D']); ALP = len(params['alpha'])
    
    fig, axs = plt.subplots(D_,ALP, sharex=False, sharey=True, figsize=(8.27, 11.69))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8, wspace=0.1, hspace=0.2)
    axs = axs.ravel(); 
  
    #para cada fila
    for k,(omega,alpha,D) in enumerate(all_combinations(params)):
            
        cond_prob_filename = data_folder +  'cond_prob_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        initial_conditions_filename = data_folder +  'initial_conditions_omega_'+str(np.round(omega,3))+'_alpha_'+str(np.round(alpha/omega,3))+'_D_'+str(D)+'.pkl'
        cond_prob = download_data(cond_prob_filename)
        initial_conditions = download_data(initial_conditions_filename)

        ax = axs[k]

        ax.plot(initial_conditions,cond_prob,linewidth=1) #,color=colors[k]
        #ax.set_ylim([-1.05,1.05]); ax.set_xlim([-1.05,1.05])
              
        if k == (D_*ALP - D_  ):
            ax.set_ylabel("frequency over 100000", fontsize=10);
            ax.set_xlabel("initial conditions", fontsize=10)
            ax.xaxis.set_label_coords(0.5, -0.1);
            ax.yaxis.set_label_coords(-0.05, 0.5)
           # set_scale(ax,[-1,1],[-1,1])
           # ax.set_xticklabels([-1,1]); ax.set_yticklabels([-1,1]); ax.tick_params(labelsize=10) 

        if k < D_:
            text = r'$\alpha = $' + str(alpha) 
            ax.text(0.7, 0.8, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)
            silent_ax(ax)


        else:
            silent_ax(ax)

        if ALP % k == 0 :
            text = 'D = ' + str(D)
            ax.text(0.7, 1.1, text , ha='center', va='center', transform=ax.transAxes, fontsize=5)

        
    plt.savefig(save_path_name + 'plotting.pdf', format='pdf')
    return(0)


