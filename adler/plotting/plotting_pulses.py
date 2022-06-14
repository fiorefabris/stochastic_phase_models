import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from adler.data_managing_functions import download_data,check_file,time
from adler.plotting.plotting_main import set_scale,create_df,download_quantifiers,mask_arr




#%%
# =============================================================================
# plotting pulses module
# =============================================================================


def plot_pulses_square(dt,beg,T,d,N,Delta,description_file,data_folder,data_folder_pulses,save_path_name):
    '''
    
    '''

######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.05,2.05] ;         
    Cols = len(ref.groupby(['D'])) ;
    Rows = len(ref.groupby(['alpha'])) ; 
    colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    


    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    #text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' 
    #axs[0].text(0,1.5, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=35)
    
    for col,(D,col_) in  enumerate(ref.groupby(['D'])):
        for row, (alpha,row_)  in  enumerate(col_.groupby(['alpha'])):
            delta= np.round(alpha/col_.omega.unique()[0],4)  
        
            order = int(row_.order); number = int(row_.number)
            file_name =  str(number)+'_'+str(order)+'.pkl'
            ax = axs[row,col]; ax.grid(False);
            
        
            ################################################
            #### download data
            ################################################
            if check_file(file_name,data_folder):            
                
                theta = download_data(data_folder + file_name) 
                t = time(dt,T+beg,d)
                end = len(t)
                beg_ = int(beg/(dt*d))
                ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color=colors[col])
            
            if check_file('max_'+file_name,data_folder_pulses):            
                
                
                MAX          = mask_arr(beg_,end, download_data(data_folder_pulses + 'max_'+file_name))
                left_minima  = mask_arr(beg_,end, download_data(data_folder_pulses + 'left_minima_'+ file_name) )
                right_minima = mask_arr(beg_,end, download_data(data_folder_pulses + 'right_minima_'+ file_name) )
                
                if len(MAX) > 0:
                    ax.plot(t[beg_:end][MAX],(1+ np.sin(theta))[beg_:end][MAX],'o',color = 'blue',markersize = 8)
                    ax.plot(t[beg_:end][left_minima],(1+ np.sin(theta))[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                    ax.plot(t[beg_:end][right_minima],(1+ np.sin(theta))[beg_:end][right_minima],'>',color='black',markersize = 8)


            
            ###############################################
            #### Plotting
            ################################################
            if row == 0:
                text = 'D = ' + str(np.round(D,5))
                ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if col == 0:
                text = 'delta = ' + str(delta)
                ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)

            ax.set_ylim(ylim);
            ax.set_xlim(xlim)
            
            if (row == Rows-1) and (col == 0): 
                ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                ax.set_xlabel('time', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.1);
                ax.yaxis.set_label_coords(-0.05, 0.5)
            
            set_scale(ax,[beg,T], [0,2])
            ax.set_xticklabels([beg,T])
            ax.set_yticklabels([0,2])
            ax.tick_params(labelsize=20)


#    for m in range((Cols*Rows - (k+1))):
#        fig.delaxes(axs[-m-1])


    plt.savefig(save_path_name + 'pulses_square.pdf', format='pdf')
    return(0)

#%%
def plot_pulses_square_ou(dt,beg,T,d,N,Delta,description_file,data_folder,data_folder_pulses,save_path_name):
    '''
    
    '''

######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    xlim = [-5+beg,T+5] ; ylim = [-0.05,2.05] ;         
    for alpha,ref in ref_.groupby(['alpha0']):
        
        Rows = len(ref.groupby(['tau'])) ; 
        Cols = len(ref.groupby(['sigma'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    


        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
        #text = r'$\omega = \frac{2\pi}{7 min}$' +' ~ ' + r'$\alpha = $' +str(alpha) + r'$ \frac{2\pi}{7 min}$' 
        #axs[0].text(0,1.5, text, ha='center', va='center', transform=axs[0].transAxes, fontsize=35)
        
        for col,(sigma,col_) in  enumerate(ref.groupby(['sigma'])):
            for row, (tau,row_)  in  enumerate(col_.groupby(['tau'])):
                delta= np.round(alpha/col_.omega.unique()[0],4)  
            
                order = int(row_.order.values[0]); number = int(row_.number.values[0])
                file_name =  str(number)+'_'+str(order)+'.pkl'
                ax = axs[row,col]; ax.grid(False);
                
            
                ################################################
                #### download data
                ################################################
                if check_file(file_name,data_folder):            
                    
                    theta = download_data(data_folder + file_name) 
                    t = time(dt,T+beg,d)
                    end = len(t)
                    beg_ = int(beg/(dt*d))
                    ax.plot(t[beg_:end:Delta],1+np.sin(theta)[beg_:end:Delta],linewidth=2,color=colors[col])
                
                if check_file('max_'+file_name,data_folder_pulses):            
                    
                    
                    MAX          = mask_arr(beg_,end, download_data(data_folder_pulses + 'max_'+file_name))
                    left_minima  = mask_arr(beg_,end, download_data(data_folder_pulses + 'left_minima_'+ file_name) )
                    right_minima = mask_arr(beg_,end, download_data(data_folder_pulses + 'right_minima_'+ file_name) )
                    
                    if len(MAX) > 0:
                        ax.plot(t[beg_:end][MAX],(1+ np.sin(theta))[beg_:end][MAX],'o',color = 'blue',markersize = 8)
                        if left_minima.size >0: ax.plot(t[beg_:end][left_minima],(1+ np.sin(theta))[beg_:end][left_minima],'<',color = 'black',markersize = 8)
                        if right_minima.size >0: ax.plot(t[beg_:end][right_minima],(1+ np.sin(theta))[beg_:end][right_minima],'>',color='black',markersize = 8)
    
    
                
                ###############################################
                #### Plotting
                ################################################
                if row == 0:
                    text = 'sigma = ' + str(np.round(sigma,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'tau = ' + str(np.round(tau,5))
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
    
                ax.set_ylim(ylim);
                ax.set_xlim(xlim)
                
                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel(r'$1 + \sin(\theta)$', fontsize=30);
                    ax.set_xlabel('time', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.1);
                    ax.yaxis.set_label_coords(-0.05, 0.5)
                
                set_scale(ax,[beg,T], [0,2])
                ax.set_xticklabels([beg,T])
                ax.set_yticklabels([0,2])
                ax.tick_params(labelsize=20)
    

    
        plt.savefig(save_path_name + 'pulses_square_ou_'+str(delta)+'.pdf', format='pdf')
    return(0)

#%%
# =============================================================================
#     pulses quantifiers plotting block
# =============================================================================

def tune_plot(ax,x_label,y_label,xlim,xscale,ylim,yscale,axlabel_fs = 10,ticklabel_fs = 8):
    ''' Esto es una funcion auxiliar para plotear. La idea es escribir una vez
    el tamano de las labels, las escalas y esas cosas'''
    
    #ax.set_xlabel(x_label, fontsize=axlabel_fs);
    #ax.set_ylabel(y_label, fontsize=axlabel_fs);
    ax.set_xlim(xlim)
    ax.set_xticklabels([np.round(item/xscale,3) for item in ax.get_xticks()],fontsize=ticklabel_fs)
    ax.set_ylim(ylim)
    ax.set_yticklabels([np.round(item*yscale,3) for item in ax.get_yticks()],fontsize=ticklabel_fs)

def compute_st_values(ax,samples,bins,scale,fs = 6):
    ''' Esto esuna función auxiliar paraplotear.
    computa la moda de samples usando bins, y los cuartiles de samples usando samples. 
    samples es el dataset que queres estudiar, y bins es el histograma que estas graficando.
    scale es por loque tengo que dividir samples para llegar a minutos. fs es el fontsize del texto. '''           
    
    if len(samples) <= 0:
        pass
    else:
        mode = (bins[1][np.argmax(bins[0])] + bins[1][np.argmax(bins[0])+1])/2 ; 
        mode = 'mode: '+str(np.round(mode/scale,2))+' min \n'
        ax.text(1, 0.75, mode, ha='right', va='center', transform=ax.transAxes, fontsize=fs) 
        ax.text(1, 0.9,r'$Q$: '+str(np.round(np.quantile(samples,0.25)/scale,2))+' ; '+str(np.round(np.quantile(samples,0.5)/scale,2))+' ; '
                +str(np.round(np.quantile(samples,0.75)/scale,2)) , ha='right', va='center', transform=ax.transAxes, fontsize=fs)
        ax.text(1, 0.7, 'total data: ' + str(len(samples)), ha='right', va='center', transform=ax.transAxes, fontsize=fs) 
    
    #chequeamos que no haya numeros raros
        if True:
            neg = sum([1 for i in samples if i < bins[1][0]]) 
            pos =  sum([1 for i in samples if i > bins[1][-1]])
            text_borders = '< lower bound :'+str(neg) + '\n > upper bound :'+str(pos)
            print(text_borders)
            #ax.text(1, 0.75, text_borders, ha='right', va='center', transform=ax.transAxes, fontsize=fs) 



#%%
def plot_dt_square(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    #xlim = [-5+beg,T+5] ; ylim = [-0.02,2.02] ;         
    Cols = len(ref.groupby(['D'])) ;
    Rows = len(ref.groupby(['alpha'])) ; 
    colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)

    
    for col,(D,col_) in  enumerate(ref.groupby(['D'])):
        for row, (alpha,row_)  in  enumerate(col_.groupby(['alpha'])):
            delta= np.round(alpha/col_.omega.unique()[0],4)  
            ax = axs[row,col]; ax.grid(False);

            if row == 0:
                text = 'D = ' + str(np.round(D,5))
                ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if col == 0:
                text = 'delta = ' + str(delta)
                ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if (row == Rows-1) and (col == 0): 
                ax.set_ylabel('dt (min)', fontsize=30);
                ax.set_xlabel('probability density (1/min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.2);
                ax.yaxis.set_label_coords(-0.1, 0.5)
            
            # download data
            DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
  
            if len(DT) > 0:
                    
                bins = ax.hist(DT,bins=np.linspace(0,20,42),density=True,alpha=1,linewidth=1,color = colors[col]); 
                #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                compute_st_values(ax,DT,bins,1,20)   
            else:
                print(delta,D,"no data")
            
            ax.set_ylim([0,0.5]);
            ax.set_xlim([0,20])
            set_scale(ax,[0,5,10,15,20], [0,0.5])
            ax.set_xticklabels([0,5,10,15,20])
            ax.set_yticklabels([0,0.5])
            ax.tick_params(labelsize=20)


    plt.savefig(save_path_name + 'dt_hist_square.pdf', format='pdf')
    plt.close()
    return(0)
    
def plot_ipi_square(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    #xlim = [-5+beg,T+5] ; ylim = [-0.02,2.02] ;         
    Cols = len(ref.groupby(['D'])) ;
    Rows = len(ref.groupby(['alpha'])) ; 
    colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)

    
    for col,(D,col_) in  enumerate(ref.groupby(['D'])):
        for row, (alpha,row_)  in  enumerate(col_.groupby(['alpha'])):
            delta= np.round(alpha/col_.omega.unique()[0],4)  
            ax = axs[row,col]; ax.grid(False);

            if row == 0:
                text = 'D = ' + str(np.round(D,5))
                ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if col == 0:
                text = 'delta = ' + str(delta)
                ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if (row == Rows-1) and (col == 0): 
                ax.set_ylabel('IPI (min)', fontsize=30);
                ax.set_xlabel('probability density (1/min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.2);
                ax.yaxis.set_label_coords(-0.1, 0.5)
            
            # download data
            DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
  
            if len(DT) > 0:
                    
                bins = ax.hist(IPI,bins=np.linspace(0,40,84),density=True,alpha=1,linewidth=1,color = colors[col]); 
                #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                compute_st_values(ax,IPI,bins,1,20)   
            else:
                print(delta,D,"no data")
            
            ax.set_ylim([0,0.2]);
            ax.set_xlim([0,40])
            set_scale(ax,[0,10,20,30,40], [0,0.2])
            ax.set_xticklabels([0,10,20,30,40])
            ax.set_yticklabels([0,0.2])
            ax.tick_params(labelsize=20)


    plt.savefig(save_path_name + 'IPI_hist_square.pdf', format='pdf')
    plt.close()
    return(0)

def plot_dm_square(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    #xlim = [-5+beg,T+5] ; ylim = [-0.02,2.02] ;         
    Cols = len(ref.groupby(['D'])) ;
    Rows = len(ref.groupby(['alpha'])) ; 
    colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)

    
    for col,(D,col_) in  enumerate(ref.groupby(['D'])):
        for row, (alpha,row_)  in  enumerate(col_.groupby(['alpha'])):
            delta= np.round(alpha/col_.omega.unique()[0],4)  
            ax = axs[row,col]; ax.grid(False);

            if row == 0:
                text = 'D = ' + str(np.round(D,5))
                ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if col == 0:
                text = 'delta = ' + str(delta)
                ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if (row == Rows-1) and (col == 0): 
                ax.set_ylabel('dm (min)', fontsize=30);
                ax.set_xlabel('probability density (1/min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.2);
                ax.yaxis.set_label_coords(-0.1, 0.5)
            
            # download data
            DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
  
            if len(DT) > 0:
                    
                bins = ax.hist(dm,bins=np.linspace(0,40,84),density=True,alpha=1,linewidth=1,color = colors[col]); 
                #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                compute_st_values(ax,dm,bins,1,20)   
            else:
                print(delta,D,"no data")
            
            ax.set_ylim([0,0.4]);
            ax.set_xlim([0,30])
            set_scale(ax,[0,10,20,30], [0,0.4])
            ax.set_xticklabels([0,10,20,30])
            ax.set_yticklabels([0,0.4])
            ax.tick_params(labelsize=20)


    plt.savefig(save_path_name + 'dm_hist_square.pdf', format='pdf')
    plt.close()
    return(0)

def plot_joint_duration_square(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    #xlim = [-5+beg,T+5] ; ylim = [-0.02,2.02] ;         
    Cols = len(ref.groupby(['D'])) ;
    Rows = len(ref.groupby(['alpha'])) ; 
    colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
    colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

    fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)

    
    for col,(D,col_) in  enumerate(ref.groupby(['D'])):
        for row, (alpha,row_)  in  enumerate(col_.groupby(['alpha'])):
            delta= np.round(alpha/col_.omega.unique()[0],4)  
            ax = axs[row,col]; ax.grid(False);

            if row == 0:
                text = 'D = ' + str(np.round(D,5))
                ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if col == 0:
                text = 'delta = ' + str(delta)
                ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
            if (row == Rows-1) and (col == 0): 
                ax.set_ylabel('joint duration (min)', fontsize=30);
                ax.set_xlabel('probability density (1/min)', fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.2);
                ax.yaxis.set_label_coords(-0.1, 0.5)
            
            # download data
            DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
  
            if len(DT) > 0:
                    
                bins = ax.hist(joint_duration,bins=np.linspace(0,20,42),density=True,alpha=1,linewidth=1,color = colors[col]); 
                #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                compute_st_values(ax,joint_duration,bins,1,20)   
            else:
                print(delta,D,"no data")
            
            ax.set_ylim([0,0.5]);
            ax.set_xlim([0,20])
            set_scale(ax,[0,5,10,15,20], [0,0.5])
            ax.set_xticklabels([0,5,10,15,20])
            ax.set_yticklabels([0,0.5])
            ax.tick_params(labelsize=20)


    plt.savefig(save_path_name + 'joint_duration_hist_square.pdf', format='pdf')
    plt.close()
    return(0)
#%%
def plot_dt_square_ou(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    for alpha,ref in ref_.groupby(['alpha0']):
        
        Rows = len(ref.groupby(['tau'])) ; 
        Cols = len(ref.groupby(['sigma'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    
        
        for col,(sigma,col_) in  enumerate(ref.groupby(['sigma'])):
            for row, (tau,row_)  in  enumerate(col_.groupby(['tau'])):
                delta= np.round(alpha/col_.omega.unique()[0],4)  
                ax = axs[row,col]; ax.grid(False);
    
                if row == 0:
                    text = 'sigma = ' + str(np.round(sigma,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'tau = ' + str(np.round(tau,5))
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel('dt (min)', fontsize=30);
                    ax.set_xlabel('probability density (1/min)', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.2);
                    ax.yaxis.set_label_coords(-0.1, 0.5)
                
                # download data
                DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
      
                if len(DT) > 0:
                        
                    bins = ax.hist(DT,bins=np.linspace(0,20,42),density=True,alpha=1,linewidth=1,color = colors[col]); 
                    #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                    compute_st_values(ax,DT,bins,1,20)   
                else:
                    print(delta,sigma,tau,"no data")
                
                ax.set_ylim([0,0.5]);
                ax.set_xlim([0,20])
                set_scale(ax,[0,5,10,15,20], [0,0.5])
                ax.set_xticklabels([0,5,10,15,20])
                ax.set_yticklabels([0,0.5])
                ax.tick_params(labelsize=20)
    
    
        plt.savefig(save_path_name + 'dt_hist_square_ou_'+str(delta)+'.pdf', format='pdf')
        plt.close()
    return(0)
    
def plot_ipi_square_ou(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    for alpha,ref in ref_.groupby(['alpha0']):
        
        Rows = len(ref.groupby(['tau'])) ; 
        Cols = len(ref.groupby(['sigma'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    
        
        for col,(sigma,col_) in  enumerate(ref.groupby(['sigma'])):
            for row, (tau,row_)  in  enumerate(col_.groupby(['tau'])):
                delta= np.round(alpha/col_.omega.unique()[0],4)  
                ax = axs[row,col]; ax.grid(False);
    
                if row == 0:
                    text = 'sigma = ' + str(np.round(sigma,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'tau = ' + str(np.round(tau,5))
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)

                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel('IPI (min)', fontsize=30);
                    ax.set_xlabel('probability density (1/min)', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.2);
                    ax.yaxis.set_label_coords(-0.1, 0.5)
                
                # download data
                DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
      
                if len(DT) > 0:
                        
                    bins = ax.hist(IPI,bins=np.linspace(0,40,84),density=True,alpha=1,linewidth=1,color = colors[col]); 
                    #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                    compute_st_values(ax,IPI,bins,1,20)   
                else:
                    print(delta,sigma,tau,"no data")
                
                ax.set_ylim([0,0.2]);
                ax.set_xlim([0,40])
                set_scale(ax,[0,10,20,30,40], [0,0.2])
                ax.set_xticklabels([0,10,20,30,40])
                ax.set_yticklabels([0,0.2])
                ax.tick_params(labelsize=20)
    
    
        plt.savefig(save_path_name + 'IPI_hist_square_ou_'+str(delta)+'.pdf', format='pdf')
        plt.close()
    return(0)

def plot_dm_square_ou(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
    for alpha,ref in ref_.groupby(['alpha0']):
        
        Rows = len(ref.groupby(['tau'])) ; 
        Cols = len(ref.groupby(['sigma'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
###############################################################################
### Figure
###############################################################################    

        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    
        
        for col,(sigma,col_) in  enumerate(ref.groupby(['sigma'])):
            for row, (tau,row_)  in  enumerate(col_.groupby(['tau'])):

                delta= np.round(alpha/col_.omega.unique()[0],4)  
                ax = axs[row,col]; ax.grid(False);
    
                if row == 0:
                    text = 'sigma = ' + str(np.round(sigma,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'tau = ' + str(np.round(tau,5))
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel('dm (min)', fontsize=30);
                    ax.set_xlabel('probability density (1/min)', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.2);
                    ax.yaxis.set_label_coords(-0.1, 0.5)
                
                # download data
                DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
      
                if len(DT) > 0:
                        
                    bins = ax.hist(dm,bins=np.linspace(0,40,84),density=True,alpha=1,linewidth=1,color = colors[col]); 
                    #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                    compute_st_values(ax,dm,bins,1,20)   
                else:
                    print(delta,sigma,tau,"no data")
                
                ax.set_ylim([0,0.4]);
                ax.set_xlim([0,30])
                set_scale(ax,[0,10,20,30], [0,0.4])
                ax.set_xticklabels([0,10,20,30])
                ax.set_yticklabels([0,0.4])
                ax.tick_params(labelsize=20)
    
    
        plt.savefig(save_path_name + 'dm_hist_square_ou_'+str(delta)+'.pdf', format='pdf')
        plt.close()
    return(0)

def plot_joint_duration_square_ou(dt,d,description_file,data_folder,save_path_name):
    '''
   plottea la grid de dt

    '''
    plt.close('all')
######## Getting data information
    ref_ = pd.read_excel(description_file,sheet_name='File_references')
    ref_.set_index('Unnamed: 0',inplace=True);
        
###############################################################################
### Plotting parameters
###############################################################################    
    for alpha,ref in ref_.groupby(['alpha0']):
        
        Rows = len(ref.groupby(['tau'])) ; 
        Cols = len(ref.groupby(['sigma'])) ; 
        colors =  sns.color_palette(sns.color_palette("viridis",Cols*1))
        colors =  colors[::1]
####################################################
### Figure
###############################################################################    

        fig, axs = plt.subplots(Rows, Cols, sharex=True, sharey=True, figsize=(8.27*5, 11.69*2))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.1, hspace=0.1)
    
        
        for col,(sigma,col_) in  enumerate(ref.groupby(['sigma'])):
            for row, (tau,row_)  in  enumerate(col_.groupby(['tau'])):
                delta= np.round(alpha/col_.omega.unique()[0],4)  
                ax = axs[row,col]; ax.grid(False);
    
                if row == 0:
                    text = 'sigma = ' + str(np.round(sigma,5))
                    ax.text(0.9, 1.05, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if col == 0:
                    text = 'tau = ' + str(np.round(tau,5))
                    ax.text(-0.2, 0.9, text , ha='center', va='center', transform=ax.transAxes, fontsize=25)
                if (row == Rows-1) and (col == 0): 
                    ax.set_ylabel('joint duration (min)', fontsize=30);
                    ax.set_xlabel('probability density (1/min)', fontsize=30)
                    ax.xaxis.set_label_coords(0.5, -0.2);
                    ax.yaxis.set_label_coords(-0.1, 0.5)
                
                # download data
                DT,IPI,joint_duration,dm = download_quantifiers(row_,data_folder,dt,d)
      
                if len(DT) > 0:
                        
                    bins = ax.hist(joint_duration,bins=np.linspace(0,20,42),density=True,alpha=1,linewidth=1,color = colors[col]); 
                    #tune_plot(ax,'dt (min)','probability density (1/min)',[0,20],1,[0,0.4],1,30,20)
                    compute_st_values(ax,joint_duration,bins,1,20)   
                else:
                    print(delta,sigma,tau,"no data")
                
                ax.set_ylim([0,0.5]);
                ax.set_xlim([0,20])
                set_scale(ax,[0,5,10,15,20], [0,0.5])
                ax.set_xticklabels([0,5,10,15,20])
                ax.set_yticklabels([0,0.5])
                ax.tick_params(labelsize=20)
    
    
        plt.savefig(save_path_name + 'joint_duration_hist_square_ou_'+str(delta)+'.pdf', format='pdf')
        plt.close()
    return(0)

#%%

def plot_2d_quantifiers(dt,T,d,description_file,data_folder,save_path_name):
######## Getting data information
    plt.rcdefaults();
    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    dt_quart = [6,8.33]; IPI_quart = [8,18.67];mean_activity_sigma = [32-.675*3,32+.675*3] #mean, sigma
    quartiles = [dt_quart,IPI_quart,IPI_quart,mean_activity_sigma] 
    
    vmM_dt = [1,10];vmM_IPI = [6,20]; vmM_activity = [0,100]
    v_values = [vmM_dt,vmM_IPI,vmM_IPI,vmM_activity]    
    
    name = ['dt','IPI','FPT','mactivity']


    
    for omega,ref_ in  ref.groupby(['omega']):

###############################################################################
### Figure
###############################################################################    

       
        df_dt,df_ipi,df_fpt,df_activity = create_df(ref_,data_folder,dt,T,d)
        
        for i,df in enumerate([df_dt,df_ipi,df_fpt,df_activity]):
            plt.close()
            fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8.27, 11.69))
            fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.3, hspace=0.3)        #df_dt,df_ipi,df_fpt =  create_df(ref,data_folder,dt,d)
 
            print(df)
            q_m_ ,q_M_= quartiles[i]
            mask = ((q_m_ <= df) & (df <= q_M_)) #.replace(False,np.nan)
            
            
            axs[1,1].imshow(mask,origin='lower',alpha=1,cmap='Greys',interpolation='none')
            axs[1,1].imshow(df,origin='lower',alpha = 0.3,vmin=v_values[i][0],vmax=v_values[i][1],cmap="viridis_r")
            axs[1,1].axhline(len(df.index)//2,linestyle='dashed',color='black')

            im = axs[1,0].imshow(df,origin='lower',alpha = 1,vmin=v_values[i][0],vmax=v_values[i][1],cmap="viridis_r")
            axs[1,0].axhline(len(df.index)//2,linestyle='dashed',color='black')
            #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            #fig.colorbar(im, cax=cbar_ax)
            cbar = plt.colorbar(im)
            
            for ax in [axs[1,0],axs[1,1]]:
                ax.tick_params(axis='both', direction='out')
                ax.set_xticks(range(0,len(df.columns),10))
                ax.set_xticklabels(df.columns[::10])
                
                ax.set_xlabel('D', fontsize=10)
                ax.set_ylabel('delta', fontsize=10)
                ax.set_yticks(range(0,len(df.index),10))
                ax.set_yticklabels([np.round(y/omega,2) for y in df.index[::10]])
                
            plt.savefig(save_path_name + str(omega)+'_'+name[i]+'_2dplot.pdf', format='pdf')
                

#%% superposition plotting
            
def plot_2d_superposition(dt,T,d,description_file,data_folder,save_path_name):
######## Getting data information
    plt.rcdefaults(); plt.close()

    ref = pd.read_excel(description_file,sheet_name='File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    
    dt_quart = [6,8.33]; IPI_quart = [8,18.67];mean_activity_sigma = [32-.675*3,32+.675*3] #mean, sigma
    quartiles = [dt_quart,IPI_quart,IPI_quart,mean_activity_sigma] 
    #vmin = [1,6,6,0]; vmax=[10,20,20,100]
    
         
    #name = ['dt','IPI','FPT','mActivity']; 
    cmaps = ['Reds','Greens','Greens', 'Blues']
    
    for omega,ref_ in  ref.groupby(['omega']):

###############################################################################
### data
###############################################################################    

        df_dt,df_ipi,df_fpt,df_activity= create_df(ref_,data_folder,dt,T,d)
        masks = []
        for i,df in enumerate([df_dt,df_ipi,df_fpt,df_activity]):
 
            q_m ,q_M= quartiles[i]
            masks.append(((q_m <= df) & (df <= q_M)))
            
###############################################################################
### figure
###############################################################################   
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8.27, 11.69))
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.1, right=0.99, wspace=0.3, hspace=0.3)        
        
        #save masks
        
        with pd.ExcelWriter("masks.xlsx") as writer:  
            (masks[0] & masks[2] & masks[3]).to_excel(writer,sheet_name='superposition')
            (masks[0]).to_excel(writer,sheet_name='duration')
            (masks[1]).to_excel(writer,sheet_name='IPI')
            (masks[2]).to_excel(writer,sheet_name='FPT')
            (masks[3]).to_excel(writer,sheet_name='mean activity')
            
        for i,mask in enumerate(masks):
            axs[1,1].imshow(mask,origin='lower',alpha=0.3,cmap=cmaps[i],interpolation='none')
        
        axs[1,1].axhline(len(df_dt.index)//2,linestyle='dashed',color='black')
          
            
        for ax in [axs[1,0],axs[1,1]]:
            ax.tick_params(axis='both', direction='out')
            ax.set_xticks(range(0,len(df.columns),10))
            ax.set_xticklabels(df.columns[::10])
            
            ax.set_xlabel('D', fontsize=10)
            ax.set_ylabel('delta', fontsize=10)
            ax.set_yticks(range(0,len(df.index),10))
            ax.set_yticklabels([np.round(y/omega,2) for y in df.index[::10]])
                
        plt.savefig(save_path_name + str(omega)+'_2dsuperposition.pdf', format='pdf')


