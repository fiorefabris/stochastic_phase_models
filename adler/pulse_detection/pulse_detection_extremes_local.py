import numpy as np
import matplotlib.pyplot as plt
from main import time_evolution_
import matplotlib.gridspec as gridspec


#%%

def search_extremes(X,TH,W):
    '''
    search_extremes(X,TH,W)
    Search for maxima and minima on a time series X. 
    
    
    Searchs for the local maxima and minima on the time series x that are grater (maxima) than a treshold or 
    smaller (minima) than - threshold. 
    The extremes have at least W separation. 
    
    
    Parameters
    ----------
    X : list
        amplitude values of the time series
    TH : float
        threshold value
    W : minimum distance between maxima / minima
    
    
    Returns
    -------
    MAX : list 
        position of the maxima
    MIN : list
        position of the minima 
    '''
    

#Comments: no se si esta bien sacarle el ultimo elemento; 
#pensar como incluir los bordes: NO LOS INCLUYO PORQUE NO TENGO ESTADISTICA DE PICOS AHI
    
### Initializing variables
    MAX = [0] ; MIN = [0]
    M = 0; m = 0; 

    
    for j, x_j in enumerate(X[1:-1]):
        #Maxima 
        if x_j > TH : 
            #si x_j supera el treshold
            for i in range(1,W):
                
                x_p = X[j+1]
                x_m = X[j-1]
            
            if (x_j >= x_p) and (x_j >= x_m ) : #Candidato a máximo
                if M == 0 or (M >=1 and ((j-MAX[M]) > W)) : 
                    #Si es el primer máximo o tiene una distancia más grande que W con
                    #el maximo anterior, anotalo
                    M = M + 1 ; MAX.append(j) 
                else: 
                    
                    if (x_j < X[MAX[M]]) : 
                        #reject new maximum
                        pass 
                    else: 
                        # overwrite the old maxima by leaving M unchanged
                        MAX[M] = j 
        #Minima
        if x_j < (-TH) :
            x_p = X[j+1]
            x_m = X[j-1]
        
            if (x_j <= x_p) and (x_j <= x_m) : #Candidato a mínimo
            
                #Now test for adjacent minimum, and delete the higher one.   
                if m == 0 or (m >=1 and ((j-MIN[m]) > W)):
                    
                    m = m + 1; MIN.append(j) 
    
                else: #old minimum exists
                    if x_j > X[MIN[m]] : #reject new minimum
                        pass
                    else:  # Overwrite old one by leaving m unchanged
                        MIN[m] = j
    MAX.pop(0);MIN.pop(0)
    return(MAX,MIN)
    
#%%
def filter_extremes(MAX,MIN,X):
    # filtra los extremos para que siempre haya un maximo antes que un minimo al principio de la serie temporal,
    # y un minimo despues de un maximo al final de la serie temporal
    if MAX[0] < MIN[0]:
        for M1,M2 in zip(MAX[:-1],MAX[1:]):
            aux_list = list(filter(lambda x: (x > M1 and x < M2), MIN))
            while len(aux_list) > 1:
                aux_list_values = [X[i] for i in aux_list]
                MIN.remove(aux_list[np.argmax(aux_list_values)])
                aux_list = list(filter(lambda x: (x > M1 and x < M2), MIN))
            if len(aux_list) == 0:
                aux_list_values = [X[i]  for i in [M1,M2]]
                MAX.remove([M1,M2][np.argmax(aux_list_values)])
    else:
        for m1,m2 in zip(MIN[:-1],MIN[1:]):
            aux_list = list(filter(lambda x: (x > m1 and x < m2), MAX))
            while len(aux_list) > 1:
                aux_list_values = [X[i]  for i in aux_list]
                MAX.remove(aux_list[np.argmax(aux_list_values)])
                aux_list = list(filter(lambda x: (x > m1 and x < m2), MAX))
            if len(aux_list) == 0:
                aux_list_values = [X[i]  for i in [m1,m2]]
                MIN.remove([m1,m2][np.argmax(aux_list_values)])    
    return(MAX,MIN)
            

#%%#
def get_fixed_points(alpha,omega):
    #calcula los puntos fijos de el sistema de adler determinista para un alpha y omega dados
    res = np.arcsin(-1/alpha)
    PFE = -res + np.pi  #np.sin(PFE)
    PFI = (2*np.pi + res)  #np.sin(PFI)
    return(PFE,PFI)


#%%

def get_left_minima(MAX,PFE,PFI,theta):
    #busca donde comienzan los pulsos
    left_minima = []; 
    n0 = theta[MAX[0]] // (2*np.pi)
    
    for n,M in zip(np.arange(n0,len(MAX)+n0),MAX):
       
        while not (PFI < (theta[M] - n*2*np.pi) < 5*np.pi/2):
            n = n - 1
            
        theta_ = theta - n*2*np.pi

        for i, ang in enumerate(theta_[:M][::-1]): #empieza uno después del maximo
            if ( ang - PFI < 0 ):
                left_minima.append(M-i); 
                break
            else: 
                pass
    return(left_minima)
            
def get_rigth_minima(MIN,PFE,PFI,theta):
    #busca donde terinan los pulsos
    right_minima = [];
    n0 = theta[MIN[0]] // (2*np.pi)
    
    for n,m in zip(np.arange(n0,len(MIN)+n0),MIN):
        
        while not (np.pi/2 < (theta[m] - n*2*np.pi) < PFE):
            n = n - 1
        
        theta_ = theta - n*2*np.pi

        for i, ang in enumerate(theta_[m:]): #empieza uno después del maximo
            #print(ang-PFE)
            if ( ang - PFE > 0 ):
                right_minima.append(m+i); 
                break
            else: 
                pass
    return(right_minima)
#%%
def clean_extremes(left_minima,right_minima,MAX,MIN):
    
    while left_minima[0] > right_minima[0]:
        right_minima.pop(0)
    while left_minima[-1] > right_minima[-1]:
        left_minima.pop(-1)
    
    while MAX[0] < left_minima[0]:
        MAX.pop(0)
    while MAX[-1] > right_minima[-1]:
        MAX.pop(-1)
    
    while MIN[0] < left_minima[0]:
        MIN.pop(0)
    while MIN[-1] > right_minima[-1]:
        MIN.pop(-1)
    return(left_minima,right_minima,MAX,MIN)
#%%
    
def get_quantifiers(left_minima,right_minima,MAX,MIN):
    IPI = []; dt = []; dm = []; joint_duration = []
    for left,right in zip(left_minima,right_minima):
        print(right-left)
        dt.append(right - left)
    for M1,M2 in zip(MAX[:-1],MAX[1:]):
        IPI.append(M2-M1)
        right_filter = list(filter(lambda x: x > M1, right_minima))[0]
        left_filter = list(filter(lambda x: x < M2, left_minima))[-1]
        joint_duration.append((right_filter-M1) + (M2 -left_filter))
    for right,left in zip(right_minima[:-1],left_minima[1:]):
        dm.append(left-right)
    return(dt,IPI,dm,joint_duration)
#%% get_rigth_minima(MIN,PFE,PFI,theta)

save_folder = '/home/fabris/Documents/Dyncode/Dyncode_simulations/pulse_detection/extremes/'
omega =2*np.pi/7
alpha = [1.01]
for a in alpha:
    theta = time_evolution_(0.001,10000,1,omega,a *omega ,0.005)
    
    fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
    gs_row = gridspec.GridSpec(nrows=3, ncols=1, figure=fig, wspace=0.4, hspace=0.3)
    
    ax = plt.subplot(gs_row[0]) 
    ax.plot(np.cos(theta))
    labels = [item/1000 for item in ax.get_xticks()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('min')
    ax.set_ylabel(r'$\cos(\theta)$')
    
    PFE , PFI = get_fixed_points(a,omega)
    MAX, MIN = search_extremes(np.cos(theta),0.9,1000)
    MAX, MIN = filter_extremes(MAX,MIN,np.cos(theta))
    left_minima = get_left_minima(MAX,PFE,PFI,theta)
    right_minima = get_rigth_minima(MIN,PFE,PFI,theta)
    
    ax = plt.subplot(gs_row[1]) 
    ax.plot(np.cos(theta))
    ax.plot(MAX,np.cos(theta)[MAX],'o',color = 'red',markersize = 2)
    ax.plot(MIN,np.cos(theta)[MIN],'o',color = 'blue',markersize = 2)
    ax.plot(left_minima,np.cos(theta)[left_minima],'o',color = 'darkgray',markersize = 2)
    ax.plot(right_minima,np.cos(theta)[right_minima],'o',color='black',markersize = 2)
    labels = [item/1000 for item in ax.get_xticks()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('min')
    ax.set_ylabel(r'$\cos(\theta)$')
    
    left_minima,right_minima,MAX,MIN = clean_extremes(left_minima,right_minima,MAX,MIN)
    
    ax = plt.subplot(gs_row[2]) 
    ax.plot(np.cos(theta))
    ax.plot(MAX,np.cos(theta)[MAX],'o',color = 'red',markersize = 2)
    ax.plot(MIN,np.cos(theta)[MIN],'o',color = 'blue',markersize = 2)
    ax.plot(left_minima,np.cos(theta)[left_minima],'o',color = 'darkgray',markersize = 2)
    ax.plot(right_minima,np.cos(theta)[right_minima],'o',color='black',markersize = 2)
    labels = [item/1000 for item in ax.get_xticks()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('min')
    ax.set_ylabel(r'$\cos(\theta)$')

    plt.savefig(save_folder + 'time_series_alpha_'+str(a)+'.pdf', format='pdf')

#%%
#1000 frames 1 minuto

dt,IPI,dm,joint_duration = get_quantifiers(left_minima,right_minima,MAX,MIN)
fig = plt.figure(constrained_layout = False,figsize=(8.27, 11.69),linewidth=0.5)
gs_row = gridspec.GridSpec(nrows=3, ncols=3, figure=fig, wspace=0.4, hspace=0.3)
ax = plt.subplot(gs_row[0,0])            
ax.hist(dt)
labels = [item/1000 for item in ax.get_xticks()]
ax.set_xticklabels(labels)
ax.set_ylabel('counts', fontsize=8);
ax.set_xlabel('dt (min)', fontsize=8)

ax = plt.subplot(gs_row[0,1])            
ax.hist(IPI)
labels = [item/1000 for item in ax.get_xticks()]
ax.set_xticklabels(labels)
ax.set_ylabel('counts', fontsize=8);
ax.set_xlabel('IPI (min)', fontsize=8)

ax = plt.subplot(gs_row[1,0])            
ax.hist(dm)
labels = [item/1000 for item in ax.get_xticks()]
ax.set_xticklabels(labels)
ax.set_ylabel('counts', fontsize=8);
ax.set_xlabel('silent intervals (min)', fontsize=8)

ax = plt.subplot(gs_row[1,1])            
ax.hist(joint_duration)
labels = [item/1000 for item in ax.get_xticks()]
ax.set_xticklabels(labels)
ax.set_ylabel('counts', fontsize=8);
ax.set_xlabel('joint duration (min)', fontsize=8)


ax = plt.subplot(gs_row[2,0])            
ax.plot(dm,joint_duration,'o')
ax.plot(np.arange(np.max(joint_duration)),np.arange(np.max(joint_duration))*2,linestyle='-',color='black')
ax.set_xlim([0,100000]);ax.set_ylim([0,100000])
labels = [item/1000 for item in ax.get_xticks()]
ax.set_xticklabels(labels)
labels = [item/1000 for item in ax.get_yticks()]
ax.set_yticklabels(labels)
ax.set_xlabel('silence (min)', fontsize=8);
ax.set_ylabel('joint duration (min)', fontsize=8)


fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.07, hspace=0.07)
fig.savefig(save_folder + 'histograms_alpha_'+str(a)+'.pdf', format='pdf')
  
#%% estas son cosas para correr con un solo d en pulse:detection_main_2

def main_pulse_detection_(data_folder,save_path_name,delta,tuple_):
    #ES UNA funcion auxiliar para paralelizar
   # (i,D,_),row = tuple_[0],tuple_[1]
    j = 0
    (i,D),row = tuple_[0],tuple_[1].iloc[[j]]
    omega =  row.omega.unique()[0]
    alpha = np.round(i/omega,4)  
    file_name =  str(int(row.number))+'_'+str(int(row.order))+'.pkl'
    theta = download_theta(file_name,data_folder)[::delta]
    while len(theta) == 0 :
        j = j + 1
        theta = download_theta(file_name,data_folder)[::delta]
    main_pulse_detection(theta,alpha,omega,D,save_path_name,file_name)
    return(1)
    
    
def compute_pulse_detection(delta,description_file,data_folder,save_path_name):
    #esta es la funcion que le calcula a cada seite temporal suspulsos
    #delta es la resolucion de la serie temporal
    ref = pd.read_excel(description_file,sheet_name= 'File_references')
    ref.set_index('Unnamed: 0',inplace=True);
    pool = mp.Pool(processes= ceil(mp.cpu_count()))
    
    tuple_ = ref.groupby(['alpha','D'])
    main_pulse_detection__ = partial(main_pulse_detection_,data_folder,save_path_name,delta)
    pool.map(main_pulse_detection__,tuple_)
    pool.close()
    pool.join()
    return (2)


