#import pickle5 as pickle 
import pickle
from sklearn.linear_model import LinearRegression
import numpy as np
from itertools import groupby

def count_iterable(iterator):
    return sum(1 for i in iterator)


# =============================================================================
# from plotting_main
# =============================================================================
def load_file(filename):
    infile = open(filename,'rb')
    output = pickle.load(infile)
    infile.close()
    return(output)



# =============================================================================
# from consecutive_main
# =============================================================================

class consecutive_cumulative:
    def __init__(self,df):
        self.df = df
        
    def joint_duration(self,data):
        values_dt_mixed = []
        cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
        return(values_dt_mixed)

    def is_consecutive_cell(self,data):
        #te devuelve par de pulsos consecutivos, nivel celula
        y = data['IPI'].dropna().values-self.joint_duration(data) #silence 
        x = self.joint_duration(data) #joint duration
        consecutive = []
        for i in range(len(x)):
            if x[i]*0.5 >= y[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        
    def is_consecutive_box(self):
        df = self.df
        #Te da estadística de pulsos en toda la poblacion PARES DE PULSOS, nivel poblacion
        box = []
        for cells, data in df.groupby(level='cell'):
            box.append(self.is_consecutive_cell(data))
        return(box)
    
    def raise_order_consecutiveness(self,box):
        #Te da lista vacia cuando no hay mas pares de pulsos. 
        # Calcula un nivel más de consecutividad
        new_box = []
        for consecutive_ in box:
            new_consecutive_ = []
            for i,j in zip(consecutive_[:-1],consecutive_[1:]):
                new_consecutive_.append(i*j)
            new_box.append(new_consecutive_)
        return(new_box)
                
    def get_consecutive_trains_of_pulses(self):
        #te devuelve un array en donde cada elemento tiene la cantidad de trenes de pulsos sobre el total de la población
        df  = self.df
        box = self.is_consecutive_box() #pares de pulsos
        box_plot_consecutive = [df.amp_peaks.count(),np.sum([np.sum(l) for l in box])]
        
        while box_plot_consecutive[-1] > 0:
            box= self.raise_order_consecutiveness(box)
            box_plot_consecutive.append(np.sum([np.sum(l) for l in box])) #cada elemento de la lista new_box(box) es una célula
        return(box_plot_consecutive[:-1]) # el último elemento es el que se va a cero, por eso lo sacamos
    
    def get_consecutive_trains_of_pulses_cells(self):
        #te devuelve un array en donde cada elemento tiene una lista donde en cada 
        # lugarcito tiene la cantidad de trenes de la célula
        df = self.df
        #cada lugarcito es la suma en una celula
        box = self.is_consecutive_box() #pares de pulsos
        box_plot_consecutive = [df.amp_peaks.groupby(level='cell').count().values,[np.sum(l) for l in box]]
        while np.sum(box_plot_consecutive[-1]) > 0:
            box = self.raise_order_consecutiveness(box)
            box_plot_consecutive.append([np.sum(l) for l in box])
            box = box
        return(box_plot_consecutive)
        
class consecutive_non_cumulative:

    def __init__(self,df):
        self.df = df
        
    def joint_duration(self,data):
        values_dt_mixed = []
        cell_t_M = data[data['amp_peaks'].notna()].FRAME.values
        cell_t_m = data[data['min_'].notna()].FRAME.values
        for (t_M_izq,t_M_der) in zip(cell_t_M[:-1],cell_t_M[1:]):
            dt_raise = t_M_der - cell_t_m[cell_t_m < t_M_der][-1]
            dt_fall = cell_t_m[cell_t_m > t_M_izq][0] - t_M_izq
            mixed_dt = dt_raise + dt_fall
            values_dt_mixed.append(mixed_dt)
        return(values_dt_mixed)
        

    def silence(self,data):
        return(data['IPI'].dropna().values-self.joint_duration(data))


    def is_isolated_cell(self,MAX,joint_duration,dm):
        # solo para una celula
        # calcula el numero de picos no consecutivis (isolated)
        # tiene en cuenta las celulas de un solo pulso
        # te devuelve la suma del total de picos menos los que son parte de intervalos consecutivos de pulsos
        consecutive = self.is_consecutive_cell(joint_duration,dm)
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(MAX.count()  - count)
                    
    def is_isolated_box(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos aislados
        box = []
        for cell,data in df.groupby(level='cell'):
            box.append(self.is_isolated_cell(data.amp_peaks,self.joint_duration(data),self.silence(data)))
        return(box)
    

    def is_consecutive_cell(self,joint_duration,dm):
        # para cada celula (data), te devuelve un array con 1 representando pares de pulsos consecutivos y
        # 0 pares de pulsos no consecutivos
        consecutive = []
        for i in range(len(joint_duration)):
            if joint_duration[i]*0.5 >= dm[i]:
                consecutive.append(1)
            else:
                consecutive.append(0) 
        return(consecutive)
        
    
    def is_consecutive_box(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos consecutivos en una lista
        box = []
        for cell,data in df.groupby(level='cell'):
                box.append(self.is_consecutive_cell(self.joint_duration(data),self.silence(data)))
        return(box)
        
    
    def raise_order_consecutiveness(self,box):
        #Te da lista vacia cuando no hay mas pares de pulsos. 
        # Calcula un nivel más de consecutividad
        new_box = []
        for consecutive_ in box:
            new_consecutive_ = []
            for i,j in zip(consecutive_[:-1],consecutive_[1:]):
                new_consecutive_.append(i*j)
            new_box.append(new_consecutive_)
        return(new_box)
    
    
    def count_consecutive_pulses(self,box):
        #te devuelve una lista.- En cada lugar de la lista está el número de unos aislados en cada célula .
        count_box = []
        for consecutive in box:
            count = 0
            for i, group in groupby(consecutive):
                if i == 1:
                    if count_iterable(group) == 1:
                        count = count + 1
                    else:
                        pass
            count_box.append(count) 
        return(count_box)
    
    
    def get_consecutive_trains_of_pulses(self):
        #te da la suma sobre todo el dataset
        
        isolated = self.is_isolated_box()
        box = self.is_consecutive_box()
        
        box_plot_consecutive = [np.sum(isolated)]
        
        while np.sum([np.sum(l) for l in box]) > 0:
            box_plot_consecutive.append(sum(self.count_consecutive_pulses(box)))
            box = self.raise_order_consecutiveness(box)
            
        return(box_plot_consecutive)
    
    def get_number_of_pulses(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos aislados
        pulses = []
        for cell,data in df.groupby(level='cell'):
            pulses.append(data.amp_peaks.count())
        return(pulses)
 
        
    
    def get_consecutive_trains_of_pulses_cells(self):
        # te da la estadística por célula
        isolated = self.is_isolated_box()
        box = self.is_consecutive_box()
        
        box_plot_consecutive = [isolated]
        
        while np.sum([np.sum(l) for l in box]) > 0:
            box_plot_consecutive.append(self.count_consecutive_pulses(box,False))
            box = self.raise_order_consecutiveness(box)
        return(box_plot_consecutive)
    
    
    
    def count_consecutive_pulses_cell_number(self,MAX,joint_duration,dm):
        # solo para una celula
        # calcula el numero de picos no consecutivis (isolated)
        # tiene en cuenta las celulas de un solo pulso
        # te devuelve la suma del total de picos menos los que son parte de intervalos consecutivos de pulsos
        consecutive = self.is_consecutive_cell(joint_duration,dm)
        count = 0
        for i, group in groupby(consecutive):
            if i == 1:
                count = count + count_iterable(group) + 1
        return(count)
                    
    def count_consecutive_pulses_number(self):
        df = self.df
        #Para cada conjunto de celulas, te da la estadistica de pulsos consecutivos
        box = []
        for cell,data in df.groupby(level='cell'):
            box.append(self.count_consecutive_pulses_cell_number(data.amp_peaks,self.joint_duration(data),self.silence(data)))
        return(box)