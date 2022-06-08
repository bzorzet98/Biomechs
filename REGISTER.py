from ezc3d import c3d
import numpy as np
import scipy.signal as signal


class REGISTER:
    
    #Datos del archivo c3d
    arch_c3d=0
    name_c3d=''
    path_c3d=''

    #Registro
    register={}

    #Lista de marcadores
    labels_markers=''
    labels_anthropometry=''
    labels_events=''

    #Longitud de marcadores
    frames=0

    #Dictionary Labels

    def __init__(self, path=''):

        #path_c3d contiene la dirección del archivo.c3d
        if path != '': 
            self.path_c3d=path
            self.arch_c3d=c3d(self.path_c3d)    
            self.__extractdata__()
        else: 
            print ('No se ingreso dirección de archivo .c3d')

    def load_c3d(self, path): 
        self.path_c3d=path
        self.arch_c3d=c3d(self.path_c3d)    
        self.__extractdata__()

    def __extractdata__(self):
        #Extraemos la cantidad de frames
        self.frames=self.arch_c3d['parameters']['POINT']['FRAMES']['value']
        #Extraemos las etiquetas de los marcadores utilizados
        self.labels_markers=self.arch_c3d['parameters']['POINT']['LABELS']['value']
        #Extraemos las etiquetas de los datos antropometricos
        self.labels_anthropometry=list(self.arch_c3d['parameters']['antropometria'])
        self.labels_events=self.arch_c3d['parameters']['EVENT']['LABELS']['value']

        #Añadimos los valores a los diccionarios
        aux_dic={} #Definimos un diccionario vacio
        data_matrix=self.arch_c3d['data']['points']#Extramos los archivos de la matriz de valores
        self.register['markers']={}#Inicializamos el diccionario vacio
        for i in range(len(self.labels_markers)):
            aux1=data_matrix[0:3,i,:]
            aux=data_matrix[0:3,i,:].T
            b=self.labels_markers[i][11:]
            aux_dic[self.labels_markers[i][11:]]=data_matrix[0:3,i,:].T#Transponemos la matriz y la guardamos en el dic vacio
        self.labels_markers=list(aux_dic)
        self.register['markers']=aux_dic
        self.register['frequency']=self.arch_c3d['parameters']['POINT']['RATE']['value']


        self.register['anthropometry']={}
        aux_dic={}
        for i in range(1,len(self.labels_anthropometry)):
            aux_dic[self.labels_anthropometry[i]]=self.arch_c3d['parameters']['antropometria'][self.labels_anthropometry[i]]['value']
        self.register['anthropometry']=aux_dic
        aux_dic={}
        aux_dic2={}
        aux_event_times=self.arch_c3d['parameters']['EVENT']['TIMES']['value']
        for i in range(len(self.labels_events)):
            aux_dic[self.labels_events[i]]=aux_event_times[1][i]
            aux_dic2[self.labels_events[i]]=int(aux_event_times[1][i]*self.register['frequency'])
        self.register['events_times']=aux_dic
        self.register['events_frames']=aux_dic2
        self.register['rot_matrix']={}
        self.register['rot_matrix']['R']=self.arch_c3d['parameters']['matriz_rotacion']['R']['value']
        self.register['rot_matrix']['L']=self.arch_c3d['parameters']['matriz_rotacion']['L']['value']
        
    def filter_markers(self, fc, N):
        N=int(N)

        #Create a butterworth filter
        b,a=signal.butter(N, fc, fs=self.register['frequency'])

        #Filter signals
        self.register['filtered_markers']={}

        for markers in self.labels_markers:
            data=self.register['markers'][markers]
            filtered_signal=signal.filtfilt(b,a, data[:,0])
            filtered_signal=np.column_stack((signal.filtfilt(b,a, data[:,0]),signal.filtfilt(b,a, data[:,1]),signal.filtfilt(b,a, data[:,2])))
            self.register['filtered_markers'][markers]=filtered_signal


    def give_data(self,data='markers'):
        if (data in list(self.register)):
            return self.register[data]
        else:
            print('La clave ingresada no existe en el diccionario')
        