from ezc3d import c3d
import copy
c = c3d('Dollyo Chagi/Izquierda/PatadaDollyoChagiPI_004.c3d')
print(c['parameters']['POINT']['USED']['value'][0]);  # Print the number of points used
point_data = c['data']['points']
points_residuals = c['data']['meta_points']['residuals']
analog_data = c['data']['analogs']

prueba={}
prueba['prueba15']={'prueba':15,'etc':51}




#Vemos las claves iniciales
c_list=list(c)

#Analizamos la organizacion de header
header=list(c['header'])

#Analizamos la organizacion de header points
header_points=list(c['header']['points'])#Contiene el tamaño de la muestra , la frecuencia de muestreo y el primer y ultimo frame

#Analizamos la organizacion de header analogs
header_analogs=list(c['header']['analogs'])#Contiene el tamaño de la muestra , la frecuencia de muestreo y el primer y ultimo frame

#Analizamos la organizacion de header events
header_events=list(c['header']['events'])#Contiene el tamaño de la muestra , el tiempo de los eventos y las etiquetas
header_events_eventstime=list(c['header']['events']['events_time'])#RARO LOS EVENTOS NO COINCIDEN CON LOS QUE SE MUESTRAN EN C[PARAMETERS]


#Analizamos la organización de data
data=list(c['data'])

#Analizamos la organizacion de data points
data_point=c['data']['points']
data_point_extract=data_point[0:3,0,:]
data_point_extract=data_point_extract.T
#DATA POINT CONTIENE UNA MATRIZ TRIDIMENSIONAL DONDE LA PRIMER DIMENSION CORRESPONDE A LAS
#COORDENADAS X,Y,Z, LA SEGUNDA DIMENSION CORRESPONDE A LOS MARCADORES Y LA TERCER DIMENSION A LOS FRAMES. 

#Analizamos data meta point

data_metapoints=list(c['data']['meta_points'])#Podemos obtener os resultados de los enmascaramientos del archivo

#Analizamos data analogs
data_analogs=c['data']['analogs']
#Tiene los datos de las señales analogicas


#Analizamos la organizacion de parameters
parameter=list(c['parameters'])
#Analizamos la organizacion del dict PARAMETERS EVENTS
parameter_event=list(c['parameters']['EVENT'])
parameter_event_labels=list(c['parameters']['EVENT']['LABELS'])
parameter_event_labels_value=c['parameters']['EVENT']['LABELS']['value']

parameter_event_times=list(c['parameters']['EVENT']['TIMES'])
parameter_event_times_value=list(c['parameters']['EVENT']['TIMES']['value'])

#Analizamos la organizacion del dict PARAMETERS ANALOG
parameter_analog=list(c['parameters']['ANALOG'])
parameter_analog_labels=list(c['parameters']['ANALOG']['LABELS'])
parameter_analog_labels_value=c['parameters']['ANALOG']['LABELS']['value']

#Analizamos la organizacion del dict PARAMETERS ANTROPOMETRIA
parameter_antropometria=list(c['parameters']['antropometria'])
parameter_antropometria_altura=list(c['parameters']['antropometria']['ALTURA'])
parameter_antropometria_altura_value= c['parameters']['antropometria']['ALTURA']['value']

parameter_matrizrot=list(c['parameters']['matriz_rotacion'])
parameter_matrizrot_r=list(c['parameters']['matriz_rotacion']['R'])
parameter_antropometria_altura_value= c['parameters']['antropometria']['ALTURA']['value']


#Analizamos la organizacion del dict PARAMETERS POINT
punto=list(c['parameters']['POINT'])

punto_rate=list(c['parameters']['POINT']['RATE'])
punto_rate_value=c['parameters']['POINT']['RATE']['value']

punto_metadata=list(c['parameters']['POINT']['__METADATA__'])
punto_metadata_desc=c['parameters']['POINT']['__METADATA__']['DESCRIPTION']

punto_used=list(c['parameters']['POINT']['USED'])
punto_used_descrip=c['parameters']['POINT']['USED']['description']
punto_used_value=c['parameters']['POINT']['USED']['value']

punto_scale=list(c['parameters']['POINT']['SCALE'])
punto_scale_descrp=c['parameters']['POINT']['SCALE']['description']
punto_scale_value=c['parameters']['POINT']['SCALE']['value']

punto_frames=list(c['parameters']['POINT']['FRAMES'])
punto_frames_desc=c['parameters']['POINT']['FRAMES']['description']
punto_frames_values=c['parameters']['POINT']['FRAMES']['value']
punto_frames_type=c['parameters']['POINT']['FRAMES']['type']

punto_label=list(c['parameters']['POINT']['LABELS'])
punto_label_desc=c['parameters']['POINT']['LABELS']['description']
punto_label_values=c['parameters']['POINT']['LABELS']['value']



# Moficando el c3d
a=c3d('Dollyo Chagi/Izquierda/PatadaDollyoChagiPI_004.c3d')

#Añadiendo el tamaño de la cintura a antropometria

#De esta forma podemos añadir lo que se nos ocurra
a.add_parameter("antropometria","cintura",150.0)
print(a['parameters']['antropometria']['cintura']['value'])

a.add_parameter("matricesderotacionpruebita","matriz_rot",[[1,1,1],[2,2,2],[3,3,3]])
print(a['parameters']['matricesderotacionpruebita']['matriz_rot']['value'])

#Copiamos el c3d
a.write("path_to_c3d.c3d")


#Probamos modificar los valores
a['parameters']['antropometria']['cintura']['value']=15
print(a['parameters']['antropometria']['cintura']['value'])
del a['parameters']['antropometria']['cintura']
print(a['parameters']['antropometria']['cintura']['value'])

del a['parameters']['matricesderotacionpruebita']
#LO que hay quea=copy.copy(c)
h=f[('tita','phi')]

f={}
f.keys


elementos=list(c)