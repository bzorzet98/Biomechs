#%% """Procesamiento de un unico registro, utilizando las clases implementadas"""
from REGISTER import REGISTER 
from BIOMECH_MODEL import BIOMECH_MODEL
import numpy as np
import matplotlib.pyplot as plt
import function_biomechs as fb
#%% Carga del registro
path='Dollyo Chagi/Izquierda/PatadaDollyoChagiPI.c3d'
patada=BIOMECH_MODEL(path)

#%% Filtrado
patada.Filter_Markers(10,4)
#%% Extraemos los eventos
list(patada.reg['events_times'])
#%% Calculamos los sistemas coordenados locales
patada.Calc_Coordinates_Systems()

#%% Graficamos los sistemas
#patada.Plot_Sist_Coor()

#%% Calculamos los sistemas anatomicos
patada.Calc_Anatomic_System()

#%% Calculamos los Angulos articulares
patada.Joint_Angle()
#patada.Plot_Joints_Angles()
# %%Calculamos las velocidades de los angulos articulares y las ploteamos
patada.Joint_Angle_Vel()
patada.Plot_Joints_Angles_Velocities()
# %%Calculamos las velocidades de los angulos articulares y ploteamos con 4 puntos de 
patada.Joint_Angle_Vel(method="4points")
patada.Plot_Joints_Angles_Velocities(method="4points")
#%% Determinamos los angulos de Euler de los segmentos
segmentos=patada.segments
if 'PI' in path:
    rotacion=['131','232','131','232','131','232','131']
elif 'PD' in path: 
    rotacion=['131','131','232','131','232','131','232']
patada.Ang_Euler(segmentos,rotacion)

patada.Plot_Euler_Angles()
a=20
