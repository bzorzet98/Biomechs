import numpy as np
import math as mt
import pandas as pd
import scipy.interpolate as interp
import matplotlib.pyplot as plt


def Euler_Angle(datos_sist, rotation='232',datos_sist_origen=[]):
    tamanio=datos_sist.shape

    if len(datos_sist_origen)==0:
        versor_I=np.zeros((tamanio[0],3));
        versor_J=np.zeros((tamanio[0],3));
        versor_K=np.zeros((tamanio[0],3));

        versor_I[:,0]=1;
        versor_J[:,1]=1;
        versor_K[:,2]=1;
    else:
        versor_I=datos_sist_origen[:,:,0]
        versor_J=datos_sist_origen[:,:,1]
        versor_K=datos_sist_origen[:,:,2]


    n=tamanio[0];
    ang_euler=np.zeros((n,3));

    if rotation=='232':
        versor_LN=np.cross(versor_J,datos_sist[:,:,2],axis=1)
        versor_LN=versor_LN/np.linalg.norm(versor_LN[:,:,1],axis=1) [:,None]


        #Calculamos el angulo de rotaci�n con respecto al eje global Y,
        #correspondiente al alfa.
        
        ang_euler[:,1]=dot(versor_I,versor_LN)
        ang_euler[:,1]=ang_euler[:,1]/np.linalg.norm(ang_euler[:,1],axis=1) [:,None]
        ang_euler[:,1]=ang_euler[:,1]*np.arccos(dot(versor_K,versor_LN))

        #Calculamos el angulo de rotacion con respecto al eje x'', obteniendo el
        #segundo sistema intermedio (xyz)'

        ang_euler[:,2]=np.arccos(dot(versor_J,datos_sist[:,:,2]));

        #Calculamos el angulo de rotacion con respecto al eje z', obteniendo la
        #rotaci�n a nuestro sistema local xyz.

        ang_euler[:,3]=dot(datos_sist[:,:,1],versor_LN)
        ang_euler[:,3]=ang_euler[:,3]/np.linalg.norm(ang_euler[:,3],axis=1)[:,None]
        ang_euler[:,3]=-1*ang_euler[:,3]*np.arccos(dot(datos_sist[:,:,3],versor_LN));
        
    elif rotation=='313':

        versor_LN=np.cross(versor_K,datos_sist[:,:,3],axis=1);
        versor_LN=versor_LN/np.linalg.norm(versor_LN,axis=1)[:,None]

        #Calculamos el angulo de rotacion con respecto al eje Z global, corresponde
        #a la primer rotacion, obteniendo el sistema intermedio (xyz)''
        aux=dot(versor_J,versor_LN)
        aux=aux/np.linalg.norm(aux,axis=1)[:,None]
        aux=aux*np.arccos(dot(versor_I,versor_LN))
        ang_euler[:,1]=aux

        #Calculamos el angulo de rotacion con respecto al eje x'', obteniendo el
        #segundo sistema intermedio (xyz)'

        aux=dot(versor_K,datos_sist[:,:,3])
        aux=np.arccos(aux)
        ang_euler[:,2]=aux
        
        #Calculamos el angulo de rotacion con respecto al eje z', obteniendo la
        #rotaci�n a nuestro sistema local xyz.
        matriz=datos_sist[:,:,2]
        aux=dot(matriz,versor_LN,2);
        aux=aux/np.linalg.norm(aux,axis=1)[:,None]
        matriz=datos_sist[:,:,1]
        aux2=dot(matriz,versor_LN)
        aux2=np.arccos(aux2);
        aux=-aux*aux2;
        ang_euler[:,3]=aux;


    elif rotation=='131':
        versor_LN=np.cross(versor_I,datos_sist[:,:,1],axis=1);
        versor_LN=versor_LN/np.linalg.norm(versor_LN,axis=1)[:,None]

        #Calculamos el angulo de rotacion con respecto al eje Z global, corresponde
        #a la primer rotacion, obteniendo el sistema intermedio (xyz)''
        
        aux=dot(versor_J,versor_LN)
        aux=aux/np.linalg.norm(aux,axis=1)[:,None]
        aux=-aux*np.arccos(dot(versor_K,versor_LN))

        ang_euler[:,1]=aux;

        #Calculamos el angulo de rotacion con respecto al eje x'', obteniendo el
        #segundo sistema intermedio (xyz)'

        aux=dot(versor_I,datos_sist[:,:,1]);
        aux=np.arccos(aux);
        ang_euler[:,2]=aux;
        
        #Calculamos el angulo de rotacion con respecto al eje z', obteniendo la
        #rotaci�n a nuestro sistema local xyz.
        ang_euler[:,3]=(dot(datos_sist[:,:,2],versor_LN));
        ang_euler[:,3]=ang_euler[:,3]/np.linalg.norm(ang_euler[:,3],axis=1)[:,None]
        ang_euler[:,3]=ang_euler[:,3]*np.arccos(dot(datos_sist[:,:,3],versor_LN))


    ang_euler=np.real(ang_euler);


    #Corregimos los saltos que puedan darse de -pi a pi 
    ang_euler[:,1]=corr_ang180Original(ang_euler[:,1]);
    ang_euler[:,2]=corr_ang180Original(ang_euler[:,2]);
    ang_euler[:,3]=corr_ang180Original(ang_euler[:,3]);


def dot(a,b):
    return np.sum(a*b,axis=1)

def corr_ang180Original(vector):
    aux=vector;
    ANGULO=vector;
    flag=-1;
    flag2=-1;
    flag3=1;

    for i in range(2,len(vector)):
        
        if ((aux[i-1]*aux[i]<0) and ((np.abs(aux[i]-aux[i-1]))>0.9*np.pi)):
            DIF=mt.ceil(abs(aux[i]-aux[i-1])/(2*np.pi));
            flag=flag*-1;
            flag2=1;
            flag3=-1;
        if (flag==1 and (aux[i-1]*aux[i]<0) ):
            flag2=flag2*-1;
        if (flag==1 and flag2==-1):
            ANGULO[i]=aux[i]-(DIF*2*np.pi*aux[i]/(np.abs(aux[i])));
        if (flag==1 and flag2==1 ):
            ANGULO[i]=aux[i]+(DIF*2*np.pi*aux[i]/(np.abs(aux[i])));    
        
    if (flag3==-1):
        ANGULO=corr_ang180Original(ANGULO);
    
def dif_finitas(v,dt,method="2points"):
    if(method=="2points"):
        """Implementacion de diferencias finitas de dos puntos"""
        dv=np.zeros(v.shape)
        for i in range(1,len(v)-1):
            dv[i,:]=0.5*(v[i+1,:]-v[i-1,:])
        dv[0,:]=v[1,:]-v[0,:]
        dv[-1,:]=v[-1,:]-v[-2,:]
    else:
        """Implementamos la derivada de 4 puntos hacia adelantes y atras"""
        ci=[1/280,-4/105,1/5,-4/5,0,4/5,-1/5,4/105,-1/280];
        dv=np.zeros(v.shape)
        for i in range(4,len(v)-4):
            dv[i,:]=ci[0]*v[i-4,:]+ci[1]*v[i-3,:]+\
                ci[2]*v[i-2,:]+ci[3]*v[i-1,:]+ci[4]*v[i,:]+\
                    ci[5]*v[i+1,:]+ci[6]*v[i+2,:]+ci[7]*v[i+3,:]+\
                        ci[8]*v[i+4,:]

        dv[3,:]=-1/60*v[0,:]+3/20*v[1,:]-3/4*v[2,:]+3/4*v[4,:]-3/20*v[5,:]+1/60*v[6,:]
        dv[2,:]=1/12*v[0,:]-2/3*v[1,:]+2/3*v[3,:]-1/12*v[4,:]
        dv[1,:]=-0.5*v[0,:]+0.5*v[2,:]
        dv[0,:]=-49/20*v[0,:]+6*v[1,:]-15/2*v[2,:]+20/3*v[3,:]-15/4*v[4,:]+6/5*v[5,:]-1/6*v[6,:]

        dv[-4,:]=-1/60*v[-7,:]+3/20*v[-6,:]-3/4*v[-5,:]+3/4*v[-3,:]-3/20*v[-2,:]+1/60*v[-1,:]
        dv[-3,:]=1/12*v[-5,:]-2/3*v[-4,:]+2/3*v[-2,:]-1/12*v[-1,:]
        dv[-2,:]=-0.5*v[-3,:]+0.5*v[-1,:]
        dv[-1,:]=-1/3*v[-4,:]+3/2*v[-3,:]-3*v[-2,:]+6/11*v[-1,:]
    
    return dv/dt    

def Calc_Ang_Art(sist_prox,sist_distal,sist_anat):
    tamanio=sist_prox.shape
    ang_art=np.zeros([tamanio[0],3])

    "Se definen todos los angulos positivos, y pensando en la cadera derecha. "

    "Angulos de flexo/extensión"
    ang_art[:,0]=(dot(sist_anat[:,:,1],sist_prox[:,:,2]))/\
        np.abs(dot(sist_anat[:,:,1],sist_prox[:,:,2]));
    ang_art[:,0]=ang_art[:,0]*np.arccos(dot(sist_anat[:,:,1],sist_prox[:,:,1]));
    "Angulos de abducci�n/aduccion"
    ang_art[:,1]=np.arcsin(dot(sist_prox[:,:,0],sist_distal[:,:,2]));

    "%Angulos de Rot int/ext"
    ang_art[:,2]=dot(sist_anat[:,:,1],sist_distal[:,:,0])/\
        np.abs(dot(sist_anat[:,:,1],sist_distal[:,:,0]))
    ang_art[:,2]=ang_art[:,2]*np.arccos(dot(sist_anat[:,:,1],sist_distal[:,:,1]));


    """ang_art=pd.DataFrame(ang_art).fillna(0, inplace=True)
    ang_art=ang_art.to_numpy()
"""
    return ang_art


def Interp_100(senial):
    n=len(senial)
    m=np.linspace(0,n,n)
    xinterp=np.linspace(0,n,1000)
    f=interp.interp1d(m,senial,kind='cubic',fill_value="extrapolate")
    return f(xinterp)

def Plot3x3AnatPlane(x,y_dic,y_labels,n0,nf,nk):
    fig_ang,ax_ang=plt.subplots(3,3,figsize=(30,30),constrained_layout=True)
    for i in range(3):
        for j in range(3):
            rigth, =ax_ang[j,i].plot(x,(180/np.pi)*y_dic[y_labels[2*i]][n0:nf,j],
            color='#77dd77')
            left, =ax_ang[j,i].plot(x,(180/np.pi)*y_dic[y_labels[2*i+1]][n0:nf,j],
            color='#ff6961')
            if i==0 or i==1:
                if (j==0):
                    ang="Flexion(+)/Extensión(-)"
                elif j==1:
                    ang="Abduction(+)/Adduction(-)"
                else:
                    ang="Internal(+)/External(-) Rotation"
            else:
                if (j==0):
                    ang="Dorsiflexion(+)/ Plantarflexion(-)"
                elif j==1:
                    ang="Inversion(+)/Eversion(-)"
                else:
                    ang="Internal(+)/External(-) Rotation"
            title_ang=y_labels[2*i].replace("_R","")+" "+ang
            ax_ang[j,i].set_title(title_ang)
            ax_ang[j,i].grid()
            ax_ang[j,i].axvspan(x[0],x[nk-n0],color="#e5dde6")
            ax_ang[j,i].axvspan(x[nk-n0],x[-1],color="#fffec2")
            ax_ang[j,i].set_xlabel('Time [s]')
            ax_ang[j,i].set_ylabel('Angle [°]')
            ax_ang[j,i].legend([rigth,left],['Rigth','Left'])

def PlotVectorSegment(x,y_dic,y_labels,n0,nf,nk,title=[],xlabel="",ylabel=""):
    nrows=3
    ncolums=len(y_labels)
    fig_seg,ax_seg=plt.subplots(3,ncolums,figsize=(30,30),constrained_layout=True)
    colors=['#ff6961','#77dd77','#84b6f4']
    for i in range(3):
        for j in range(ncolums):
            ax_seg[i,j].plot(x,y_dic[y_labels[j]][n0:nf,i],color=colors[i])
            ax_seg[i,j].set_title(y_labels[j]+" "+title[i])
            ax_seg[i,j].grid()
            ax_seg[i,j].set_xlabel(xlabel)
            ax_seg[i,j].set_ylabel(ylabel)
            ax_seg[i,j].axvline(x[nk],ymax=np.max(y_dic[y_labels[j]][n0:nf,i]),ymin=np.min(y_dic[y_labels[j]][n0:nf,i]))

        

