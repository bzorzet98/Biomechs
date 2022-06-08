from REGISTER import REGISTER
import numpy as np
import function_biomechs as fb
import matplotlib.pyplot as plt

class BIOMECH_MODEL:
    data={}
    reg={}
    limb="R"
    label_markers=[]
    segments=['Pelvis','Tigh_R','Tigh_L','Leg_R','Leg_L','Foot_R','Foot_L','Foot_R_aux','Foot_L_aux']
    joints=['Hip_R','Hip_L','Knee_R','Knee_L','Ankle_R','Ankle_L']

    def __init__(self,path=''):
        if path =='':
            print('No se ha seleccionado una ruta')
        else:
            self.reg= REGISTER(path)
            if 'PI'in path:
                self.limb='L'
    def __extract_info_register__(self):
        aux={}
        aux['markers']=self.reg.give_data(data='filtered_markers')
        aux['anthropometry']=self.reg.give_data(data='anthropometry')
        aux['events_times']=self.reg.give_data(data='events_times')
        aux['events_frames']=self.reg.give_data(data='events_frames')
        aux['rot_matrix']=self.reg.give_data(data='rot_matrix')
        aux['frequency']=self.reg.give_data(data='frequency')
        aux['frames']=self.reg.frames
        self.reg=aux
        self.label_markers=list(aux['markers'])
    def Return_Data(self, info=''):
        if info == '': 
            print('No existe la información')
        else:
            return self.data[info]
    def Filter_Markers(self, fc, N):
        self.reg.filter_markers(fc, N)
        self.__extract_info_register__()
    def Calc_Virtual_Markers(self): 
        #Ankle Joint
        self.reg['markers']['L_AJC']=\
            (self.reg['markers']['L_MAL_L']+\
            self.reg['markers']['L_MAL_M'])/2
        self.reg['markers']['R_AJC']=\
            (self.reg['markers']['R_MAL_L']+ \
            self.reg['markers']['R_MAL_M'])/2
        #Knee Joint
        self.reg['markers']['L_KJC']=\
            (self.reg['markers']['L_KNEE_EL']+\
            self.reg['markers']['L_KNEE_EM'])/2
        self.reg['markers']['R_KJC']=\
            (self.reg['markers']['R_KNEE_EL']+\
            self.reg['markers']['R_KNEE_EM'])/2
        #Virtual markers between metatarsals
        self.reg['markers']['L_META_I']=\
            (self.reg['markers']['L_META_1']+\
            self.reg['markers']['L_META_5'])/2
        self.reg['markers']['R_META_I']=\
            (self.reg['markers']['R_META_1']+\
            self.reg['markers']['R_META_5'])/2
        #Virtual marker between ASIS
        self.reg['markers']['ASIS_I']=\
            (self.reg['markers']['L_ASIS']+\
            self.reg['markers']['R_ASIS'])/2
    def Calc_Coordinates_Systems(self):

        self.Calc_Virtual_Markers()
        self.data['CorSys']={}
        n=int(self.reg['frames'])
        Sys=np.zeros([n,3,3])

        #Coordinate System of pelvis
        #X
        Sys[:,:,0]=\
            self.reg['markers']['R_ASIS']-self.reg['markers']['L_ASIS']
        Sys[:,:,0]=(Sys[:,:,0])/np.linalg.norm(Sys[:,:,0],axis=1)[:,None]

        #Z
        Sys[:,:,2]=\
            np.cross((self.reg['markers']['R_ASIS']-self.reg['markers']['SACRUM']),\
                (self.reg['markers']['L_ASIS']-self.reg['markers']['SACRUM']),axis=1)
        Sys[:,:,2]=(Sys[:,:,2])/np.linalg.norm(Sys[:,:,2],axis=1) [:,None]

        #Y
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)

        self.data['CorSys']['Pelvis']=Sys
        Sys=np.zeros([n,3,3])
        #Coordinate System of leg
        #Left

        #Y
        Sys[:,:,1]=np.cross((self.reg['markers']['L_MAL_M']-\
            self.reg['markers']['L_KJC']),\
                (self.reg['markers']['L_MAL_L']-\
            self.reg['markers']['L_KJC']),axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]

        #X
        Sys[:,:,0]=self.reg['markers']['L_MAL_M']-\
            self.reg['markers']['L_MAL_L']
        Sys[:,:,0]=Sys[:,:,0]/np.linalg.norm(Sys[:,:,0],axis=1) [:,None]    

        #Z
        Sys[:,:,2]=np.cross(Sys[:,:,0],Sys[:,:,1],axis=1)

        self.data['CorSys']['Leg_L']=Sys
        Sys=np.zeros([n,3,3])
        #Right

        #Y
        Sys[:,:,1]=np.cross((self.reg['markers']['R_MAL_L']-\
            self.reg['markers']['R_KJC']),\
                (self.reg['markers']['R_MAL_M']-\
            self.reg['markers']['R_KJC']),axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]

        #X
        Sys[:,:,0]=self.reg['markers']['R_MAL_L']-\
            self.reg['markers']['R_MAL_M']
        Sys[:,:,0]=Sys[:,:,0]/np.linalg.norm(Sys[:,:,0],axis=1) [:,None]    

        #Z
        Sys[:,:,2]=np.cross(Sys[:,:,0],Sys[:,:,1],axis=1)

        self.data['CorSys']['Leg_R']=Sys
        Sys=np.zeros([n,3,3])
        #Feet Coordinates Systems
        #Rigth
        #Y
        Sys[:,:,1]=self.reg['markers']['R_META_I']-\
            self.reg['markers']['R_HEEL']
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        #X
        Sys[:,:,0]=np.cross(Sys[:,:,1],\
            (self.reg['markers']['R_AJC']-\
                self.reg['markers']['R_HEEL']),axis=1)
        Sys[:,:,0]=Sys[:,:,0]/np.linalg.norm(Sys[:,:,0],axis=1) [:,None]
        #Z
        Sys[:,:,2]=np.cross(Sys[:,:,0],Sys[:,:,1],axis=1)
        self.data['CorSys']['Foot_R_aux']=Sys

        rot_matrix=self.reg['rot_matrix']['R']
        vec=np.zeros([n,3,3])
        for i in range(n):
            vec[i,:,0]=(rot_matrix@Sys[i,:,0])
            vec[i,:,1]=(rot_matrix@Sys[i,:,1])
            vec[i,:,2]=(rot_matrix@Sys[i,:,2])
        Sys=vec

        self.data['CorSys']['Foot_R']=Sys
        Sys=np.zeros([n,3,3])
        #Left
        #Y
        Sys[:,:,1]=self.reg['markers']['L_META_I']-\
            self.reg['markers']['L_HEEL']
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        #X
        Sys[:,:,0]=np.cross(Sys[:,:,1],\
            (self.reg['markers']['L_AJC']-\
                self.reg['markers']['L_HEEL']),axis=1)
        Sys[:,:,0]=Sys[:,:,0]/np.linalg.norm(Sys[:,:,0],axis=1) [:,None]
        #Z
        Sys[:,:,2]=np.cross(Sys[:,:,0],Sys[:,:,1],axis=1)
        self.data['CorSys']['Foot_L_aux']=Sys
        rot_matrix=self.reg['rot_matrix']['L']

        vec=np.zeros([n,3,3])
        for i in range(n):
            vec[i,:,0]=(rot_matrix@Sys[i,:,0])
            vec[i,:,1]=(rot_matrix@Sys[i,:,1])
            vec[i,:,2]=(rot_matrix@Sys[i,:,2])
        Sys=vec

        self.data['CorSys']['Foot_L']=Sys
        Sys=np.zeros([n,3,3])
        #Davis Method for calculus of HJC
        tita=np.radians(28.4)
        beta=np.radians(18)
        d_asis=self.reg['anthropometry']['D_ASIS']/100
        long_mi_R=self.reg['anthropometry']['LONG_MIEM_INF_R']/100
        long_mi_L=self.reg['anthropometry']['LONG_MIEM_INF_L']/100
        asis_sag_R=self.reg['anthropometry']['ASIS_SAGITAL_R']/100
        asis_sag_L=self.reg['anthropometry']['ASIS_SAGITAL_L']/100

        C_R=0.115*long_mi_R-0.0153
        C_L=0.115*long_mi_L-0.0153
            
        yh_R=(-asis_sag_R)*np.cos(beta)+C_R*np.cos(tita)*np.sin(beta)
        yh_L=(-asis_sag_L)*np.cos(beta)+C_L*np.cos(tita)*np.sin(beta)
        xh_R=(d_asis/2)-C_R*np.sin(tita)
        xh_L=(d_asis/2)-C_L*np.sin(tita)
        zh_R=(-asis_sag_R)*np.sin(beta)-C_R*np.cos(tita)*np.cos(beta)
        zh_L=(-asis_sag_L)*np.sin(beta)-C_L*np.cos(tita)*np.cos(beta)

        #Hip Joint
        self.reg['markers']['R_HJC']=\
            self.reg['markers']['ASIS_I']+\
                yh_R*self.data['CorSys']['Pelvis'][:,:,1]+\
                    xh_R*self.data['CorSys']['Pelvis'][:,:,0]+\
                        zh_R*self.data['CorSys']['Pelvis'][:,:,2]
        self.reg['markers']['L_HJC']=\
            self.reg['markers']['ASIS_I']+\
                yh_L*self.data['CorSys']['Pelvis'][:,:,1]-\
                    xh_L*self.data['CorSys']['Pelvis'][:,:,0]+\
                        zh_L*self.data['CorSys']['Pelvis'][:,:,2]

        #Coordinate System of tigh
        #Rigth
        #Z
        Sys[:,:,2]=\
            self.reg['markers']['R_HJC']-\
                self.reg['markers']['R_KJC']
        Sys[:,:,2]=(Sys[:,:,2])/np.linalg.norm(Sys[:,:,2],axis=1)[:,None]   

        #Y
        Sys[:,:,1]=\
            np.cross((self.reg['markers']['R_KNEE_EL']-\
                self.reg['markers']['R_HJC']),\
                (self.reg['markers']['R_KNEE_EM']\
                    -self.reg['markers']['R_HJC']),axis=1)
        Sys[:,:,1]=(Sys[:,:,1])/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]

        #X
        Sys[:,:,0]=np.cross(Sys[:,:,1],Sys[:,:,2],axis=1)
        self.data['CorSys']['Tigh_R']=Sys

        Sys=np.zeros([n,3,3])
        #Left
        #Z
        Sys[:,:,2]=\
            self.reg['markers']['L_HJC']-\
                self.reg['markers']['L_KJC']
        Sys[:,:,2]=(Sys[:,:,2])/np.linalg.norm(Sys[:,:,2],axis=1)[:,None]   

        #Y
        Sys[:,:,1]=\
            np.cross((self.reg['markers']['L_KNEE_EM']-\
                self.reg['markers']['L_HJC']),\
                (self.reg['markers']['L_KNEE_EL']\
                    -self.reg['markers']['L_HJC']),axis=1)
        Sys[:,:,1]=(Sys[:,:,1])/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]

        #X
        Sys[:,:,0]=np.cross(Sys[:,:,1],Sys[:,:,2],axis=1)
        self.data['CorSys']['Tigh_L']=Sys
        Sys=np.zeros([n,3,3])
    def Calc_Anatomic_System (self):
        n=self.reg['frames']
        Sys=np.zeros((n[0],3,3))
        self.data['SysAnat']={}
        #HIP Joint Coordinate System
        #Righ
        Sys[:,:,0]=self.data['CorSys']['Pelvis'][:,:,0]
        Sys[:,:,2]=self.data['CorSys']['Tigh_R'][:,:,2]
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        self.data['SysAnat']['HJC_R']=Sys
        #Left
        Sys=np.zeros((n[0],3,3))
        Sys[:,:,0]=self.data['CorSys']['Pelvis'][:,:,0]
        Sys[:,:,2]=self.data['CorSys']['Tigh_L'][:,:,2]
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        self.data['SysAnat']['HJC_L']=Sys
        #KNEE Joint Coordinate System
        Sys=np.zeros((n[0],3,3))
        #Righ
        Sys[:,:,0]=self.data['CorSys']['Tigh_R'][:,:,0]
        Sys[:,:,2]=self.data['CorSys']['Leg_R'][:,:,2]
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        self.data['SysAnat']['KJC_R']=Sys
        #Left
        Sys=np.zeros((n[0],3,3))
        Sys[:,:,0]=self.data['CorSys']['Tigh_L'][:,:,0]
        Sys[:,:,2]=self.data['CorSys']['Leg_L'][:,:,2]
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        self.data['SysAnat']['KJC_L']=Sys
        #ANKLE Joint Coordinate System
        Sys=np.zeros((n[0],3,3))
        #Righ
        Sys[:,:,0]=self.data['CorSys']['Leg_R'][:,:,0]
        Sys[:,:,2]=self.data['CorSys']['Foot_R'][:,:,2]
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        self.data['SysAnat']['AJC_R']=Sys
        #Left
        Sys=np.zeros((n[0],3,3))
        Sys[:,:,0]=self.data['CorSys']['Leg_L'][:,:,0]
        Sys[:,:,2]=self.data['CorSys']['Foot_L'][:,:,2]
        Sys[:,:,1]=np.cross(Sys[:,:,2],Sys[:,:,0],axis=1)
        Sys[:,:,1]=Sys[:,:,1]/np.linalg.norm(Sys[:,:,1],axis=1) [:,None]
        self.data['SysAnat']['AJC_L']=Sys
    def Joint_Angle(self):
        self.data['JointAng']={}
        #HIP
        self.data['JointAng']['Hip_R']=fb.Calc_Ang_Art(self.data['CorSys']['Pelvis'],\
            self.data['CorSys']['Tigh_R'],self.data['SysAnat']['HJC_R'])
        self.data['JointAng']['Hip_R'][:,1]=-self.data['JointAng']['Hip_R'][:,1]    
        self.data['JointAng']['Hip_L']=fb.Calc_Ang_Art(self.data['CorSys']['Pelvis'],\
            self.data['CorSys']['Tigh_L'],self.data['SysAnat']['HJC_L'])
        self.data['JointAng']['Hip_L'][:,2]=-self.data['JointAng']['Hip_L'][:,2]    
        #KNEE
        self.data['JointAng']['Knee_R']=fb.Calc_Ang_Art(self.data['CorSys']['Tigh_R'],\
            self.data['CorSys']['Leg_R'],self.data['SysAnat']['KJC_R'])
        self.data['JointAng']['Knee_R'][:,0]=-self.data['JointAng']['Knee_R'][:,0]
        self.data['JointAng']['Knee_R'][:,1]=-self.data['JointAng']['Knee_R'][:,1]
        self.data['JointAng']['Knee_L']=fb.Calc_Ang_Art(self.data['CorSys']['Tigh_L'],\
            self.data['CorSys']['Leg_L'],self.data['SysAnat']['KJC_L'])
        self.data['JointAng']['Knee_L'][:,0]=-self.data['JointAng']['Knee_L'][:,0]
        self.data['JointAng']['Knee_L'][:,2]=-self.data['JointAng']['Knee_L'][:,2]
        #Ankle
        self.data['JointAng']['Ankle_R']=fb.Calc_Ang_Art(self.data['CorSys']['Leg_R'],\
            self.data['CorSys']['Foot_R'],self.data['SysAnat']['AJC_R'])
        self.data['JointAng']['Ankle_R'][:,1]=-self.data['JointAng']['Ankle_R'][:,1]
        self.data['JointAng']['Ankle_L']=fb.Calc_Ang_Art(self.data['CorSys']['Leg_L'],\
            self.data['CorSys']['Foot_L'],self.data['SysAnat']['AJC_L'])
        self.data['JointAng']['Ankle_L'][:,2]=-self.data['JointAng']['Ankle_L'][:,2]
    def Joint_Angle_Vel(self,method='2points'):
        self.data['JointAngVel']={}
        aux={}
        dt=1/self.reg['frequency']
        for joint in self.joints:
            aux[joint]=fb.dif_finitas(self.data['JointAng'][joint],dt,method=method)
        if method=="2points":
            self.data['JointAngVel']=aux
        else:
            self.data['JointAngVel_4p']=aux
    def Ang_Euler(self, sys, rot,sys_o=[]):
        """Con esta funcion podemos elegir si se desea calcular 
        los angulos de Euler y la razon de cambio de estos angulos, 
        en base al sistema global o en otro sistema inercial"""
        aux={}
        der={}
        dt=1/self.reg['frequency']
        if len(sys_o)==0:
            for i in range(len(sys)):
                aux[sys[i]]=fb.Euler_Angle(self.data['CorSys'][sys[i]], rot[i])
                der[sys[i]]=fb.dif_finitas(aux[sys[i]],dt)
            self.data['AngEuler']=aux
            self.data['DerAngEuler']=der
        else: 
            for i in range(len(sys)):
                aux[sys[i]]=fb.Euler_Angle(self.data[sys[i]], rot[i], self.data[sys_o[i]])
                der[sys[i]]=fb.dif_finitas(aux[sys[i]],dt)
            self.data['AngEulerRelative']=aux
            self.data['DerAngEulerRelative']=der

    def Plot_Markers(self,ind=0):
        figure = plt.figure(figsize=[100,100])
        ax = figure.add_subplot(projection = '3d')

        for i in range(len(self.label_markers)):
            ax.scatter(self.reg['markers'][self.label_markers[i]][ind,0],
                            self.reg['markers'][self.label_markers[i]][ind,1],
                            self.reg['markers'][self.label_markers[i]][ind,2],
                            marker='o',color='black',
                            label = self.label_markers[i])
            ax.text(self.reg['markers'][self.label_markers[i]][ind,0],
                            self.reg['markers'][self.label_markers[i]][ind,1],
                            self.reg['markers'][self.label_markers[i]][ind,2],self.label_markers[i], color='black',fontsize='xx-small')
        ax.quiver(0,0,0,1,0,0,color='black',length=1)
        ax.quiver(0,0,0,0,1,0,color='black',length=1)
        ax.quiver(0,0,0,0,0,1,color='black',length=1)
        ax.set_title('Puntos ')
        ax.set_adjustable("box")
        ax.set(xlim3d=(0, 2), xlabel='X')
        ax.set(ylim3d=(-1, 1), ylabel='Y')
        ax.set(zlim3d=(0, 2), zlabel='Z')
        
        
        #ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,fontsize='small',ncol=3)
        plt.show()
    def Plot_Sist_Coor(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        system=['Pelvis','Tigh_R','Tigh_L','Leg_R','Leg_L','Foot_R','Foot_L','Foot_R_aux','Foot_L_aux']
        orig=['ASIS_I','R_HJC','L_HJC','R_KJC','L_KJC','R_AJC','L_AJC','R_META_I','L_META_I']
        i=0
        for j in range(len(system)):
            sys=self.data['CorSys'][system[j]]
            ori=self.reg['markers'][orig[j]]
            ax.quiver(ori[i,0],ori[i,1],ori[i,2],sys[i,0,0],sys[i,1,0],sys[i,2,0],color='red',length=0.1)
            ax.quiver(ori[i,0],ori[i,1],ori[i,2],sys[i,0,1],sys[i,1,1],sys[i,2,1],color='green',length=0.1)
            ax.quiver(ori[i,0],ori[i,1],ori[i,2],sys[i,0,2],sys[i,1,2],sys[i,2,2],color='blue',length=0.1)
            ax.scatter(ori[i,0],ori[i,1],ori[i,2],color='black')
        ax.quiver(0,0,0,1,0,0,color='black',length=1)
        ax.quiver(0,0,0,0,1,0,color='black',length=1)
        ax.quiver(0,0,0,0,0,1,color='black',length=1)
        ax.set_title('Sistemas Coordenados')
        ax.set(xlim3d=(0, 2), xlabel='X')
        ax.set(ylim3d=(-1, 1), ylabel='Y')
        ax.set(zlim3d=(0, 2), zlabel='Z')
        ax.legend(['x','y','z'])
        plt.show()
    def Plot_Joints_Angles(self,norm=None):
        if norm==None:
            t0=self.reg['events_times']['Inicio']
            tf=self.reg['events_times']['Fin']
            n0=self.reg['events_frames']['Inicio']
            nf=self.reg['events_frames']['Fin']
            nk=self.reg['events_frames']['Golpe']
            t=np.linspace(t0,tf,nf-n0)
            fb.Plot3x3AnatPlane(t,self.data['JointAng'],self.joints,n0,nf,nk)
    def Plot_Joints_Angles_Velocities(self,norm=None,method="2points"):
        if norm==None:
            t0=self.reg['events_times']['Inicio']
            tf=self.reg['events_times']['Fin']
            n0=self.reg['events_frames']['Inicio']
            nf=self.reg['events_frames']['Fin']
            nk=self.reg['events_frames']['Golpe']
            t=np.linspace(t0,tf,nf-n0)
            if method=="2points":
                fb.Plot3x3AnatPlane(t,self.data['JointAngVel'],self.joints,n0,nf,nk)
            else:
                fb.Plot3x3AnatPlane(t,self.data['JointAngVel_4p'],self.joints,n0,nf,nk)
    def Plot_Euler_Angles(self,norm=None):
        if norm==None:
            t0=self.reg['events_times']['Inicio']
            tf=self.reg['events_times']['Fin']
            n0=self.reg['events_frames']['Inicio']
            nf=self.reg['events_frames']['Fin']
            nk=self.reg['events_frames']['Golpe']
            t=np.linspace(t0,tf,nf-n0)
            title=['Alfa', 'Beta', 'Gamma']
            fb.PlotVectorSegment(t,self.data['AngEuler'],self.segments,n0,nf,nk,title=title,xlabel="Tiempo [s]",ylabel="Angle [rads]") 

