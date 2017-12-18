# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:48:19 2017

@author: YDH
"""
#this is base line
import matplotlib.pylab as plt
import tkinter as tk
from tkinter import ttk
from scipy import interpolate
import time as tm
from scipy.signal import fftconvolve 
import numpy as np

from scipy.constants.constants import c,pi,kilo,milli,micro,nano,mega

LARGE_FONT=("맑은고딕", 12)

class MMC:
    def __init__(self):
        self.parent = tk.Tk()
        self.parent.title("Single Lens Mode-matching Calculator")
        
        menubar = tk.Menu(self.parent)
        plotmenu = tk.Menu(menubar, tearoff=0)
        plotmenu.add_command(label="Full-scope Beamsize", command=self.vis_full_w)
        plotmenu.add_command(label="Cavity-scope Beamsize", command=self.vis_cav_w)
        plotmenu.add_separator()
        plotmenu.add_command(label="Full-scope Wavefront", command=self.vis_full_R)
        plotmenu.add_command(label="Cavity-scope Wavefront", command=self.vis_cav_R)
        plotmenu.add_separator()
        plotmenu.add_command(label="Outplut Intensity", command=self.vis_CE)
        plotmenu.add_separator()
        plotmenu.add_command(label="Cavity Modes", command=self.vis_transmission)        
        
        menubar.add_cascade(label="Plot", menu=plotmenu)
        
        self.leftf = tk.Frame(self.parent)
        self.leftf.pack(side=tk.LEFT)
        self.rightf = tk.Frame(self.parent)
        self.rightf.pack(side=tk.RIGHT)
        
        "여기서부터 왼쪽 프레임"
        
        introf=tk.Frame(self.leftf)
        tk.Label(introf, text="Single Lens Mode-matching Calculator", font=LARGE_FONT).pack()
        tk.Label(introf, text="Developed by YDH").pack(side=tk.RIGHT,pady=10)
        introf.pack() 
        
        self.mainf=tk.Frame(self.leftf)
        self.mainf.pack()
        
        tk.Label(self.mainf, text="ROC of Input Mirror, R1 :").grid(row=0, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="Reflectivity of Input Mirror :").grid(row=1, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="ROC of Output Mirror, R2 :").grid(row=2, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="Reflectivity of Output Mirror :").grid(row=3, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="Length of Cavity, L :").grid(row=4, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="frquency of PZT :").grid(row=5, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="Incidenting wavelength of Laser :").grid(row=6, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="Incidenting FWHM of Laser :").grid(row=7, pady=3, sticky=tk.E)
        tk.Label(self.mainf, text="Focal length of Modematching lens :").grid(row=8, pady=3,sticky=tk.E)
        tk.Label(self.mainf, text="Beam waist size of Starting point :").grid(row=9, pady=3, sticky=tk.E)
        
        
        tk.Label(self.mainf, text="mm").grid(row=0, column=2, sticky=tk.E)
        tk.Label(self.mainf, text="%").grid(row=1, column=2, sticky=tk.E)        
        tk.Label(self.mainf, text="mm").grid(row=2, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="%").grid(row=3, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="mm").grid(row=4, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="Hz").grid(row=5, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="nm").grid(row=6, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="kHz").grid(row=7, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="mm").grid(row=8, column=2,sticky=tk.E)
        tk.Label(self.mainf, text="μm").grid(row=9, column=2,sticky=tk.E)
               
        self.e1 = ttk.Entry(self.mainf)
        self.e2 = ttk.Entry(self.mainf)
        self.e3 = ttk.Entry(self.mainf)
        self.e4 = ttk.Entry(self.mainf)
        self.e5 = ttk.Entry(self.mainf)
        self.e6 = ttk.Entry(self.mainf)
        self.e7 = ttk.Entry(self.mainf)
        self.e8 = ttk.Entry(self.mainf)
        self.e9 = ttk.Entry(self.mainf)
        self.e10 = ttk.Entry(self.mainf)
        
        
        self.e1.grid(row=0, column=1,sticky=tk.W)
        self.e2.grid(row=1, column=1,sticky=tk.W)
        self.e3.grid(row=2, column=1,sticky=tk.W)
        self.e4.grid(row=3, column=1,sticky=tk.W)        
        self.e5.grid(row=4, column=1,sticky=tk.W)        
        self.e6.grid(row=5, column=1,sticky=tk.W)
        self.e7.grid(row=6, column=1,sticky=tk.W)
        self.e8.grid(row=7, column=1,sticky=tk.W)
        self.e9.grid(row=8, column=1,sticky=tk.W)
        self.e10.grid(row=9, column=1,sticky=tk.W)
        
        
        f=tk.Frame(self.mainf) #계산하기 버튼 프레임
        ttk.Button(f, text='Calculate', command=self.setPARA).pack(side=tk.LEFT) 
        
        f.grid(row=10, column=1)
        
        self.statusf=tk.Frame(self.leftf) #상태메시지 프레임
        self.statusf.pack()
        
        self.status=tk.Label(self.statusf, text=" ")
        self.status.grid(row=1)
        self.status2=tk.Label(self.statusf, text="  ")
        self.status2.grid(row=2)
        self.status3=tk.Label(self.statusf, text="  ")
        self.status3.grid(row=3)
        
        f2=tk.Frame(self.leftf) #결과 프레임
        tk.Label(f2, text="Beam waist size in resonator : ").grid(row=0, sticky=tk.E)
        tk.Label(f2, text="from Input mirror to Beam waist distance : ").grid(row=1, sticky=tk.E)
        tk.Label(f2, text="from Laser to lens distance : ").grid(row=2, sticky=tk.E)
        tk.Label(f2, text="from lens to Input mirror distance : ").grid(row=3, sticky=tk.E)
        tk.Label(f2, text="characteristic focal length of lens : ").grid(row=4, sticky=tk.E)
        tk.Label(f2, text="                      ").grid(row=5, sticky=tk.E)
        tk.Label(f2, text="Finesse : ").grid(row=6, sticky=tk.E)
        tk.Label(f2, text="Free Spectral Range(frequency) : ").grid(row=7, sticky=tk.E)
        tk.Label(f2, text="Free Spectral Range(cavity length) : ").grid(row=8, sticky=tk.E)
        tk.Label(f2, text="Full width at half maximum(frequency) : ").grid(row=9, sticky=tk.E)
        tk.Label(f2, text="Full width at half maximum(cavity length) : ").grid(row=10, sticky=tk.E)
        
        self.Result1=tk.Label(f2, text="                                       ")
        self.Result2=tk.Label(f2, text="                                       ")
        self.Result3=tk.Label(f2, text="                                       ")
        self.Result4=tk.Label(f2, text="                                       ")
        self.Result5=tk.Label(f2, text="                                       ")
        self.Result6=tk.Label(f2, text="                                       ")
        self.Result7=tk.Label(f2, text="                                       ")
        self.Result8=tk.Label(f2, text="                                       ")
        self.Result9=tk.Label(f2, text="                                       ")
        self.Result10=tk.Label(f2, text="                                       ")
        self.Result11=tk.Label(f2, text="                                       ")
        
        self.Result1.grid(row=0,column=1,sticky=tk.W)
        self.Result2.grid(row=1,column=1,sticky=tk.W)
        self.Result3.grid(row=2,column=1,sticky=tk.W)
        self.Result4.grid(row=3,column=1,sticky=tk.W)
        self.Result5.grid(row=4,column=1,sticky=tk.W)
        self.Result6.grid(row=5,column=1,sticky=tk.W)
        self.Result7.grid(row=6,column=1,sticky=tk.W)
        self.Result8.grid(row=7,column=1,sticky=tk.W)
        self.Result9.grid(row=8,column=1,sticky=tk.W)
        self.Result10.grid(row=9,column=1,sticky=tk.W)
        self.Result11.grid(row=10,column=1,sticky=tk.W)
        
        
        f2.pack()
        
        "여기까지 왼쪽 프레임 밑으로는 오른쪽 프레임"
        
        self.parent.config(menu=menubar)
        self.parent.mainloop()
        
    def Result(self,n):
        self.Result1.configure(text="%s μm                        "%round(self.w2/micro,2))
        self.Result2.configure(text="%s mm"%round(self.t/milli,2))
        self.Result3.configure(text="%s mm"%round(self.d1/milli,2))
        self.Result4.configure(text="%s mm"%round(self.d2/milli,2))
        self.Result5.configure(text="%s mm"%round(self.f0/milli,2))
        
        self.Result7.configure(text="%s "%self.Finesse)
        self.Result8.configure(text="%s MHz"%round(self.f_FSR/mega,2))
        self.Result9.configure(text="%s nm"%round(self.l_FSR/nano,2))
        self.Result10.configure(text="%s kHz"%round(self.f_FWHM/kilo,2))
        self.Result11.configure(text="%s nm"%round(self.l_FWHM/nano,3))
        
        
    def setPARA(self):
        self.R1 = float(self.e1.get())*milli
        self.Ref1 = float(self.e2.get())/100
        self.r1 = np.sqrt(self.Ref1)
        self.Ref2 = float(self.e4.get())/100
        self.r2 = np.sqrt(self.Ref2)
        self.Refl = np.sqrt(self.Ref1*self.Ref2)
        self.R2 = float(self.e3.get())*milli
        self.pzt_freq = float(self.e6.get())
        self.L=float(self.e5.get())*milli
        self.I0= 1
        
        if self.stability(self.R1,self.R2,self.L)==1:
            self.redmessage1("This resonator is unstabled.")
        
        self.wl=float(self.e7.get())*nano
        self.w_L=float(self.e8.get())*kilo
        self.fl=float(self.e9.get())*milli
        self.w1=float(self.e10.get())*micro
             
        Reson=Resonator(self.R1,self.R2,self.L,self.wl)
        param1=Reson.get()
        self.w2=param1[3]
        self.t=param1[4]
        
        Mode=Modematching(self.w1,self.w2,self.fl,self.wl)
        param2=Mode.get()
        self.f0=param2[4]
        if self.f0 > self.fl:
            self.redmessage3("focal length of lens must bigger than %s mm." % (round(self.f0/milli,2)))
        else:
            self.bluemessage3("This lens is able to use to mode-match.")
        self.d1=param2[5]
        self.d2=param2[6]
        
        self.Finesse = int(np.pi*np.sqrt(abs(self.r1*self.r2))/(1-abs(self.r1*self.r2)))
        self.f_FSR = int(c/(2*self.L))
        self.l_FSR = float(self.wl/2)
        self.f_FWHM = int(self.f_FSR/self.Finesse)
        self.l_FWHM = float(self.l_FSR/self.Finesse)
        
        self.Result(1)
        
    def stability(self,R1,R2,L):
        g1=round(1-L/R1,2)
        g2=round(1-L/R2,2)
        stable=round(g1*g2,2)
        print(" ")
        if stable >= 1 or stable <= 0:
            self.redmessage2("g1: %s , g2: %s, g1g2: %s" % (g1,g2,stable)) 

            return 1
        self.bluemessage1("  This Resonator is stabled.  ")
        self.bluemessage2("g1: %s , g2: %s, g1g2: %s" % (g1,g2,stable))        
        return 0
        
    def bluemessage1(self,m):
        self.status.configure(text=m, fg='blue')
        
    def bluemessage2(self,m):
        self.status2.configure(text=m, fg='blue')
        
    def bluemessage3(self,m):
        self.status3.configure(text=m, fg='blue')
        
    def redmessage1(self,m):
        self.status.configure(text=m, fg='red')
        
    def redmessage2(self,m):
        self.status2.configure(text=m, fg='red')
        
    def redmessage3(self,m):
        self.status3.configure(text=m, fg='red')    
        
    def vis_modematching(self):
        a=Listizing(self.R1,self.R2,self.L,self.wl, self.fl,self.w1, self.w2,self.t,self.d1,self.d2)
        
        wlist=np.array(a.get_cavity_w())
        wlist_rev=-1*self.wlist

        cavR = a.get_cavity_R()
        fullw = a.get_full_w()
        fullR = a.get_full_R()
        
        #plotting graphs
        fig = plt.figure(figsize=(15,10))
        fig.subplots_adjust(hspace=0.8)
    
        data1 = fig.add_subplot(3,2,1)
        plt.figure('Resonator')
        plt.title('Beam Contour in Resonator')
        plt.xlabel('distance(mm)')
        plt.ylabel('beam size(microm)')
        plt.plot(self.wlist, '-r')
        plt.plot(self.wlist_rev, '-r')
        plt.axhline(0, color='yellow')
        
        plt.xlabel('Time(µs)')
        plt.ylabel('Amplitude')
        plt.title('Signal of monochromatic wave')
        data1.plot(time,result, color='red', label='Monochromatic wave')
        plt.legend()
        plt.minorticks_on()
        
        data2 = fig.add_subplot(3,2,2)
        plt.xlabel('Time(µs)')
        plt.ylabel('Amplitude')
        plt.title('Signal with %s kHz Lineshape' %(self.e8.get()))
        data2.plot(time,result2, color='blue', linestyle='solid', label='with Lineshape', linewidth=1.5)
        plt.legend()
        plt.minorticks_on()
        
        data3 = fig.add_subplot(3,1,2)
        plt.xlabel('Time(µs)')
        plt.ylabel('Amplitude')
        plt.title('Build up, Ringdown Signal')
        data2.plot(time,result3, color='green', linestyle='solid', label='ringdown signal', linewidth=1.5)
        plt.legend()
        plt.minorticks_on()
        
        data4 = fig.add_subplot(3,1,3)


        plt.show()
        
    def vis_cav_R(self):
        a=Listizing(self.R1,self.R2,self.L,self.wl, self.fl,self.w1, self.w2,self.t,self.d1,self.d2)
        Rlist=a.get_cavity_R()


        plt.figure('Resonator Wavefront')
        plt.title('Beam wavefront in Resonator')
        plt.xlabel('distance(mm)')
        plt.ylabel('size(mm)')
        plt.axis([0,len(Rlist),-3*self.R1*kilo,3*self.R1*kilo])
        plt.plot(Rlist)   
        plt.show()
        
    def vis_cav_w(self):
        a=Listizing(self.R1,self.R2,self.L,self.wl, self.fl,self.w1, self.w2,self.t,self.d1,self.d2)
        self.wlist=np.array(a.get_cavity_w())
        self.wlist_rev=-1*self.wlist

        plt.figure('Resonator')
        plt.title('Beam Contour in Resonator')
        plt.xlabel('distance(mm)')
        plt.ylabel('beam size(microm)')
        plt.plot(self.wlist, '-r')
        plt.plot(self.wlist_rev, '-r')
        plt.axhline(0, color='yellow')
        plt.show()
        
    def vis_full_w(self):
        a=Listizing(self.R1,self.R2,self.L,self.wl, self.fl,self.w1, self.w2,self.t,self.d1,self.d2)
        wlist=np.array(a.get_full_w()) #플롯용 전체 w 리스트를 어레이화 함
        wlist_rev=-1*wlist # -쪽 플롯
        
        cwlist=a.get_cavity_w()
        ccwlist=np.array([np.nan]*(len(wlist)-len(cwlist))+cwlist)
        ccwlist_rev=-1*ccwlist
        
        plt.figure('Beam Contour_Normal')
        plt.title('Beam Contour')
        plt.xlabel('distance(mm)')
        plt.ylabel('beam size(microm)')
        plt.plot(wlist, '-b', label = 'Beam Propagation')
        plt.plot(wlist_rev, '-b')
        plt.plot(ccwlist,'-r',linewidth=1.5, label = 'Cavity')
        plt.plot(ccwlist_rev,'-r',linewidth=1.5)
        plt.axhline(0, color='yellow')
        plt.legend()
        plt.show()
        
    def vis_full_R(self):
        a=Listizing(self.R1,self.R2,self.L,self.wl, self.fl,self.w1, self.w2,self.t,self.d1,self.d2)
        Rlist=a.get_full_R()
        cRlist=a.get_cavity_R()
        ccRlist=[None]*(len(Rlist)-len(cRlist))+cRlist
        
        plt.figure('Beam Wavefront_fullscope')
        plt.title('Beam Wavefront-Full range')
        plt.xlabel('distance(mm)')
        plt.ylabel('size(mm)')
        plt.axis([0,len(Rlist), -5*self.R1*kilo, 5*self.R1*kilo])
        plt.plot(Rlist, '-b', label = 'Wavefront Propagation')
        plt.plot(ccRlist, '-r',linewidth=1.5, label = 'Cavity')
        plt.legend()
        plt.show()
        
    def vis_CE(self):
        if self.Refl > 0.99999:
            self.redmessage3("Above 99.999%, the reflectivity is too large to be calculated.")
        else:
            a=CouplingEfficiency(self.Refl,self.wl,self.L,self.pzt_freq,self.I0,self.w_L)
            time = a.get_axis_time()
            result = a.get_mono()
            result2 = a.get_lineshape()
            result3 = a.get_ringdown()
            
            #plotting graphs
            fig = plt.figure()
            fig.subplots_adjust(hspace=0.8)
        
            data1 = fig.add_subplot(3,1,1)
            plt.xlabel('Time(µs)')
            plt.ylabel('Amplitude')
            plt.title('Signal of monochromatic wave')
            data1.plot(time,result, color='red', label='Monochromatic wave')
            plt.legend()
            plt.minorticks_on()
            data2 = fig.add_subplot(3,1,2)
            plt.xlabel('Time(µs)')
            plt.ylabel('Amplitude')
            plt.title('Signal with %s kHz Lineshape' %(self.e8.get()))
            data2.plot(time,result2, color='blue', linestyle='solid', label='with Lineshape', linewidth=1.5)
            plt.legend()
            plt.minorticks_on()
            data2 = fig.add_subplot(3,1,3)
            plt.xlabel('Time(µs)')
            plt.ylabel('Amplitude')
            plt.title('Build up, Ringdown Signal')
            data2.plot(time,result3, color='green', linestyle='solid', label='ringdown signal', linewidth=1.5)
            plt.legend()
            plt.minorticks_on()
    
    
            plt.show()
        
    def vis_transmission(self):
        a=transmission(self.Ref1, self.Ref2, self.L)
        Transmission_xplot = a.get_Transmission_xplot()
        Transmission_y = a.get_Transmission_y()
        plt.figure('Transmission Peak')
        plt.title('Transmission Peaks')
        plt.xlabel('frequency (MHz)')
        plt.ylabel('Transmission rate')
        
        plt.plot(Transmission_xplot,Transmission_y, label='Transmission peak')
        plt.minorticks_on()
        plt.legend

        plt.show()
        
        
class Resonator:
    def __init__(self, R1, R2, L, wl):
        self.R1=R1
        self.R2=R2
        self.L=L
        self.wl=wl
        self.ABCD()
        
    def ABCD(self):
        #"1.resonance_w2 함수를 이용하여 공진조건에서의 cavity내 beam waist size를 구한다."   
        ABCD=self.ABCD_lin(self.R1,self.R2,self.L)
        A=ABCD[0,0]
        B=ABCD[0,1]
        C=ABCD[1,0]
        D=ABCD[1,1]
        self.w2=self.resonance_w2(A,B,C,D,self.wl,self.R1) #Modeling 시 R1 바로 앞의 조건으로 계산했으므로 R1만기입
        self.t=abs((A-D)/(2*C))
        
    #linear Cavity의 lens wave guide 로 해석한 ABCD MATRIX
    def ABCD_lin(self,R1,R2,L):
        return self.lens(R1/2)*self.free(L)*self.lens(R2/2)*self.free(L)
    
    #Free space propagation 계산기
    def free(self,d):
        return np.matrix([[1,d],[0,1]])
    
    #thin lens 계산기
    def lens(self,fl):
        return np.matrix([[1,0],[-1/fl,1]])
    
        
    """★ Cavity 공진 Self Consistency 조건의  구하는 함수 ★"""
    def resonance_w2(self,A,B,C,D,wl,R):
        w2=np.sqrt(abs((wl*R*(A-D)/(4*np.pi*B*C))*np.sqrt(abs(4-(A+D)**2)))) #(58)식과 (23)식을 조합해서 만든 함수
        return w2
        
    def get(self):
        return [self.R1,self.R2,self.L,self.w2,self.t,self.wl]
    
    def show(self):
        print(self.R1,self.R2,self.L,self.w2,self.t,self.wl)
        
        
class Modematching:
    def __init__(self,w1, w2, fl, wl):
        self.w1=w1
        self.w2=w2
        self.fl=fl
        self.wl=wl
        
        self.f0=np.pi*self.w1*self.w2/self.wl  # 특성 focal length.
        self.d1=self.fl+self.w1/self.w2*np.sqrt(self.fl**2-self.f0**2)  # 초점1과 렌즈 사이거리
        self.d2=self.fl+self.w2/self.w1*np.sqrt(self.fl**2-self.f0**2) 
        
    def get(self):
        return [self.w1, self.w2, self.fl, self.wl, self.f0, self.d1, self.d2]
    
    
class Listizing:
    
    def __init__(self,R1,R2,L,wl,fl,w1,w2,t,d1,d2):
        self.R1=R1
        self.R2=R2
        self.L=L
        self.wl=wl
        self.fl=fl
        self.w1=w1
        self.w2=w2
        self.t=t
        self.d1=d1
        self.d2=d2
        self.fullscope_R()
        self.fullscope_w()
        self.Cavity_R()
        self.Cavity_w()
        
    def fullscope_w(self):
        distance1=np.arange(0,self.d1,milli)
        data1=[]
        for x in distance1:
            y=self.ABCD_d1(x,self.wl,self.w1)*mega
            data1.append(y)
    
        distance2=np.arange(self.d1,self.d1+self.d2+self.L-self.t,milli)
        data2=[]
        for x in distance2:
            y=self.ABCD_d2(self.d1,self.fl,x,self.wl,self.w1)*mega
            data2.append(y)
        
        self.fullwdata=data1+data2

    def fullscope_R(self):
        distance1=np.arange(0,self.d1,milli)
        data1=[]
        for x in distance1:
            y=self.ABCD_d1_R(x,self.wl,self.w1)*kilo
            data1.append(y)
        
        distance2=np.arange(self.d1,self.d1+self.d2+self.L-self.t,milli)
        data2=[]
        for x in distance2:
            y=self.ABCD_d2_R(self.d1,self.fl,x,self.wl,self.w1)*kilo
            data2.append(y)
        
        self.fullRdata=data1+data2
        

    def Cavity_w(self): #곡률반지름1,곡률반지름2,캐비티길이,waist size,파장
        distance1=np.arange(-self.t,self.L-self.t,milli)
        data1=[]
        for x in distance1:
            y=self.ABCD_d1(x,self.wl,self.w2)*mega
            data1.append(y)
        
        self.cavitywdata=data1

    def Cavity_R(self):
        distance1=np.arange(-self.t,self.L-self.t,milli)
        data1=[]
        for x in distance1:
            y=self.ABCD_d1_R(x,self.wl,self.w2)/milli
            data1.append(y)
            
        self.cavityRdata=data1
        
    def get_cavity_R(self):
        return self.cavityRdata
                  
    def get_cavity_w(self):
        return self.cavitywdata
    
    def get_full_R(self):
        return self.fullRdata
    
    def get_full_w(self):
        return self.fullwdata

    "출발점부터 렌즈까지"   
    def ABCD_d1(self,d1,wl,w1):
        ABCD=self.free(d1) 
        A=ABCD[0,0]
        B=ABCD[0,1]
        C=ABCD[1,0]
        D=ABCD[1,1]
        return self.contour(A,B,C,D,wl,w1)
    
    def ABCD_d1_R(self,d1,wl,w1):
        ABCD=self.free(d1) 
        A=ABCD[0,0]
        B=ABCD[0,1]
        C=ABCD[1,0]
        D=ABCD[1,1]
        return self.wavefront(A,B,C,D,wl,w1)
    
    "렌즈부터 beam waist까지"    
    def ABCD_d2(self,d1,fl,d2,wl,w1):
        ABCD=self.free(d2-d1)*self.lens(fl)*self.free(d1) 
        A=ABCD[0,0]
        B=ABCD[0,1]
        C=ABCD[1,0]
        D=ABCD[1,1]
        return self.contour(A,B,C,D,wl,w1)    
    
    def ABCD_d2_R(self,d1,fl,d2,wl,w1):
        ABCD=self.free(d2-d1)*self.lens(fl)*self.free(d1) 
        A=ABCD[0,0]
        B=ABCD[0,1]
        C=ABCD[1,0]
        D=ABCD[1,1]
        return self.wavefront(A,B,C,D,wl,w1)    

    """ABCD 매트릭스 에서 w 구하는 함수"""     
    def contour(self,A,B,C,D,wl,w1):
        q1=np.complex(0,np.pi*w1**2/wl)  #(18)식
        q2=(A*q1+B)/(C*q1+D)  #(43)식
        q2_inverse=q2**(-1)  #(43)식 응용
        w=abs(q2_inverse.imag/wl*np.pi)**(-1/2) #(17)식 응용
        return w

    """"ABCD 매트릭스에서 R구하는 함수"""
    def wavefront(self,A,B,C,D,wl,w1):
        q1=np.complex(0,np.pi*w1**2/wl)
        q2=(A*q1+B)/(C*q1+D)
        q2_inverse=q2**(-1)
        R=q2_inverse.real**(-1)
        return R

    """ABCD Matrix 함수""" 
    #Free space propagation 계산기
    def free(self,d):
        return np.matrix([[1,d],[0,1]])
    
    #thin lens 계산기
    def lens(self,fl):
        return np.matrix([[1,0],[-1/fl,1]])
              
    
class CouplingEfficiency:
    
    def __init__(self,Refl,wl,L,pzt_freq,I0,w_L):
        self.Refl = Refl  #Reflectivity
        self.wl = wl  #wavelength of beam
        self.L = L  #length of cavity
        self.pzt_freq = pzt_freq  #cavity PZT frequency
        self.I0 = I0  #Initial intensity I0를 상수로하지 않고, 함수로 만들어서 곱해주면 이야기는 끝난다. 어레이를 하나 만들어라.
        self.k = 2*np.pi/wl  #wavevector
        self.T = 1-Refl  #Transmitance
        self.v = wl/2*pzt_freq  #cavity meter velocity
        self.FSR_f = c/(2*L)  #Cavity FSR [Hz]
        self.tr = 2*L/c  #round-trip time
        self.w_L = w_L   #FWHM of beam
        self.v_freq = self.FSR_f*self.pzt_freq   #Cavity frequency velocity   
        self.mono()  # DO IT!! 
        self.lineshape() # DO IT!!
        self.ringdown() # DO IT!!

    def mono(self):
        l = np.vectorize(lambda t: int(t/self.tr))
        Finesse = int(np.pi*np.sqrt(self.Refl)/(1-self.Refl))

        self.number = np.arange(0,Finesse,1)
        self.R_num = self.Refl**self.number
        
        x = [0, 0.99, 0.999,0.9995, 0.9999,0.99993, 0.99995,0.99997, 0.99998, 0.99999] #반사율에 따른 축
        y = [15, 20, 30, 50, 80, 120, 160, 200, 200, 400]  # 반사율에 따른 감쇠시간 축 (positive 영역)
        
        x_new = np.linspace(0,1,1000000) # 반사율 0부터 1까지를 1000000으로 나눔

        func1 = abs(x_new-self.Refl) # 입력값의 x_new 내에서의 index를 찾기위해 뺐음
        func2 = np.where(func1==min(func1)) # numpy에서 해당 index를 찾아내는 스킬임. 

        #interpolate 함수를 통해서, 해당 점사이를 직선적으로 이은 function을 만든다.
        fnc = interpolate.interp1d(x,y,fill_value='extrapolate')  
    
        y_new = fnc(x_new) #만들어진 fnc 함수를 이용해서 아까 만든 x_new를 채운다.


        self.start_time = -20*micro
        self.end_time = int(y_new[func2])*micro
        self.time_split = 1000
        
        time_t = np.linspace(self.start_time,self.end_time,self.time_split)
        self.time = time_t/micro
        self.time_l = l(time_t) #시간 t에 따른 l 생성
        time_l_real = np.array([self.time_l]*len(self.number)).T  #시간별로 무한번 반복되는 요소 추가
        del(time_t)
        

        partition = time_l_real.reshape((20,int(self.time_split/20),Finesse))
        result = []

        for cal_time in partition:
            cal_1 = 2*cal_time-self.number #같은 시간에 몇번 튕김? [0, -9993, ... ] 같은 시간 안에서 0번, 1번, 2번 튕기는거
            cal_2 = self.number*cal_1 #number와 곱했음 단순 계산식
            del(cal_1)
            cal_3 = 1j*self.k*self.v*2*self.L/c*cal_2 #단순계산식임
            del(cal_2)
            cal_4 = np.exp(cal_3)
            del(cal_3)
            cal_5 = (self.R_num*cal_4).sum(axis=1)
            del(cal_4)
            cal_6 = abs(cal_5)**2
            del(cal_5)
            cal_result=self.T**2*self.I0*cal_6
            
            result.append(cal_result)
    
        result = np.array(result)
        self.result_mono = np.reshape(result, (1000,))
    
    def get_mono(self):
        return self.result_mono
    
    
    
    def g(self, omega): #lorentzian
        return 1/pi*self.w_L/2/((omega)**2+(self.w_L/2)**2)
    
    def lineshape(self):#window time -> frequency    
        start_f=-0.5*self.v_freq*(self.end_time-self.start_time)
        end_f=0.5*self.v_freq*(self.end_time-self.start_time)
        lorentz_split=self.time_split
        dw= (end_f-start_f)/lorentz_split
        xline = np.linspace(start_f,end_f,lorentz_split)
        lorentz = [self.g(x)*dw for x in xline]

        #convolution   
        convol=fftconvolve(lorentz,self.result_mono)

        #modify of convolution result
        self.result_lineshape=convol[500:1500]
        
    def get_lineshape(self):
        return self.result_lineshape
    
    def get_axis_time(self):
        return self.time
    
    
    
    def ringdown(self):
        "링다운 시그널이 시작되는 지점"
        cutpoint = np.where(self.result_lineshape==max(self.result_lineshape))[0][0]-10

        "링다운 시그널 계산"
        down_1 = 1j*self.k*self.v*self.tr*self.number**2 #단순계산식임
        down_2 = np.exp(down_1)
        down_3 = (self.R_num*down_2).sum(axis=0)
        down_4 = abs(down_3)**2
        del(down_1,down_2,down_3)

        "링다운 시그널-본식"
        ringdown_result = self.T**2*self.Refl**(2*self.time_l)*self.I0*down_4

        "Buildup 그래프와 ringdown 그래프의 혼합을 위한 과정"
        match_cal=abs(ringdown_result-self.result_lineshape[cutpoint])
        matchpoint = np.where(match_cal==min(match_cal))[0][0]
        del(down_4)
        ringdown = ringdown_result[matchpoint:]
        buildup = self.result_lineshape[:cutpoint]

        "최종 빌드업-링다운 시그널"
        result3 = np.append(buildup,ringdown)[:len(self.result_lineshape)]

        "링다운 시그널 최종 검사"
        if len(result3) != len(self.result_lineshape):
            dummy = np.array([result3[-1]]*(len(self.result_lineshape)-len(result3)))
            self.result_ringdown = np.append(result3,dummy)
    
        else:
            self.result_ringdown = result3
            
        del(self.number, self.R_num)
            
    def get_ringdown(self):
        return self.result_ringdown

class transmission:
    def __init__(self, Refl1, Refl2, L):
        "Cavity의 주파수에 따른 Transmission rate"
        self.Refl1 = Refl1
        self.Refl2 = Refl2
        self.L = L
        self.r1=np.sqrt(self.Refl1)
        self.r2=np.sqrt(self.Refl2)
        self.Finesse=int(pi*np.sqrt(abs(self.r1*self.r2))/(1-abs(self.r1*self.r2)))
        self.T_max=(1-abs(self.r1)**2)*(1-abs(self.r2)**2)/(1-abs(self.r1*self.r2))**2
        self.f_FSR = int(c/(2*self.L))
        self.FWHM = int(self.f_FSR/self.Finesse)
        self.Trans_x()
        self.Trans_xplot()
        self.Trans_y()
        
    def Trans_x(self):
        self.Transmission_x=np.arange(0,3*self.f_FSR,100)
        
    def get_Transmission_x(self):
        return self.Transmission_x
    
    def Trans_xplot(self):
        self.Transmission_xplot=np.linspace(0,3*self.f_FSR/mega,len(self.Transmission_x))
    
    def get_Transmission_xplot(self):
        return self.Transmission_xplot
    
    def Trans_y(self):
        self.Transmission_y = self.T_max/(1+(2*self.Finesse/np.pi)**2*np.sin(np.pi*self.Transmission_x/self.f_FSR)**2)  
   
    def get_Transmission_y(self):
        return self.Transmission_y
        
        
MMC()
