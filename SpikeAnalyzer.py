#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# Edited in 2016.01.04
###############################################

import stfio
import numpy as np
import matplotlib.pyplot as plt
import math
import os.path


def raster(event_times_list, color='k'):
    """
    Creates a raster plot
    Parameters
    ----------
    event_times_list : iterable
                       a list of event time iterables
    color : string
            color of vlines
    Returns
    -------
    ax : an axis containing the raster plot
    """
    ax = plt.gca()
    for ith, trial in enumerate(event_times_list):
        plt.vlines(trial, ith + .5, ith + 1.5,linewidth=1.0, color=color)
    plt.ylim(.5, len(event_times_list) + .5)
    return ax


class spikedetector:
    stim_threshold = 0.2
    pos_threshold = 0.1
    neg_threshold = -0.2
    #neg_spikes = []
    #pos_spikes = []
    FILENAME = ""
    SAVEFILE = ""
    def __init__(self,_FILENAME_=None):
        print "** Spike Detector is Loaded **"
        print "FILE : ",_FILENAME_
        self.volt = [] # cclamp data
        self.stim = [] # IN2 data
        self.time = []
        self.stim_timing = [] # Start time of stimulus and End time of stimulus will be appended
        self.steps= -1
        self.dt   = -1
        self.pos_threshold=0.1
        self.neg_threshold=-0.2
        #self.neg_spikes = np.array[])
        #self.pos_spikes = np.array([])
        self.neg_spikes = []
        self.pos_spikes = []

        if(_FILENAME_!=None):
            self.FILENAME = _FILENAME_
            self.SAVEFILE = self.FILENAME.rsplit('.',1)
            self.SAVEFILE = self.SAVEFILE[0]
            self.read(_FILENAME_)
        print "dt : ",self.dt 

    def read(self, _FILENAME_):
        self.rec = stfio.read(_FILENAME_)
        self.volt = self.rec[0][0]
        self.volt = np.array(self.volt)
        self.stim = self.rec[0][2]
        self.steps= self.rec[0][0].asarray().shape[0]
        self.dt   = self.rec.dt
        self.time = np.arange(0.0,self.steps,1.0)*self.dt
        """
        print rec
        print rec.comment
        print rec.date
        print rec.dt
        print rec.time
        #print rec.xunints
        print rec[0].name
        print rec[0][0]
        #print rec[0].name
        print rec[0][2]
        """
    def get_stim_offset(self):
        flg=0
        for i in range(self.steps):
            if i==0:
                pass
            elif(self.stim[i-1]-self.stim_threshold)*(self.stim[i]-self.stim_threshold)<0:
                if flg==0:
                    start = i
                    flg +=1
                elif flg==1:
                    end = i
                    flg = 0
                    self.stim_timing.append([start,end])
        print "Stim offset : ",self.stim_timing

    def forward_detector(self,figure=True):
        cnt = 0
        spike_flg = False
        neg_peak = 1
        neg_timing = -1
        for i in range(self.steps):
            if i==0:
                pass
            elif(self.volt[i-1] - self.neg_threshold)*(self.volt[i] - self.neg_threshold)<0:
                if(spike_flg == True)&(self.volt[i]>self.neg_threshold):
                    spike_flg= False
                    #print "*******OUT",i,self.volt[i], self.neg_threshold, self.volt[i-1]
                    if(neg_timing>0):
                        result = self.backward_detector(neg_timing)
                        #print "R",result
                        if(result[0]>0):
                            if result[0] not in self.pos_spikes:
                                self.pos_spikes.append(result[0])
                                self.neg_spikes.append(neg_timing)
                                print result[0]*self.dt, neg_timing*self.dt
                                if(figure==True):
                                    #self.drawGraph_of_Spike(result[0],neg_timing)
                                    self.drawGraph_of_Spike_forThesis(result[0],neg_timing)
                    neg_peak = 1
                    neg_timing = -1
                elif(spike_flg == False)&(self.volt[i]<self.neg_threshold):
                    spike_flg = True
                    #print "*******IN",i,self.volt[i], self.neg_threshold, self.volt[i-1]
                    #print "F",i, spike_flg
            if(spike_flg == True) & (i != self.steps-1):
                #if(self.volt[i-1]>=self.volt[i])&(self.volt[i]<=self.volt[i+1]):
                #print "NP",i,self.volt[i]
                if(self.volt[i]<neg_peak):
                    neg_peak = self.volt[i]
                    neg_timing = i
                    #print "NP",i,neg_peak
            """
            elif(spike_flg == False):
                result = self.backward_detector(neg_timing)
                print "R",result
                if(result[0]>0):
                   pos_spikes.append(result[0])
                   neg_spikes.append(neg_timing)
                   print result[0]*self.dt, neg_timing*self.dt
            """
    def backward_detector(self,n):
        backsteps= 25
        spike_flg_= False
        pos_peak= -1
        pos_timing=-1
        #print "BD"
        for i in range(backsteps):
            if i==0:
                pass
            elif(self.volt[n-i-1]-self.pos_threshold)*(self.volt[n-i]-self.pos_threshold)<0:
                #print "POSITIVIE",n-i,self.volt[n-i-1]
                if spike_flg_ == True:
                    spike_flg_ = False
                    #print "CH TtoF"
                    return [pos_timing,pos_peak]
                elif spike_flg_ == False:
                    spike_flg_ = True
                    #print "CH FtoT"
                    
            if(spike_flg_==True)&(n-i-1 != 0):
                #print "TR"
                #if(self.volt[n-i-1]<=self.volt[n-i])&(self.volt[n-i]>=self.volt[n-i+1]):
                #print "PEAK"
                if(self.volt[n-i-1]>pos_peak):
                    pos_peak=self.volt[n-i-1]
                    pos_timing= n-i-1
                    #print "POS",pos_peak,pos_timing
            #elif(spike_flg_==False):
            #return [pos_timing,pos_peak]   
        return [pos_timing,pos_peak]
    
    def drawGraph_of_Spike(self,p_peak,n_peak):
        fig = plt.figure()
        if(p_peak-100<0):
            _min_ = 0
        else:
            _min_ = p_peak-100
        if(n_peak+100>self.steps-1):
            _max_ = self.steps-1
        else:
            _max_ =n_peak+100
        #print "Min Max",_min_,_max_
        #print self.time[_min_:_max_]
        #print self.volt
        #print self.volt[_min_:_max_]
        if not os.path.exists("./%s"%self.SAVEFILE):
            os.mkdir("./%s"%self.SAVEFILE)
        plt.plot(self.time[_min_:_max_],self.volt[_min_:_max_],'b')
        plt.plot([self.time[_min_],self.time[_max_]],[self.pos_threshold,self.pos_threshold],'g-')
        plt.plot([self.time[_min_],self.time[_max_]],[self.neg_threshold,self.neg_threshold],'g-')
        plt.plot(self.time[p_peak],self.volt[p_peak],"rx")
        plt.text(self.time[p_peak],self.volt[p_peak]+0.1,"%f,%f"%(self.volt[p_peak],self.time[p_peak]),fontsize=15)
        plt.plot(self.time[n_peak],self.volt[n_peak],"rx")
        plt.text(self.time[n_peak],self.volt[n_peak]-0.1,"%f"%self.volt[n_peak],fontsize=15)
        plt.ylim(-0.8,0.8)
        plt.xlabel("Time[s]")
        plt.ylabel("Voltage[mV]")
        plt.title("%f to %f"%(self.time[p_peak],self.time[n_peak]))
        plt.savefig("./%s/%d.png"%(self.SAVEFILE,p_peak))
        plt.close()

    def drawGraph_of_Spike_forThesis(self,p_peak,n_peak):
        #fig = plt.figure()
        #fig = plt.figure()
        fig = plt.figure(figsize=(10,8),dpi=400)
        #fig.subplots_adjust(bottom=0.2)
        #ax = fig.add_subplot(111)
        #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
        plt.rcParams['font.size'] = 20 #フォントサイズを設定
        #plt.rcParams['axes.linewidth'] = 1.5 #軸の太さを設定。目盛りは変わらない
        #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                         
        #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
        
        if(p_peak-100<0):
            _min_ = 0
        else:
            _min_ = p_peak-100
        if(n_peak+100>self.steps-1):
            _max_ = self.steps-1
        else:
            _max_ =n_peak+100
        #print "Min Max",_min_,_max_
        #print self.time[_min_:_max_]
        #print self.volt
        #print self.volt[_min_:_max_]
        if not os.path.exists("./%s"%self.SAVEFILE):
            os.mkdir("./%s"%self.SAVEFILE)
        plt.plot(self.time[_min_:_max_],self.volt[_min_:_max_],'k',linewidth=2.0)
        plt.plot([self.time[_min_],self.time[_max_]],[self.pos_threshold,self.pos_threshold],'r-',linewidth=2.0,label='threshold')
        plt.plot([self.time[_min_],self.time[_max_]],[self.neg_threshold,self.neg_threshold],'r-',linewidth=2.0)
        plt.plot(self.time[p_peak],self.volt[p_peak],"rx")
        #plt.text(self.time[p_peak],self.volt[p_peak]+0.1,"%f,%f"%(self.volt[p_peak],self.time[p_peak]),fontsize=15)
        plt.plot(self.time[n_peak],self.volt[n_peak],"rx")
        #plt.text(self.time[n_peak],self.volt[n_peak]-0.1,"%f"%self.volt[n_peak],fontsize=15)
        plt.ylim(-0.8,0.8)
        plt.legend(frameon=False)
        plt.xlabel("Time[s]")
        plt.ylabel("Voltage[mV]")
        #plt.title("%f to %f"%(self.time[p_peak],self.time[n_peak]))
        plt.savefig("./%s/%d_ft.png"%(self.SAVEFILE,p_peak))
        plt.close()
    
    def drawGraph_per_sec(self):
        x_range10 = int(1.0/self.dt)
        x_range12 = int(1.2/self.dt)
        n_graphs  = np.ceil(self.steps/float(x_range10))
        print "steps", self.steps
        if not os.path.exists("./%s"%self.SAVEFILE):
            os.mkdir("./%s"%self.SAVEFILE)
        for i in range(int(n_graphs)):
            fig = plt.figure()
            x_min = i*x_range10
            x_max = x_min+x_range12
            if(x_max>=self.steps):
                x_max=self.steps-1
            plt.plot(self.time[x_min:x_max],self.volt[x_min:x_max],'b')
            plt.plot([self.time[x_min],self.time[x_max]],[self.pos_threshold,self.pos_threshold],'g-')
            plt.plot([self.time[x_min],self.time[x_max]],[self.neg_threshold,self.neg_threshold],'g-')
            for sp in self.pos_spikes:
                if(sp>=x_min)&(sp<x_max):
                    plt.plot(self.time[sp],self.volt[sp],"rx")
                    #plt.text(self.time[sp],self.volt[sp]+0.1,"%f"%self.volt[sp],fontsize=15)
            for sp in self.neg_spikes:
                if(sp>=x_min)&(sp<x_max):
                    plt.plot(self.time[sp],self.volt[sp],"rx")
                    #plt.text(self.time[sp],self.volt[sp]-0.1,"%f"%self.volt[sp],fontsize=15)
            plt.ylim(-0.8,0.8)
            plt.xlim(1.0*i,1.0*i+1.2)
            plt.xlabel("Time[s]")
            plt.ylabel("Voltage[mV]")
            plt.title("%.2fs to %.2fs"%(x_min*self.dt,x_max*self.dt))
            plt.savefig("./%s/%dsec.png"%(self.SAVEFILE,i))
            plt.close()

    def drawGraph_per_sec_forThesis(self):
        x_range10 = int(1.0/self.dt)
        x_range12 = int(1.2/self.dt)
        n_graphs  = np.ceil(self.steps/float(x_range10))
        print "steps", self.steps
        if not os.path.exists("./%s"%self.SAVEFILE):
            os.mkdir("./%s"%self.SAVEFILE)
        for i in range(int(n_graphs)):
            #fig = plt.figure()
            fig = plt.figure(figsize=(10,8),dpi=400)
            #fig.subplots_adjust(bottom=0.2)
            #ax = fig.add_subplot(111)
            #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
            plt.rcParams['font.size'] = 20 #フォントサイズを設定
            #plt.rcParams['axes.linewidth'] = 1.5 #軸の太さを設定。目盛りは変わらない
            #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                         
            #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     

            x_min = i*x_range10
            x_max = x_min+x_range12
            if(x_max>=self.steps):
                x_max=self.steps-1
            plt.plot(self.time[x_min:x_max],self.volt[x_min:x_max],'b')
            plt.plot([self.time[x_min],self.time[x_max]],[self.pos_threshold,self.pos_threshold],'r-',linewidth=2.0)
            plt.plot([self.time[x_min],self.time[x_max]],[self.neg_threshold,self.neg_threshold],'r-',linewidth=2.0)
            for sp in self.pos_spikes:
                if(sp>=x_min)&(sp<x_max):
                    plt.plot(self.time[sp],self.volt[sp],"rx",linewidth=3.0)
                    #plt.text(self.time[sp],self.volt[sp]+0.1,"%f"%self.volt[sp],fontsize=15)
            for sp in self.neg_spikes:
                if(sp>=x_min)&(sp<x_max):
                    plt.plot(self.time[sp],self.volt[sp],"rx",linewidth=3.0)
                    #plt.text(self.time[sp],self.volt[sp]-0.1,"%f"%self.volt[sp],fontsize=15)
            plt.ylim(-0.8,0.8)
            plt.xlim(1.0*i,1.0*i+1.2)
            plt.xlabel("Time[s]")
            plt.ylabel("Voltage[mV]")
            plt.title("%.2fs to %.2fs"%(x_min*self.dt,x_max*self.dt))
            plt.savefig("./%s/%dsec_forThesis.png"%(self.SAVEFILE,i))
            plt.close()

    def save_spt(self):
        #Save spike timing file
        HEADER = "Num : %d"%(len(self.pos_spikes))
        np.savetxt("%s_spt.txt"%self.SAVEFILE,np.array(self.pos_spikes)*self.dt,fmt='%.6f',header=HEADER)
    def save_stim_info(self):
        HEADER = "STIM START\tEND"
        #print self.SAVEFILE
        np.savetxt("%s_stim.txt"%self.SAVEFILE,np.array(self.stim_timing)*self.dt,fmt='%.6f',delimiter='\t',header=HEADER)

class analyzer:
    data_dir = './bombykol_200ms/analyzed_data/'
    def __init__(self):
        self.stim_start = -1
        self.stim_end   = -1
        self.spikes     = []

    def __init__(self,Dict):
        self.stim_start = -1
        self.stim_end   = -1
        self.spikes     = []
        self.name       = Dict["name"]
        self.ID         = Dict["id"]
        self.dose       = Dict["dose"]
        self.nstims     = Dict["nstims"]
        self.name2, ext = os.path.splitext(self.name)

    def load_stim_info(self,_LOAD_):
        tmp = np.loadtxt(_LOAD_,delimiter='\t')
        self.stim_start = tmp[0]
        self.stim_end   = tmp[1]
        
    def load_spt(self,_LOAD_):
        spikes = np.loadtxt(_LOAD_)
        #print spikes
        if(self.stim_start==-1):
            print "load_stim_info should be loaded"
            return
        self.spikes = spikes-self.stim_start
        self.write_adjusted_spt()

        #print self.spikes
    def write_adjusted_spt(self):
        # spike timing was adjusted by offset of stim
        SAVE_DIR = self.data_dir+'adjusted_spt/'
        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        HEADER="STIMULUS Start: %.6f, End: %.6f"%(self.stim_start,self.stim_end)
        #print "spikes :",self.spikes
        try:
            np.savetxt("%s%s_adjusted_spt.txt"%(SAVE_DIR,self.name2),self.spikes,fmt='%.6f',delimiter='\t',header=HEADER)
        except IndexError:
            spikearray = np.array([self.spikes])
            np.savetxt("%s%s_adjusted_spt.txt"%(SAVE_DIR,self.name2),spikearray,fmt='%.6f',delimiter='\t',header=HEADER)

    def psth(self,x_range=10):
        #x_range  = 10.0 # Second
        self.bin_size = 0.1  # Second
        self.time     = np.arange(0.0,x_range,self.bin_size)
        nspikes  = [0 for i in range(len(self.time))]
        self.freqs    = [0 for i in range(len(self.time))]
        try:
            for sp in self.spikes:
                x = int(np.floor(sp/self.bin_size))
                #print sp,x
                if(x<(len(self.time)))&(x>=0):
                    nspikes[x]+=1
        except TypeError:
            x = int(np.floor(self.spikes/self.bin_size))
            if(x<(len(self.time)))&(x>=0):
                nspikes[x]+=1
        #print nspikes

        for i in range(len(nspikes)):
            self.freqs[i] = nspikes[i]/self.bin_size

        SAVE_DIR = self.data_dir+'psth/'
        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        self.writePSTH(self.freqs,self.dose,SAVE_DIR+self.name2+'_psth.dat')
        fig = plt.figure()
        plt.bar(self.time,self.freqs,width=0.1,color='b')
        plt.xlabel('Time[s]')
        plt.ylabel('Frequency[Hz]')
        plt.xlim(0,x_range)
        plt.xticks(np.arange(10))
        plt.title('NAME : %s DOSE : %d'%(self.name2,self.dose))
        plt.savefig(SAVE_DIR+self.name2+'_psth.png')
        plt.close()

    def psth_forRaster(self,x_range=10):
        #x_range  = 10.0 # Second
        self.bin_size = 0.1  # Second
        self.time     = np.arange(-1.0,x_range,self.bin_size)
        nspikes  = [0 for i in range(len(self.time))]
        self.freqs    = [0 for i in range(len(self.time))]
        try:
            for sp in self.spikes:
                x = int(np.floor(sp/self.bin_size))
                #print sp,x
                if(x<(len(self.time)-1/self.bin_size))&(x>=-10):
                    nspikes[x+int(1/self.bin_size)]+=1
                #print nspikes
        except TypeError:
            x = int(np.floor(self.spikes/self.bin_size))
            if(x<(len(self.time)-1/self.bin_size))&(x>=-10):
                nspikes[x+int(1/self.bin_size)]+=1
        #print nspikes

        for i in range(len(nspikes)):
            self.freqs[i] = nspikes[i]/self.bin_size

        SAVE_DIR = self.data_dir+'psth/'
        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        self.writePSTH(self.freqs,self.dose,SAVE_DIR+self.name2+'_psth_raster.dat')
        fig = plt.figure()
        plt.bar(self.time,self.freqs,width=0.1,color='b')
        plt.xlabel('Time[s]')
        plt.ylabel('Frequency[Hz]')
        plt.xlim(-1.0,x_range)
        plt.xticks(np.arange(-1,10))
        plt.title('NAME : %s DOSE : %d'%(self.name2,self.dose))
        plt.savefig(SAVE_DIR+self.name2+'_psth_raster.png')
        plt.close()

    def writePSTH(self,Freqs,Dose,Name=None):
        SAVE_DIR = self.data_dir+'psth/'
        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        #print self.bin_size
        #print self.time+self.bin_size/2.0
        psth_dat = np.array([self.time+self.bin_size/2.0,Freqs])
        psth_dat = psth_dat.T
        if(Name==None):
            Name='%sDose%d_psth.dat'%(SAVE_DIR,Dose)
        print Name
        np.savetxt(Name,psth_dat,fmt='%.3f',delimiter=',')
        
    def writePSTH_ALL(self,Freqs,Dose,Name=None):
        SAVE_DIR = self.data_dir+'psth/'
        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        #print self.bin_size
        #print self.time+self.bin_size/2.0
        psth_dat = np.array([self.time+self.bin_size/2.0])
        print psth_dat, len(psth_dat)
        for freqs in Freqs:
            print freqs, len(freqs)
            psth_dat=np.append(psth_dat,[freqs], axis=0)
        print psth_dat
        psth_dat = psth_dat.T
        if(Name==None):
            Name='%sDose%d_psth_all.dat'%(SAVE_DIR,Dose)
        print Name
        np.savetxt(Name,psth_dat,fmt='%.3f',delimiter=',')

    def drawPSTH(self,Freqs,Dose):
        fig = plt.figure()
        print Freqs
        plt.bar(self.time,Freqs,width=0.1,color='b')
        plt.xlabel('Time[s]')
        plt.ylabel('Frequency[Hz]')
        plt.xlim(0,10)
        plt.xticks(np.arange(10))
        plt.title('DOSE : %d'%(Dose))
        SAVE_DIR = self.data_dir+'psth/'
        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        plt.savefig('%sDose%d_psth.png'%(SAVE_DIR,Dose))
        plt.close()

    def drawPSTH_withRaster(self,Freqs,Dose,spt):
        fig = plt.figure(figsize=(10,8),dpi=400)
        fig = plt.figure()
        print Freqs
        print spt

        #fig.subplots_adjust(bottom=0.2)
        #ax = fig.add_subplot(111)
        #plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
        plt.rcParams['font.size'] = 20 #フォントサイズを設定
        #plt.rcParams['axes.linewidth'] = 1.5 #軸の太さを設定。目盛りは変わらない
        #plt.rcParams['xtics.major.size'] = 10 #x軸目盛りの長さ                         
        #plt.rcParams['xtics.major.width'] = 1.5 #x軸目盛りの太さ     
    
        ax1 = fig.add_axes((0.1, 0.75, 0.8, 0.15))
        ax1 = raster(spt,color='b')
        plt.ylabel('RNs')
        plt.title('DOSE : %d [ng]'%(Dose))
        ax2 = fig.add_axes((0.1, 0.15, 0.8, 0.60), sharex=ax1)

        # 散布図のx軸のラベルとヒストグラムのy軸のラベルを非表示
        ax1.tick_params(labelbottom="off")
        ax1.tick_params(labelleft="off")
        #ax2.tick_params(labelleft="off")

        plt.bar(self.time,Freqs,width=0.1,color='b')
        plt.xlabel('Time[s]')
        plt.ylabel('Frequency[Hz]')
        #plt.xticks(np.arange(10))
        plt.plot([0,0.200],[58,58],'k',linewidth=4.0)
        plt.xlim(-1,5)
        plt.ylim(0,60)
        plt.yticks(np.arange(0,60,10))
        SAVE_DIR = self.data_dir+'psth/'

        if not os.path.exists(SAVE_DIR):
            os.mkdir(SAVE_DIR)
        plt.savefig('%sDose%d_PSTH_withRaster.png'%(SAVE_DIR,Dose))
        plt.close()

class fitting:
    def __init__(self):
        self.psth_t = [] #Time
        self.psth_f = []#Frequency
        self.psth_f_all = [[]]
    def readPSTH(self,fname):
        tmp = np.loadtxt(fname,delimiter=',')
        tmp = tmp.T
        self.psth_t=tmp[0]
        self.psth_f=tmp[1]
        self.psth_t = np.array(self.psth_t)
        self.psth_t = self.psth_t*1000 #単位変換　sec t ms
        #print self.psth_t
        #print self.psth_f

    def readPSTH_ALL(self,fname):
        tmp = np.loadtxt(fname,delimiter=',')
        tmp = tmp.T
        self.psth_f_all=[[] for i in range(len(tmp)-1)]
        self.psth_t_all=[[] for i in range(len(tmp)-1)]
        #print "len", len(tmp)
        time = tmp[0]
        for (i,t) in enumerate(tmp):
            if(i==0):
                pass
            else:
                self.psth_f_all[i-1]=t
                self.psth_t_all[i-1]=time
        self.psth_f_all=np.array(self.psth_f_all)
        self.psth_t_all = np.array(self.psth_t_all)
        self.psth_t_all = self.psth_t_all*1000 #単位変換　sec t ms
        #print self.psth_t_all
        #print self.psth_f_all
