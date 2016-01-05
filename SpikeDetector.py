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

class spikedetector:
    stim_threshold = 0.2
    pos_threshold = 0.1
    neg_threshold = -0.2
    neg_spikes = []
    pos_spikes = []
    FILENAME = ""
    def __init__(self,_FILENAME_):
        print "** Spike Detector is Loaded **"
        self.volt = [] # cclamp data
        self.stim = [] # IN2 data
        self.time = []
        self.stim_timing = [] # Start time of stimulus and End time of stimulus will be appended
        self.steps= -1
        self.dt   = -1
        self.FILENAME = _FILENAME_
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
        print self.stim_timing

    def forward_detector(self):
        cnt = 0
        spike_flg = False
        neg_peak = 1
        neg_timing = -1
        for i in range(self.steps):
            if i==0:
                pass
            elif(self.volt[i-1] - self.neg_threshold)*(self.volt[i] - self.neg_threshold)<0:
               if spike_flg == True:
                   spike_flg = False
                   #print "TR"
                   if(neg_timing>0):
                       result = self.backward_detector(neg_timing)
                       #print "R",result
                       if(result[0]>0):
                           self.pos_spikes.append(result[0])
                           self.neg_spikes.append(neg_timing)
                           print result[0]*self.dt, neg_timing*self.dt
                           self.drawGraph_of_Spike(result[0],neg_timing)
                   neg_peak = 1
                   neg_timing = -1
               elif spike_flg == False:
                   spike_flg = True

               #print "F",i, spike_flg
            if(spike_flg == True) & (i != self.steps-1):
                if(self.volt[i-1]>self.volt[i])&(self.volt[i]<self.volt[i+1]):
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
        backsteps = 20
        spike_flg = False
        pos_peak = -1
        pos_timing=-1
        #print "BD"
        for i in range(backsteps):
            if i==0:
                pass
            elif(self.volt[n-i-1]-self.pos_threshold)*(self.volt[n-i]-self.pos_threshold)<0:
                #print "POSITIVIE",n-i,self.volt[n-i-1]
                if spike_flg == True:
                    spike_flg = False
                    #print "CH TtoF"
                    return [pos_timing,pos_peak]
                elif spike_flg == False:
                    spike_flg = True
                    #print "CH FtoT"

            if(spike_flg==True):
                #print "TR"
                if(self.volt[n-i-1]<self.volt[n-i])&(self.volt[n-i]>self.volt[n-i+1]):
                    #print "PEAK"
                    if(self.volt[n-i]>pos_peak):
                        pos_peak=self.volt[n-i]
                        pos_timing= n-i
                        #print "POS",pos_peak,pos_timing
            #elif(spike_flg==False):
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
        _tmp_ = self.FILENAME.rsplit('.',1)
        if not os.path.exists("./%s"%_tmp_[0]):
            os.mkdir("./%s"%_tmp_[0])
        plt.plot(self.time[_min_:_max_],self.volt[_min_:_max_])
        plt.plot(self.time[p_peak],self.volt[p_peak],"rx")
        plt.text(self.time[p_peak],self.volt[p_peak]+0.2,"%f"%self.volt[p_peak],fontsize=15)
        plt.plot(self.time[n_peak],self.volt[n_peak],"rx")
        plt.text(self.time[n_peak],self.volt[n_peak]-0.2,"%f"%self.volt[n_peak],fontsize=15)

        plt.xlabel("Time[s]")
        plt.ylabel("Voltage[mV]")
        plt.title("%d to %d"%(p_peak,n_peak))
        plt.savefig("./%s/%d.png"%(_tmp_[0],p_peak))
        plt.close()
