#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# Edited in 2015.03.23
###############################################

import stfio
import numpy as np
import matplotlib.pyplot as plt
import math
import os.path
#import seaborn as sns
from SpikeDetector import spikedetector

#readfilelist = open("ListOfData.txt",'r')
readfilelist = open("listoffiles2.txt",'r')
filenames = readfilelist.readlines()
print filenames
readfilelist.close()

RANGE = 170
w_step_time = []
w_step_freq = [0 for _ in range(RANGE)]
DATABIN     = 0.1
AVER_STIM   = 6.0
#splist:spikelist2 ,dt = rec.dt, databin = timestep? timeWidth = limit
def analysis_by_bin(wname,splist2,Bin, steps, delay, dt, timeWidth,sstart,send) :
    step_freq = [0 for _ in range(steps+1)]
    step_num = [0 for _ in range(steps+1)]
    step_time = [0 for _ in range(steps+1)]
    l_step_freq = [0 for _ in range(RANGE)]
    for j in range(steps+1):
        step_time[j] = j*Bin*dt
        if (j+1)*Bin >= timeWidth:
            for i in range(j*Bin, timeWidth):
                if splist2[i] == 1:
                    step_num[j] +=1
        else:
            for i in range(j*Bin, (j+1)*Bin):
                if splist2[i] == 1:
                    step_num[j] +=1

    for j in range(steps+1):
        step_freq[j] = step_num[j]/DATABIN
    
    #fig2 =plt.figure()
    fig2 =plt.figure(figsize=(12,8),dpi=400)
    plt.bar(step_time, step_freq,width = 0.1,color='0.2')
    plt.plot([sstart,send],[175,175],linewidth=8,color='k')
    plt.xlabel('Time[s]',fontsize=20)
    plt.xlim(0,15)
    plt.ylim(0,180)
    plt.ylabel('Frequency[Hz]',fontsize=20)
    plt.title('PSTH',fontsize=20)
    plt.savefig("./graph/"+wname+"bin.png")
    plt.close()
    bin_file = open("./analyzed_data/"+wname+"bin.dat",'w')
    bin_file.writelines(repr(steps)+' '+repr(delay)+'\n')
    for j in range(steps+1):
#        print 'j '+repr(j)
        w_step_freq[j+10+delay] +=step_freq[j]
        l_step_freq[j+10+delay] = step_freq[j]

    for k in range(RANGE):
        bin_file.writelines(repr(w_step_time[k]) + '\t' + repr(l_step_freq[k])+'\n')
    bin_file.close()

wholedata = open('./analyzed_data/ifreq_whole.dat','w')
def data_analysis(string) :
    #刺激スタート時間の補正が必要(とりあえず　全部0にした)
    filename, starttime = string.split(' ')
    rec = stfio.read(filename)
    #rec = stfio.read("300ng.atf")
    print 'filename '+filename+' start time '+ starttime
    filename = os.path.basename(filename)
    writename = filename.rstrip('.atf')
    print writename, filename
    f = open("./analyzed_data/"+writename + '.dat','w')
    histogram = open("./analyzed_data/"+writename + 'histogram.dat','w') # To make histogram, it will record top point of spike
    freqdata = open("./analyzed_data/"+writename + 'ifreq.dat','w')
    
    clamp = rec[0][0]
    voltage = rec[0][2]
    limit = rec[0][0].asarray().shape[0]

    flg=0
    for i in range(limit):
        if i==0:
            pass
        elif(voltage[i-1]-0.2)*(voltage[i]-0.2)<0:
            if flg==0:
                stim_start = i*rec.dt
                flg +=1
            elif flg==1:
                stim_end = i*rec.dt
                flg = 0
    flg=0
    
    time = np.arange(0.0, limit, 1.0) * rec.dt
    threshold_c = 0.75 #for cclamp 
    threshold_v = 0.2 #for voltage IN2
    
    spike_cnt = 0 #spike number counter
    bs_cnt = 0 #before stimulus spike number counter
    as_cnt = 0 #after stimulus spike number counter
    
    s_start = 0 # start of stimulation by Voltage
    s_end = 0 # end of stimulation by Voltage
    s_flg = 0 
    
    in_spike = False
    spikelist = [] #time of the spike will be recorded. spikelist.x[i]*dt = Spike time
    spikelist2 = [] #length of this list is same as clamp and voltage, if spike: spikelist2[i]=1; if not spike:spikelist2[i]=0

    rs_start = 0 # real start of stimulation calcuate by spike frequency
    rs_start_flg = 0

    if_2 = [] #instantaneus frequency, 2 spike
    if_2t = [] #instantaneus frequency time

    if_3 = [] #instantaneus frequency, 3 spike
    if_3t = [] #instantaneus frequency time

    if_4 = [] #instantaneus frequency, 4 spike
    if_4t = [] #instantaneus frequency time
    
    sec1spike = 0
    sec2spike = 0
    
    print rec.dt
    print "time "
    print rec.dt * limit 
    f.writelines('#Whole time :'+repr(rec.dt * limit)+'\n')
    for i in range(limit):
        if i ==0:
            pass
        elif (voltage[i-1]-threshold_v)*(voltage[i]-threshold_v)<0:
            if s_flg == 0:
 #               print(i)
                s_start = i
                s_flg +=1
            elif s_flg > 0:
                s_end = i
            
    print 's_start : ' + repr(s_start) + ' s_end : ' +repr(s_end) + ' i : ' + repr(i)
    print (limit)
    i = 0
    spikelist2 =[0 for _ in range(limit)]
          
    for i in range(limit):
        if i == 0:
            pass
        elif (clamp[i-1] - threshold_c) * (clamp[i] - threshold_c) < 0:
            spike_cnt +=1            
            if in_spike == True:
                in_spike = False
            elif in_spike == False:
                in_spike = True

        if (in_spike == True) & (i != limit-1):
            if(clamp[i-1]<clamp[i]) & (clamp[i]>clamp[i+1]):
                spikelist.append(i)
                spikelist2[i] = 1
#                print(i)
                #            print (i)

    ##Print spike timing
    spfile = open("./analyzed_data/"+writename + 'spiketiming.dat','w')
    for i in range(len(spikelist)):
        spfile.writelines(repr(spikelist[i]*rec.dt)+'\n')
    spfile.close()
    ##Print spike timing end


    i=0
    for i in range(limit):
        if (i <limit-1) &(i>0) :
            if(clamp[i-1]-clamp[i]) * (clamp[i]-clamp[i+1]) < 0:
                histogram.writelines(repr(clamp[i])+'\n')

    histogram.close()

    print s_start
    print spike_cnt
    f.writelines('#Number of Spike :'+repr(spike_cnt/2)+'\n')
  
    for j in range(len(spikelist)):
        if spikelist[j]>s_start:
            if j< len(spikelist)-3 :
                if (spikelist[j+3]-spikelist[j])<10000:
                    if rs_start_flg ==0 :
                        rs_start = spikelist[j]
                        rs_start_flg +=1
                
    print 'rs_start : ' + repr(rs_start) +'  '+ repr(rs_start*rec.dt)+'s'
    f.writelines('#Start time of stimulation  : ' + repr(rs_start)+'\n')

    i=0
    for i in range(len(spikelist)):
        if i<len(spikelist)-1:
            if spikelist[i]-rs_start==0:
                bs_cnt = i
 
   
    f.writelines('#bs_cnt :' + repr(bs_cnt)+'\n')
    for j in range(len(spikelist)):
        if (spikelist[j]-rs_start>=0)&(spikelist[j]-rs_start<(1/rec.dt)):
            sec1spike += 1

    for j in range(len(spikelist)):
        if (spikelist[j]-rs_start>=0)&(spikelist[j]-rs_start<(2/rec.dt)):
            sec2spike += 1
    f.writelines('#1sec spike : ' + repr(sec1spike)+'\n')
    f.writelines('#2sec spike : ' + repr(sec2spike)+'\n')
    
    i=0
    freqdata.writelines(repr(len(spikelist)-1-bs_cnt)+'\n')
    for i in range(len(spikelist)):
        if i< len(spikelist)-1:
            a = 1 / ((spikelist[i+1]-spikelist[i])*rec.dt)
            at = (spikelist[i+1]+spikelist[i])*0.5*rec.dt
            #print repr(a)
            """
            if a>200:
                a=200
            """
            f.writelines(repr(i)+'\t'+repr(spikelist[i])+'\t'+repr(spikelist[i+1])+'\t')
            f.writelines(repr(at)+ '\t'+repr(a)+'\n')
            if bs_cnt<=i:
                freqdata.writelines(repr(at)+ '\t'+repr(a)+'\n')
            if_2.append(a)
            if_2t.append(at)

    i=0
    for i in range(len(spikelist)):
        if i< len(spikelist)-2:
            b = 2/ ((spikelist[i+2] - spikelist[i])*rec.dt)
            bt = (spikelist[i+2]+spikelist[i+1]+spikelist[i])/3*rec.dt
            """
            if b>200:
                b=200
            """
            f.writelines(repr(bt)+'\t'+repr(b)+'\n')
            if_3.append(b)
            if_3t.append(bt)

    f.writelines('\n\n')
    i=0
    for i in range(len(spikelist)):
        if i< len(spikelist)-3:
            c = 3/ ((spikelist[i+3] - spikelist[i])*rec.dt)
            ct = (spikelist[i+3]+spikelist[i+2]+spikelist[i+1]+spikelist[i])/4*rec.dt
            """
            if c>200:
                c=200
            """
            f.writelines(repr(ct)+'\t'+repr(c)+'\n')
            if_4.append(c)
            if_4t.append(ct)

    for i in range(len(spikelist)):
        if i< len(spikelist)-5:
            d = 5/ ((spikelist[i+5] - spikelist[i])*rec.dt)
            dt = (spikelist[i+5]+spikelist[i+4]+spikelist[i+3]+spikelist[i+2]+spikelist[i+1]+spikelist[i])/6*rec.dt
            if bs_cnt<=i:
                wholedata.writelines(repr(dt)+ '\t'+repr(d)+'\n')
   
            
    Bin   = int(DATABIN/rec.dt)
    steps = int(limit/Bin)
    delay = int(starttime.strip('\n'))
    print 'delay ' +repr(delay)+'steps '+repr(steps)

    analysis_by_bin(writename, spikelist2,Bin,steps,delay,rec.dt, limit,stim_start,stim_end)

    f.close()
    flg = plt.figure()
    plt.subplot(311)
    plt.plot(if_2t,if_2)
    plt.ylim(0,300)
    plt.subplot(312)
    plt.plot(if_3t,if_3)
    plt.ylim(0,300)
    plt.subplot(313)
    plt.plot(if_4t,if_4)
    plt.ylim(0,300)
    plt.savefig("./graph/"+writename+"instantaneusFrequency.png")
    plt.close()

    #flg2__ = plt.figure()
    fig2__ =plt.figure(figsize=(12,8),dpi=400)
    plt.plot([stim_start,stim_end],[295,295],linewidth=8,color='k')
    plt.plot(if_2t,if_2,color='0.2')
    plt.ylim(0,300)
    plt.xlim(0,15)
    plt.xlabel('Time[s]',fontsize=20)
    plt.ylabel('Frequency[Hz]',fontsize=20)
    plt.title('Instantaneus Frequency',fontsize=20)
    plt.savefig("./graph/"+writename+"if.png")
    plt.close()
    #plt.show()
#end of data_analysis
        


w_step_time =[0 for _ in range(RANGE)]
for i in range(RANGE):
    w_step_time[i] = (i-10)*0.1

i=0
for i in range(len(filenames)):
    data_analysis(filenames[i])

for j in range(RANGE):
    w_step_freq[j] = w_step_freq[j]/6.0

aver_f = open('./analyzed_data/average.dat','w')

for i in range(RANGE):
    aver_f.writelines(repr(w_step_time[i])+'\t'+repr(w_step_freq[i])+'\n')

fig3 = plt.figure(figsize=(10,8),dpi=400)
fig3.patch.set_alpha(0.0)
plt.rcParams['font.size']=20
plt.bar(np.array(w_step_time)-6.0, w_step_freq,width = 0.1,color="#3b3b3b")
plt.ylim(0,120)
plt.xlim(-6,8)
plt.xlabel('Time[s]',fontsize=30)
plt.ylabel('Frequency[Hz]',fontsize=30)
plt.savefig("./graph/average.png")
#plt.close()

#### Graph for Presentation
#plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 2.0
#plt.rcParams['xtics.major.size'] = 10
#plt.rcParams['xtics.major.width'] = 1.5

fig = plt.figure(figsize=(5,4),dpi=250)
fig.subplots_adjust(bottom=0.2,left =0.20)
ax = fig.add_subplot(111)
fig.patch.set_alpha(0.0)

plt.bar(np.array(w_step_time)-6.0, w_step_freq,width = 0.1,color="#3b3b3b")
plt.ylim(0,120)
plt.xlim(-6,8)
plt.xlabel('Time[s]')
plt.ylabel('Frequency[Hz]')
plt.savefig("./graph/average_fp.png")
plt.close()
