import stfio
import numpy as np
import matplotlib.pyplot as plt

readfilelist = open("listoffiles.txt",'r')
filenames = readfilelist.readlines()
print filenames
readfilelist.close()

def data_analysis(filename) :
    rec = stfio.read(filename.strip('\n'))
    #rec = stfio.read("300ng.atf")
    print filename
    writename = filename.rstrip('.atf\n')
    print writename, filename
    f = open(writename + '.dat','w')
    histogram = open(writename + 'histogram.dat','w') # To make histogram, it will record top point of spike

    clamp = rec[0][0]
    voltage = rec[0][2]
    limit = rec[0][0].asarray().shape[0]
    
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
  
    rs_start = 0 # real start of stimulation calcuate by spike frequency
    rs_start_flg = 0

    if_2 = [] #instantaneus frequency, 2 spike
    if_2t = [] #instantaneus frequency time

    if_3 = [] #instantaneus frequency, 3 spike
    if_3t = [] #instantaneus frequency time

    if_4 = [] #instantaneus frequency, 4 spike
    if_4t = [] #instantaneus frequency time
    
    
    print rec.dt
    
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
           
    for i in range(limit):
        if i == 0:
            pass
        elif (clamp[i-1] - threshold_c) * (clamp[i] - threshold_c) < 0:
            spike_cnt +=1
            if in_spike == True:
                in_spike = False
            else:
                in_spike = True

        if (in_spike == True) & (i != limit-1):
            if(clamp[i-1]-clamp[i]) * (clamp[i]-clamp[i+1]) < 0:
                spikelist.append(i)
#                print(i)
                #            print (i)

    i=0
    for i in range(limit):
        if (i <limit-1) &(i>0) :
            if(clamp[i-1]-clamp[i]) * (clamp[i]-clamp[i+1]) < 0:
                histogram.writelines(repr(clamp[i])+'\n')

    histogram.close()

    print s_start
    f.writelines('#Number of Spike :'+repr(spike_cnt/2)+'\n')
            
    print len(spikelist)
    print spikelist
    j=0
    for j in range(len(spikelist)):
        if s_start < spikelist[j]:
            s_start_spikelist = j

    for s_start_spikelist in range(len(spikelist)):
        i=s_start
        if i< len(spikelist)-3 :
            if (spikelist[i+3]-spikelist[i])<10000:
                if rs_start_flg ==0 :
                    rs_start = spikelist[i]
                    rs_start_flg +=1
                
    print 'rs_start : ' + repr(rs_start)
    f.writelines('#Start time of stimulation  : ' + repr(rs_start)+'\n')

    i=0
    for i in range(len(spikelist)):
        if i<len(spikelist)-1:
            if spikelist[i]-rs_start==0:
                bs_cnt = i
 
   
    f.writelines('#bs_cnt :' + repr(bs_cnt)+'\n')

    i=0
    for i in range(len(spikelist)):
        if i< len(spikelist)-1:
            a = 1 / ((spikelist[i+1]-spikelist[i])*rec.dt)
            at = (spikelist[i+1]+spikelist[i])*0.5*rec.dt
#            print repr(a)
            """
            if a>200:
                a=200
            """
            f.writelines(repr(i)+'\t'+repr(spikelist[i])+'\t'+repr(spikelist[i+1])+'\t')
            f.writelines(repr(at)+ '\t'+repr(a)+'\n')
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
    plt.savefig(writename+"instantaneusFrequency.png")
    #plt.show()
#end of data_analysis

i=0
for i in range(len(filenames)):
    data_analysis(filenames[i])

