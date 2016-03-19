#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# Edited in 2015.08.12
###############################################

import stfio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import math
import os.path

from SpikeDetector import spikedetector
"""
def drawGraph(filename):
    data = stfio.read(filename)
    img_filename = filename.rsplit('.',1)
    print img_filename
    clamp   = data[0][0]
    voltage = data[0][2]
    limit   = data[0][0].asarray().shape[0]
    time    = np.arange(0.0, limit, 1.0) * data.dt
    STIM   = -1
    stim   = []
    stim_t = []
    for i in range(len(voltage)):
        if(voltage[i]>0.2):
            stim.append(STIM)
            stim_t.append(time[i])
            #print time[i]
    fig = plt.figure(figsize=(10,4),dpi=400)
    plt.plot(time,clamp,'k',linewidth=0.5,color="#3d3d3d")
    #plt.plot(time,voltage,'r')
    fig.patch.set_alpha(0.0)
    plt.axis('off')
    plt.plot(stim_t,stim,'k',linewidth=8)
    img_filename="%s.png"%(img_filename[0])
    print img_filename
    plt.savefig(img_filename)
    #plt.show()
drawGraph("./original_data/3000ng.atf")     
"""
"""
if __name__=="__main__":
    if len(sys.argv) is 1:
        print "NO FILENAME"
    elif len(sys.argv) is 2:
        if(os.path.isfile(sys.argv[1])):
            sd_detect(sys.argv[1])
            #analyze(sys.argv[1])
        elif(os.path.isdir(sys.argv[1])):
            print "%s is directory"%sys.argv[1]
            target_dir = os.path.normpath(sys.argv[1])
            for fname in os.listdir(target_dir):
                full_dir = os.path.join(target_dir,fname)
                if(os.path.isfile(full_dir)):
                    #ext = os.path.splitext(full_dir)
                    sp_detect(full_dir)
        else:
            print "Wrong directory or filename"
    else:
        print "Wrong input"
"""

data = []
title = []
label = []
def ReadData(_FILENAME_):
    sd = spikedetector(_FILENAME_)
    return sd

d = ReadData("./original_data/30ng.atf")
t = "30ng"
l = "A"
label.append(l)
title.append(t)
d.get_stim_offset()
data.append(d)
d = ReadData("./original_data/100ng.atf")
t = "100ng"
l = "B"
label.append(l)
title.append(t)
d.get_stim_offset()
data.append(d)
d = ReadData("./original_data/300ng.atf")
t = "300ng"
l = "C"
label.append(l)
title.append(t)
d.get_stim_offset()
data.append(d)
d = ReadData("./original_data/1000ng.atf")
t = "1000ng"
l = "D"
label.append(l)
title.append(t)
d.get_stim_offset()
data.append(d)
d = ReadData("./original_data/3000ng.atf")
t = "3000ng"
l = "E"
label.append(l)
title.append(t)
d.get_stim_offset()
data.append(d)
d = ReadData("./original_data/08812025.atf")
t = "10000ng"
l = "F"
label.append(l)
title.append(t)
d.get_stim_offset()
data.append(d)

i=0
fig = plt.figure(figsize=(8,10),dpi=400)
font0 = FontProperties()
for i in range(6):
    plt.subplot(3,2,i+1)
    plt.plot(data[i].time,data[i].volt,'k',linewidth=0.5,color="#3d3d3d")
    #plt.plot(time,voltage,'r')
    fig.patch.set_alpha(0.0)
    plt.axis('off')
    plt.ylim(-1.5,2.0)
    plt.xlim(-1.0,15)
    """
    if(i==0):
        X = 0,0
        Y = -1.5,-1.5
        U = 1.0,0
        V = 0.0,0.5
        plt.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=1)
    """
    if(i==4):
        plt.plot([0,0],[-1.3,-0.8],'k',linewidth=3)
        plt.plot([0.0,1.0],[-1.3,-1.3],'k',linewidth=3)
        plt.text(-4.0,-1.2,"0.5mV",fontsize=15)
        plt.text(0,-1.7,"1s",fontsize=15)
    font = font0.copy()
    font.set_weight('bold')
    font.set_size(20)
    plt.text(-1.0,2,label[i],fontproperties=font)
    plt.title(title[i])
    plt.plot(data[i].dt*np.array(data[i].stim_timing[0]),[-1.2,-1.2],'k',linewidth=4)
plt.savefig("RawdataGraph.png")

#plt.show()
