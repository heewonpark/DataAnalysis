#! /usr/bin/python
# coding: UTF-8

##############################################
# This file is written by Park
# Edited in 2015.08.12
###############################################

import stfio
import numpy as np
import matplotlib.pyplot as plt
import math
import os.path

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
