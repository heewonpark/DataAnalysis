import numpy as np
import sys
import os.path
import csv
from SpikeAnalyzer import spikedetector,analyzer

def analyze(_FILENAME_):
    sd = spikedetector()
    print _FILENAME_
    if '_spt.txt' in _FILENAME_:
        sd.load_spt(_FILENAME_)
    if '_stim.txt' in _FILENAME_:
        sd.load_stim_info(_FILENAME_)

data=[]
def load_datalist(_FILENAME_):
    reader = csv.reader(open(_FILENAME_,'r'))
    for row in reader:
        d = {"name":row[0],"dose":int(row[1]),"nstims":int(row[2]),"id":row[3]}
        data.append(d)
    #print data

def main():
    print "*** Spike Analyzer ***"
   
    load_datalist('./bombykol_200ms/analyzed_data/datalist.csv')    
    data_dir = './bombykol_200ms/analyzed_data/'
    anas = []
    for D in data:
        name, ext = os.path.splitext(D["name"])
        print name
        spt = "%s%s_spt.txt"%(data_dir,name)
        stim= "%s%s_stim.txt"%(data_dir,name)
        ana = analyzer(D)
        ana.load_stim_info(stim)
        ana.load_spt(spt)
        ana.psth_forRaster()
        anas.append(ana)

    d10   = [0.0 for i in range(110)]
    d100   = [0.0 for i in range(110)]
    d1000   = [0.0 for i in range(110)]
    d2000   = [0.0 for i in range(110)]
    d5000   = [0.0 for i in range(110)]
    d10000   = [0.0 for i in range(110)]
    d10 = np.array(d10)
    d100 = np.array(d100)
    d1000 = np.array(d1000)
    d2000 = np.array(d2000)
    d5000 = np.array(d5000)
    d10000= np.array(d10000)
    d10all   =[]
    d100all  =[]
    d1000all =[]
    d2000all =[]
    d5000all =[]
    d10000all=[]
    
    s10all = []
    s100all = []
    s1000all = []
    s2000all = []
    s5000all = []
    s10000all = []

    c10  =0
    c100 =0
    c1000=0
    c2000=0
    c5000=0
    c10000=0

    for a in anas:
        if(a.dose==10):
            d10 = d10 + a.freqs
            d10all.append(a.freqs)
            s10all.append(a.spikes)
            c10+=1
        elif(a.dose==100):
            d100 = d100 + a.freqs
            d100all.append(a.freqs)
            s100all.append(a.spikes)
            c100+=1
        elif(a.dose==1000):
            d1000 = d1000 + a.freqs
            d1000all.append(a.freqs)
            s1000all.append(a.spikes)
            c1000+=1
        elif(a.dose==2000):
            d2000 = d2000 + a.freqs
            d2000all.append(a.freqs)
            s2000all.append(a.spikes)
            c2000+=1
        elif(a.dose==5000):
            d5000 = d5000 + a.freqs
            d5000all.append(a.freqs)
            s5000all.append(a.spikes)
            c5000+=1
        elif(a.dose==10000):
            d10000 = d10000 + a.freqs
            d10000all.append(a.freqs)
            s10000all.append(a.spikes)
            c10000+=1

    print ana.bin_size
    ana.drawPSTH_withRaster(d10/c10,10,s10all)
    ana.drawPSTH_withRaster(d100/c100,100,s100all)
    ana.drawPSTH_withRaster(d1000/c1000,1000,s1000all)
    ana.drawPSTH_withRaster(d2000/c2000,2000,s2000all)
    ana.drawPSTH_withRaster(d5000/c5000,5000,s5000all)
    ana.drawPSTH_withRaster(d10000/c10000,10000,s10000all)
    print "*** END ***"

def main_30stims():
    print "*** Spike Analyzer ***"
    D = {"name":"15d29000.abf","dose":3000,"nstims":30,"id":"1"}
    #load_datalist('./atf/analyzed_data/datalist.csv')    
    data_dir = './bombykol_200ms/'
    anas = []
    #for D in data:
    name, ext = os.path.splitext(D["name"])
    print name
    spt = "%s%s_spt.txt"%(data_dir,name)
    #stim= "%s%s_stim.txt"%(data_dir,name)
    ana = analyzer(D)
    #ana.load_stim_info(stim)
    ana.stim_start=0
    ana.load_spt(spt)
    ana.psth(x_range=40)
    anas.append(ana)
    
    d3000   = [0.0 for i in range(100)]
    d3000 = np.array(d3000)
    
    c3000  =0
    """
    for a in anas:
        if(a.dose==3000):
            d3000 = d3000 + a.freqs
            c3000+=1
    ana.drawPSTH(d3000/c3000,3000)
    """
    print "*** END ***"

if __name__=="__main__":
    main()
    #main_30stims()

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
