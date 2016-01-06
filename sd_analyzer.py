import numpy as np
import sys
import os.path
import csv
from SpikeDetector import spikedetector,analyzer

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
    
    load_datalist('./atf/analyzed_data/datalist.csv')    
    data_dir = './atf/analyzed_data/'
    anas = []
    for D in data:
        name, ext = os.path.splitext(D["name"])
        print name
        spt = "%s%s_spt.txt"%(data_dir,name)
        stim= "%s%s_stim.txt"%(data_dir,name)
        ana = analyzer(D)
        ana.load_stim_info(stim)
        ana.load_spt(spt)
        ana.psth()
        anas.append(ana)

    d10   = [0.0 for i in range(100)]
    d100   = [0.0 for i in range(100)]
    d1000   = [0.0 for i in range(100)]
    d2000   = [0.0 for i in range(100)]
    d5000   = [0.0 for i in range(100)]
    d10000   = [0.0 for i in range(100)]
    d10 = np.array(d10)
    d100 = np.array(d100)
    d1000 = np.array(d1000)
    d2000 = np.array(d2000)
    d5000 = np.array(d5000)
    d10000= np.array(d10000)
    
    c10  =0
    c100 =0
    c1000=0
    c2000=0
    c5000=0
    c10000=0

    for a in anas:
        if(a.dose==10):
            d10 = d10 + a.freqs
            c10+=1
        elif(a.dose==100):
            d100 = d100 + a.freqs
            c100+=1
        elif(a.dose==1000):
            d1000 = d1000 + a.freqs
            c1000+=1
        elif(a.dose==2000):
            d2000 = d2000 + a.freqs
            c2000+=1
        elif(a.dose==5000):
            d5000 = d5000 + a.freqs
            c5000+=1
        elif(a.dose==10000):
            d10000 = d10000 + a.freqs
            c10000+=1
    ana.drawPSTH(d10/c10,10)
    ana.drawPSTH(d100/c100,100)
    ana.drawPSTH(d1000/c1000,1000)
    ana.drawPSTH(d2000/c2000,2000)
    ana.drawPSTH(d5000/c5000,5000)
    ana.drawPSTH(d10000/c10000,10000)
    print "*** END ***"

if __name__=="__main__":
    main()

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