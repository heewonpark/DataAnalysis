import sys
import os.path

from SpikeDetector import spikedetector

"""
sd = spikedetector("./08816035.atf")
sd.get_stim_offset()
sd.forward_detector()

print "** END DETECTING **"
print sd.pos_spikes
print sd.neg_spikes
"""
def sp_detect(_FILENAME_):
    if '.atf' not in _FILENAME_:
        return
    #print _FILENAME_
    sd = spikedetector(_FILENAME_)
    sd.pos_threshold=0.1
    sd.neg_threshold=-0.2
    sd.get_stim_offset()
    sd.save_stim_info()
    #sd.forward_detector()
    #sd.drawGraph_per_sec()
    #sd.save_spt()
    print "** END DETECTING **"
    #print sd.pos_spikes
    #print sd.neg_spikes

def analyze(_FILENAME_):
    sd = spikedetector()
    print _FILENAME_
    if '_spt.txt' in _FILENAME_:
        sd.load_spt(_FILENAME_)
    if '_stim.txt' in _FILENAME_:
        sd.load_stim_info(_FILENAME_)

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

