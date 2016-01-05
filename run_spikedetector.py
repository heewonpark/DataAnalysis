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

def main(_FILENAME_):
    if '.atf' not in _FILENAME_:
        return
    print _FILENAME_
    sd = spikedetector(_FILENAME_)
    sd.get_stim_offset()
    sd.forward_detector()
    print "** END DETECTING **"
    print sd.pos_spikes
    print sd.neg_spikes


if __name__=="__main__":
    if len(sys.argv) is 1:
        print "NO FILENAME"
    elif len(sys.argv) is 2:
        if(os.path.isfile(sys.argv[1])):
            main(sys.argv[1])
        elif(os.path.isdir(sys.argv[1])):
            print "%s is directory"%sys.argv[1]
            target_dir = os.path.normpath(sys.argv[1])
            for fname in os.listdir(target_dir):
                full_dir = os.path.join(target_dir,fname)
                if(os.path.isfile(full_dir)):
                    #ext = os.path.splitext(full_dir)
                    main(full_dir)
        else:
            print "Wrong directory or filename"
    else:
        print "Wrong input"

