#!/usr/bin/env python

"""Select lowest energy frames out of a SimRNA trajectory file. This code uses heavily the SimRNATrajectory class. Be default 100 lowest energy frames is exported."""

from .simrna_trajectory import *
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('trafl', help="SimRNA trafl file")
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()    
    
    s = SimRNATrajectory()
    fn = args.trafl

    s.load_from_file(fn, top_level=True)

    #for f in s.frames:
    #    print f
    
    sorted_frames = s.sort()

    for c, f in enumerate(sorted_frames[:100]):
        print((c+1,f))
        #print f.header
        #print f.coords

    s2 = SimRNATrajectory()
    s2.load_from_list(sorted_frames[:100])
    s2.plot_energy('subset.png')
    s2.save(fn.replace('.trafl','') + '_100low.trafl')
