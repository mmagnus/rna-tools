#!/usr/bin/env python

"""Select lowest energy frames out of a SimRNA trajectory file."""

from simrna_trajectory import *
import argparse

def get_parser():
        """Get parser of arguments"""
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--file', required=True)
        return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()    
    
    s = SimRNATrajectory()
    fn = args.file

    s.load_from_file(fn, top_level=True)

    #for f in s.frames:
    #    print f
    
    sorted_frames = s.sort()

    for c, f in enumerate(sorted_frames[:100]):
        print c+1,f
        #print f.header
        #print f.coords

    s2 = SimRNATrajectory()
    s2.load_from_list(sorted_frames[:100])
    s2.save(fn.replace('.trafl','') + '_100low.trafl')
