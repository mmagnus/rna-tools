#!/usr/bin/env python
"""THIS IS NO MY CODE, I FIXED SO IT WORKS ON MY MAC @MMAGNUS
IT'S A PART OF ROSETTA RNA-TOOLS.
"""

import argparse
#import optparse
#from parse_options import get_ints

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

parser = argparse.ArgumentParser(description='Sort a slient file according to its scores and output selected ones')
parser.add_argument('file', help='Input silent file')
parser.add_argument('-select', help='Fix just type: 10, will give you models from 1 to 10. Before: Selected output file index ranked by scores, ex. 1-9 for lowest score 9 decoys') #  nargs='+',
parser.add_argument('-o', default='silent_sort.out', help='Filename of output silent file')
parser.add_argument('-term', default='score', help='Option to input a specific column to sort by (e.g. fa_atr, rms); default is score')
args = parser.parse_args()

select_idx = None
if args.select is not None:
    select_idx = []
    for i in args.select:
        #print(i)
        #get_ints(i, select_idx)
        pass
    select_idx = range(0, int(args.select))

#Reads the silent file
header = ''
scores_data = []
new_score_data = None
is_header_over = False

terms = []
for line in open(args.file):
    if not is_header_over:
        if len(line) > 5 and line[:5] == 'SCORE':
            if is_number(line.split()[1]):
                is_header_over = True
            else:
                term_col = line.split().index(args.term)

    if not is_header_over:
        header += line
    else:
        if len(line) > 5 and line[:5] == 'SCORE':
            if new_score_data is not None:
                scores_data.append(new_score_data)
            new_score_data = [0, '']
            term = float(line.split()[term_col])
            new_score_data[0] = term
            terms.append(term)
            new_score_data[1] = line
        else:
            new_score_data[1] += line

            
else:

    scores_data.append(new_score_data)
    
#Sort the data
sorted_data = sorted(scores_data, key = lambda x : x[0])

#Output
terms.sort(reverse=True)
for i, s in enumerate(terms[:int(args.select)]):
    print(i + 1, s)
                  
print(terms[:int(args.select)])
print(args.o)
out = open(args.o, 'w')
out.write(header)
if select_idx is None: #Output everything
    for data in sorted_data:
        out.write(data[1])
else: #Output selected data
    for i in select_idx:
        out.write(sorted_data[i-1][1])
out.close()
