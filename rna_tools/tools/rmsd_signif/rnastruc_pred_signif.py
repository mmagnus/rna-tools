#!/usr/bin/python

from subprocess import Popen, PIPE
from os import remove, path, readlink

PATH = path.abspath(__file__)
if path.islink(PATH):
    PATH = path.dirname(readlink(PATH))
else:
    PATH = path.dirname(path.abspath(__file__))

BINARY_PATH = PATH + '/opt/RNA_PredictionSignificance.app'

def get_p_value(rmsd, length, verbose=False):
    cmd = BINARY_PATH + ' ' + str(length) + ' ' + str(rmsd) 
    out = Popen([cmd], stderr=PIPE, stdout=PIPE, shell=True)
    
    stdout = out.stdout.read()
    outerr = out.stderr.read()

    if verbose: print(stdout)

    pvalues = []
    for l in stdout.split('\n'):
        if l.startswith(' p-value of the prediction: '):
            pvalue = l.replace(' p-value of the prediction: ', '')
            if pvalue.startswith('<'):
                pvalue = float(pvalue[2:])
            else:
                pvalue = float(pvalue)
            pvalues.append(pvalue)
    return pvalues

if __name__ == '__main__':
    print((get_p_value(1,1000)))
    print((get_p_value(1,100)))
    print((get_p_value(1,10)))
    print((get_p_value(1,1)))
