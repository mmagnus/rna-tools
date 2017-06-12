#!/usr/bin/env python

import subprocess
import os
import tempfile
from rpt_config import RFAM_DB_PATH
from Seq import RNASequence

class RfamSearchError(Exception):
    pass

class RfamSearch():
    """RfamSearch (local).

    Infernal cmscan is used to search the CM-format Rfam database.

    Set up ``RFAM_DB_PATH``

    Install http://eddylab.org/infernal/

    Cite: Nawrocki and S. R. Eddy, Infernal 1.1: 100-fold faster RNA homology searches, Bioinformatics 29:2933-2935 (2013). """
    def __init__(self):
        pass

    def cmscan(self, seq):
        """Run cmscan on the seq.

        :param seq: string
        :returns: result
        :rtype: string """
        print(seq)

        # make tmp file
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.name += '.fa'
        with open(tf.name, 'w') as f:
            f.write('>target\n')
            f.write(seq.seq + '\n')

        # make output file
        of = tempfile.NamedTemporaryFile(delete=False)

        # run cmscan
        cmd = 'cmscan -E 1 ' + RFAM_DB_PATH + ' ' + tf.name + '  > ' + of.name
        o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = o.stdout.read().strip()
        err = o.stderr.read().strip()
        if err: raise RfamSearchError(err)
        self.output = open(of.name).read()
        #os.chdir(old_pwd)
        return self.output
#main
if __name__ == '__main__':
    seq = RNASequence("GGCGCGGCACCGUCCGCGGAACAAACGG")
    rs = RfamSearch()
    hit = rs.cmscan(seq)
    print hit
