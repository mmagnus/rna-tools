#!/usr/bin/env python

"""Seq and secondary structure prediction.

Installation:

- ContextFold needs java. Try this on Ubuntu 14-04 https://askubuntu.com/questions/521145/how-to-install-oracle-java-on-ubuntu-14-04 

Q: does it work for more than one chain??? Hmm.. I think it's not.
"""

import subprocess
import tempfile

class Seq:
    """Seq.

    Usage::

        >>> seq = Seq("CCCCUUUUGGGG")
        >>> seq.name = 'RNA03'
        >>> print(seq.predict_ss("RNAfold", constraints="((((....))))"))
        >RNA03
        CCCCUUUUGGGG
        ((((....)))) ( -6.40)
    """
    def __init__(self, seq):
        self.seq = seq
        self.ss = ''
        self.ss_log = ''
        self.name = 'rna_seq'
        
    def __repr__(self):
        return self.seq

    def predict_ss(self, method="RNAfold", constraints='', verbose=False):
        """Predict secondary structure of the seq.

        It creates a seq fasta file and runs various methods for secondary structure
        prediction. You can provide also a constraints file for RNAfold and RNAsubopt.
        
        ContextFold::
        
            $ java -cp bin contextFold.app.Predict in:CCCCUUUGGGGG
            CCCCUUUGGGGG
            ((((....))))
            
        It seems that a seq has to be longer than 9. Otherwise::
        
            $ java -cp bin contextFold.app.Predict in:UUUUUUGGG
            Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 10

            # this is OK
            $ java -cp bin contextFold.app.Predict in:CCCCUUUGGG
            CCCCUUUGGG
            .(((...)))
            
        """
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.name += '.fa'
        with open(tf.name,'w') as f:
            f.write('>' + self.name + '\n')
            f.write(self.seq + '\n')
            if constraints:
                f.write(constraints)    
            
        # run prediction
        if method == "RNAfold" and constraints:
            cmd = 'RNAfold -C < ' + tf.name
            if verbose: print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True)
            return '\n'.join(self.ss_log.strip().split('\n')[:])

        if method == "RNAsubopt" and constraints:
            cmd = 'RNAsubopt -C < ' + tf.name
            if verbose: print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True)
            return '\n'.join(self.ss_log.split('\n')[:])

        ## if method == "RNAsubopt":
        ##     from cogent.app.vienna_package import RNAfold, RNAsubopt
        ##     r = RNAsubopt(WorkingDir="/tmp")
        ##     res = r([self.seq])
        ##     return str(res['StdOut'].read()).strip()

        ## if method == 'RNAfold':
        ##     from cogent.app.vienna_package import RNAfold, RNAsubopt
        ##     r = RNAfold(WorkingDir="/tmp")
        ##     res = r([self.seq])
        ##     self.ss_log = res['StdOut'].read()
        ##     return self.ss_log.strip().split('\n')[-1].split()[0]

        if method == "ipknot":
            self.ss_log = subprocess.check_output('ipknot ' + tf.name, shell=True)
            return '\n'.join(self.ss_log.split('\n')[2:])

        if method == "contextfold":
            ### !!! hardcoded path ### to fix ###
            cmd = "cd /home/magnus/work/opt/ContextFold_1_00 && java -cp bin contextFold.app.Predict in:" + self.seq
            self.ss_log = subprocess.check_output(cmd, shell=True)
            return '\n'.join(self.ss_log.split('\n')[1:])
        
        if method == "centroid_fold":
            self.ss_log = subprocess.check_output('centroid_fold ' + tf.name, shell=True)
            return '\n'.join(self.ss_log.split('\n')[2:])

    def get_ss():
        if self.ss:
            return self.ss
        else:
            return self.predict_ss()
        
#main
if __name__ == '__main__':
    seq = Seq("CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG")
    seq.name = 'RNA01'
    print(seq.predict_ss("RNAfold", constraints="((((...............................................................))))"))
    seq = Seq("CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG")
    seq.name = 'RNA02'
    print(seq.predict_ss("RNAsubopt", constraints="((((...............................................................))))"))
    #print seq.predict_ss(method="ipknot")
