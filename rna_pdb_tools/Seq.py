#!/usr/bin/env python

"""RNA Sequence with secondary structure prediction methods.

This tool takes a given sequence and returns the secondary structure prediction provided by 5 different tools: RNAfold, RNAsubopt, ipknot, contextfold and centroid_fold. You must have these tools installed. You don't have to install all tools if you want to use only one of the methods.

It's easy to add more methods of your choince to this class.

Installation:

ContextFold (https://www.cs.bgu.ac.il/~negevcb/contextfold/) needs Java. Try this on Ubuntu 14-04 https://askubuntu.com/questions/521145/how-to-install-oracle-java-on-ubuntu-14-04 Single chain only!

ViennaRNA (https://www.tbi.univie.ac.at/RNA/)

ipknot OSX (https://github.com/satoken/homebrew-rnatools)

RNAStructure (http://rna.urmc.rochester.edu/)

Works with 5.8.1; Jun 16, 2016.

Download http://rna.urmc.rochester.edu/RNAstructureDownload.html and untar it in ``<RNA_PDB_TOOLS>/opt/RNAstructure/``; run make, the tools will be compiled in a folder ``exe``. Set up ``DATPATH`` in your bashrc to ``<RNA_PDB_TOOLS>/opt/RNAstructure/data_tables`` ``DATAPATH=/home/magnus/work/src/rna-pdb-tools/opt/RNAstructure/data_tables/`` (read more http://rna.urmc.rochester.edu/Text/Thermodynamics.html). RNAstructure can be run with SHAPE restraints, read more http://rna.urmc.rochester.edu/Text/File_Formats.html#Constraint about the format. The file format for SHAPE reactivity comprises two columns. The first column is the nucleotide number, and the second is the reactivity. Nucleotides for which there is no SHAPE data can either be left out of the file, or the reactivity can be entered as less than -500. Columns are separated by any white space.

FAQ:

- Does it work for more than one chain??? Hmm.. I think it's not. You have to check on your own. --magnus

TIPS:

Should you need to run it on a list of sequences, use the following script::

    from rna_pdb_tools import Seq
    f = open("listOfSequences.fasta")
    for line in f:
     if line.startswith('>'):
       print line,
     else:
       print line,
       s = Seq.Seq(line.strip()) # module first Seq and class second Seq #without strip this has two lines
       print s.predict_ss(method="contextfold"),
       #print s.predict_ss(method="centroid_fold")

@todo should be renamed to RNASeq, and merged with RNASeq class from RNAalignment.
"""  # noqa
import os

try:
    RPT_PATH = os.environ['RNA_PDB_TOOLS']
except KeyError:
    print ('Set up RNA_PDB_TOOLS, see Installation note')

import subprocess
import tempfile
import sys


from rna_pdb_tools.rpt_config import CONTEXTFOLD_PATH, RNASTRUCTURE_PATH


class MethodNotChosen(Exception):
    pass


class RNASequence(object):
    """RNASequence.

    Usage::

        >>> seq = RNASequence("CCCCUUUUGGGG")
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

    def predict_ss(self, method="RNAfold", constraints='', shapefn='', verbose=0):
        """Predict secondary structure of the seq.

        :param method:
        :param constraints:
        :param shapefn: path to a file with shape reactivites
        :param verbose:

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


        RNAstructure::

            >>> seq = RNASequence("GGGGUUUUCCC")
            >>> print(seq.predict_ss("rnastructure"))
            >  ENERGY = -4.4  rna_seq
            GGGGUUUUCCC
            ((((...))))

        and with the shape data::

            >>> print(seq.predict_ss("rnastructure", shapefn="data/shape.txt"))
            >  ENERGY = -0.2  rna_seq
            GGGGUUUUCCC
            .(((....)))

        the shape data::

            1 10
            2 1
            3 1

        You can easily see that the first G is unpaired right now! The reactivity of this G was
        set to 10. Worked!
        """
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.name += '.fa'
        with open(tf.name, 'w') as f:
            f.write('>' + self.name + '\n')
            f.write(self.seq + '\n')
            if constraints:
                f.write(constraints)

        # run prediction
        if method == "RNAfold" and constraints:
            cmd = 'RNAfold -C < ' + tf.name
            if verbose:
                print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True).decode()
            return '\n'.join(self.ss_log.strip().split('\n')[:])

        elif method == "RNAsubopt" and constraints:
            cmd = 'RNAsubopt -C < ' + tf.name
            if verbose:
                print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True).decode()
            return '\n'.join(self.ss_log.split('\n')[:])

        # if method == "RNAsubopt":
        #     from cogent.app.vienna_package import RNAfold, RNAsubopt
        #     r = RNAsubopt(WorkingDir="/tmp")
        #     res = r([self.seq])
        #     return str(res['StdOut'].read()).strip()

        # if method == 'RNAfold':
        #     from cogent.app.vienna_package import RNAfold, RNAsubopt
        #     r = RNAfold(WorkingDir="/tmp")
        #     res = r([self.seq])
        #     self.ss_log = res['StdOut'].read()
        #     return self.ss_log.strip().split('\n')[-1].split()[0]

        elif method == "ipknot":
            self.ss_log = subprocess.check_output('ipknot ' + tf.name, shell=True)
            return '\n'.join(self.ss_log.split('\n')[2:])

        elif method == "contextfold":
            if not CONTEXTFOLD_PATH:
                print('Set up CONTEXTFOLD_PATH in configuration.')
                sys.exit(0)
            cmd = "cd " + CONTEXTFOLD_PATH + \
                  " + && java -cp bin contextFold.app.Predict in:" + self.seq
            if verbose:
                print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True).decode()
            return '\n'.join(self.ss_log.split('\n')[1:])

        elif method == "centroid_fold":
            self.ss_log = subprocess.check_output('centroid_fold ' + tf.name, shell=True)
            return '\n'.join(self.ss_log.split('\n')[2:])

        elif method == "rnastructure":
            cmd = RNASTRUCTURE_PATH + '/exe/Fold ' + tf.name + ' ' + tf.name + '.out '
            if shapefn:
                cmd += ' -sh ' + shapefn
            if verbose:
                print(cmd)
            o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stderr = o.stderr.read().strip()
            if stderr:
                print(stderr)

            cmd = RNASTRUCTURE_PATH + '/exe/ct2dot ' + tf.name + '.out 1 ' + \
                             tf.name + '.dot'
            if verbose:
                print(cmd)
            o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stderr = o.stderr.read().strip()
            if not stderr:
                with open(tf.name + '.dot') as f:
                    return f.read().strip()
        else:
            raise MethodNotChosen('You have to define a correct method to use.')


# main
if __name__ == '__main__':
    import doctest
    doctest.testmod()

    seq = RNASequence("CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG")
    seq.name = 'RNA01'
    print(seq.predict_ss("RNAfold",
                          constraints="((((...............................................................))))"))  # noqa
    seq = RNASequence("CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG")
    seq.name = 'RNA02'
    print(seq.predict_ss("RNAsubopt",
                         constraints="((((...............................................................))))"))  # noqa

    print(seq.predict_ss("contextfold"))
    # print seq.predict_ss(method="ipknot")

    verbose = False
    seq = RNASequence("GGGGUUUUCCC")
    print(seq.predict_ss("rnastructure", verbose=verbose))
    print(seq.predict_ss("rnastructure", shapefn="data/shape.txt", verbose=verbose))

    seq = RNASequence("CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGgUGCAA")
    print(seq.predict_ss("rnastructure", shapefn="data/shape.txt", verbose=verbose))


    # test of MethodNotChose
    # print(seq.predict_ss("test"))
