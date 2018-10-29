#!/usr/bin/env python
#-*- coding: utf-8 -*-
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

    def eval(self, no_dangling_end_energies=True, verbose=False):
        """Evaluate energy of RNA sequence.

        Args:
            no_dangling_end_energies (Boolean)
            verbose (Boolean)

        Returns:
            Energy (flaot)

        The RNAeval web server calculates the energy of a RNA sequence on a given secondary structure. You can use it
 to get a detailed thermodynamic description (loop free-energy decomposition) of your RNA structures.

        Simply paste or upload your sequence below and click Proceed. To get more information on the meaning of the options click the help symbols. You can test the server using this sample sequence/structure pair.

        An equivalent RNAeval command line call would have been

        RNAeval -v -d0 < input.txt

        Read more: http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAeval.cgi
        """
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.name += '.fa'
        with open(tf.name, 'w') as f:
            f.write('>' + self.name + '\n')
            f.write(self.seq + '\n')
            f.write(self.ss + '\n')

        dopt = ''
        if no_dangling_end_energies:
            dopt = ' -d0 '

        cmd = 'RNAeval ' + dopt + ' < ' + tf.name
        if verbose:
            print(cmd)
        self.ss_log = subprocess.check_output(cmd, shell=True).decode()
        # [u'>rna_seq\nGGCAGGGGCGCUUCGGCCCCCUAUGCC\n((((((((.((....)).)))).))))', u'(-13.50)']
        return float(self.ss_log.strip().split(' ')[-1].replace('(','').replace(')', ''))


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

        # check for seq and constraints
        if constraints:
            if len(self.seq) != len(constraints):
                raise Exception('The seq and constraints should be of the same length: %i %s %i %s' % (len(self.seq), self.seq, len(constraints), constraints))

        # run prediction
        # rnafold without contraints
        if method == "RNAfold" and constraints:
            cmd = 'RNAfold -C < ' + tf.name
            if verbose:
                print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True).decode()
            return '\n'.join(self.ss_log.strip().split('\n')[:])

        elif method == "RNAfold":
            cmd = 'RNAfold < ' + tf.name
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

        elif method == "RNAsubopt":
            cmd = 'RNAsubopt < ' + tf.name
            if verbose:
                print(cmd)
            self.ss_log = subprocess.check_output(cmd, shell=True).decode()
            return '\n'.join(self.ss_log.split('\n')[:])

        elif method == "mcfold":
            if constraints:
                cmd = "curl -Y 0 -y 300 -F \"pass=lucy\" -F mask=\"" + constraints + "\"" + \
                " -F sequence=\"" + self.seq + "\" http://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi"
            else:
                cmd = "curl -Y 0 -y 300 -F \"pass=lucy\" -F sequence=\"" + self.seq + "\" http://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi"
            if verbose:
                print(cmd)
            o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = o.stdout.read().strip()
            err = o.stderr.read().strip()

            # If the structure can't be find, detect this statement and finish this routine.
            if 'Explored 0 structures' in out:
                return 0.00, ''

            energy = ''
            out = out.split('\n')
            for l in out :
                # first you will find the best dynamic energy, and in the next loop
                # it will be used to search for lines with this energy and secondary
                # structure

                # (((..)))  -5.43
                if energy:  # if energy is set
                    if energy in l:
                        if verbose: print(l)
                        ss = l.split()[0]

                # Performing Dynamic Programming...
                # Best Dynamic Programming Solution has Energy:  -5.43
                if l.startswith('Best Dynamic Programming Solution has Energy:'):
                    energy = l.split(':')[1].strip()
                    if verbose:
                        print ('mcfold::energy: ' + energy)

            # Ok, for whatever reason Best DP energy might not be exactly the same as and
            # the energy listed later for secondary structure. So this code finds this secondary
            # structure and gets again the energy for this secondary structure,
            # and overwrites the previous energy.
            # In this case:
            # Best Dynamic Programming Solution has Energy:  -5.46
            # ...
            # CUCUCGAAAGAUG
            # (((.((..)))))  -5.44 ( +0.00)
            # (((.((..)))))  BP >=  50%

            for l in out:
                if 'target="_blank">MARNA</a>-formatted:<P><P><P></H2><pre>' in l:
                    index = out.index(l)
                    ss_line = out[index + 2]
                    ss, energy = ss_line.split()[0:2]  # '(((.((..)))))  -5.44 ( +0.00)'

                    # if there is
                    # UUGCCGUAAGACA
                    # .............  BP >=  50%
                    # then finish with energy 0.00, and empty ss
                    if energy == 'BP':
                        return 0.00, ''
                    break

            # prepare outputs, return and self-s
            self.log = out
            self.ss = ss
            return float(energy), ss

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

        # (-51.15, '.(.(((((((((((((((..))))))))))))))))(..((((((((....)))).))))).')
        elif method == "rnastructure_CycleFold":
            cmd = RNASTRUCTURE_PATH + '/exe/CycleFold ' + tf.name + ' > ' + tf.name + '.ct '
            if verbose:
                print(cmd)
            o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stderr = o.stderr.read().strip()
            if stderr:
                print(stderr)

            # get energy
            energy = float(open(tf.name + '.ct').readline().split("energy:")[1].strip())  # >rna_seq	energy: -51.1500

            # get ss in dot-bracket notation
            cmd = RNASTRUCTURE_PATH + '/exe/ct2dot ' + tf.name + '.ct 1 ' + \
                             tf.name + '.dot'
            if verbose:
                print(cmd)
            o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stderr = o.stderr.read().strip()
            if not stderr:
                with open(tf.name + '.dot') as f:
                    # (-51.15, '.(.(((((((((((((((..))))))))))))))))(..((((((((....)))).))))).')
                    return energy, f.read().strip().split('\n')[2]


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
