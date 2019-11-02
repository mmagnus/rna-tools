#!/usr/bin/env python3
#-*- coding: utf-8 -*-
"""RNA Sequence with secondary structure prediction methods.

This tool takes a given sequence and returns the secondary structure prediction provided by 5 different tools: RNAfold, RNAsubopt, ipknot, contextfold and centroid_fold. You must have these tools installed. You don't have to install all tools if you want to use only one of the methods.

It's easy to add more methods of your choince to this class.

Installation
~~~~~~~~~~~~~
Depends on what tools you want to use, follow the instructions below.

ContextFold
^^^^^^^^^^^^^^^^^^^^^
https://www.cs.bgu.ac.il/~negevcb/contextfold/

needs Java. Try this on Ubuntu 14-04 https://askubuntu.com/questions/521145/how-to-install-oracle-java-on-ubuntu-14-04 Single chain only!

ViennaRNA
^^^^^^^^^^^^^^
https://www.tbi.univie.ac.at/RNA/

For OSX install from the binary Installer from the page.

ipknot OSX
^^^^^^^^^^^^^
https://github.com/satoken/homebrew-rnatools

If one encounters a problem::

    [mm] Desktop$ /usr/local/opt/bin/ipknot
    dyld: Library not loaded: /usr/local/opt/glpk/lib/libglpk.40.dylib
      Referenced from: /usr/local/opt/bin/ipknot
      Reason: image not found
    [1]    51654 abort      /usr/local/opt/bin/ipknot

the solution is::

     brew install glpk # on OSX

RNA Structure
^^^^^^^^^^^^^
http://rna.urmc.rochester.edu/

Works with 5.8.1; Jun 16, 2016.

Download http://rna.urmc.rochester.edu/RNAstructureDownload.html and untar it in ``<RNA_PDB_TOOLS>/opt/RNAstructure/``; run make, the tools will be compiled in a folder ``exe``. Set up ``DATPATH`` in your bashrc to ``<RNA_PDB_TOOLS>/opt/RNAstructure/data_tables`` ``DATAPATH=/home/magnus/work/src/rna-pdb-tools/opt/RNAstructure/data_tables/`` (read more http://rna.urmc.rochester.edu/Text/Thermodynamics.html). RNAstructure can be run with SHAPE restraints, read more http://rna.urmc.rochester.edu/Text/File_Formats.html#Constraint about the format. The file format for SHAPE reactivity comprises two columns. The first column is the nucleotide number, and the second is the reactivity. Nucleotides for which there is no SHAPE data can either be left out of the file, or the reactivity can be entered as less than -500. Columns are separated by any white space.

MC-Sym
^^^^^^^^^^^^^

FAQ
~~~~~~~~~~~~~

- Does it work for more than one chain??? Hmm.. I think it's not. You have to check on your own. --magnus

TIPS
~~~~~~~~~~~~~

Should you need to run it on a list of sequences, use the following script::

    from rna_tools import Seq
    f = open("listOfSequences.fasta")
    for line in f:
     if line.startswith('>'):
       print line,
     else:
       print line,
       s = Seq.Seq(line.strip()) # module first Seq and class second Seq #without strip this has two lines
       print s.predict_ss(method="contextfold"),
       #print s.predict_ss(method="centroid_fold")

TODO
~~~~~~~~~~~~~

- This calss should be renamed to RNASeq and merged with RNASeq class from RNAalignment

"""  # noqa
import os

try:
    RPT_PATH = os.environ['RNA_TOOLS_PATH']
except KeyError:
    print ('Set up RNA_TOOLS_PATH, see Installation note')

import subprocess
import tempfile
import sys

from rna_tools.SecondaryStructure import parse_vienna_to_pairs
from rna_tools.rna_tools_config import CONTEXTFOLD_PATH, RNASTRUCTURE_PATH, ENTRNA_PATH


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

    def __init__(self, seq, ss='', name='rna_seq'):
        self.seq = seq
        self.ss = ss
        self.ss_log = ''
        self.name = name

    def __repr__(self):
        return self.name + '\n' + self.seq + '\n' + self.ss

    def eval(self, ss='', no_dangling_end_energies=False, verbose=False):
        """Evaluate energy of RNA sequence.

        Args:
            ss (optional), if not set, then self.ss is taken for calc
            no_dangling_end_energies (Boolean)
            verbose (Boolean)

        Returns:
            Energy (float)

        The RNAeval web server calculates the energy of a RNA sequence on a given secondary structure.
        You can use it to get a detailed thermodynamic description (loop free-energy decomposition) of your RNA structures.

        Simply paste or upload your sequence below and click Proceed. To get more information on the meaning of the options click the help symbols. You can test the server using this sample sequence/structure pair.

        An equivalent RNAeval command line call would have been::

            RNAeval -v -d0 < input.txt

        Read more: <http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAeval.cgi>

        """
        tf = tempfile.NamedTemporaryFile(delete=False)
        if not ss:
            ss = self.ss

        tf.name += '.fa'
        with open(tf.name, 'w') as f:
            f.write('>' + self.name + '\n')
            f.write(self.seq + '\n')
            f.write(ss + '\n')

        dopt = ' -d2 '
        if no_dangling_end_energies:
            dopt = ' -d0 '

        cmd = 'RNAeval ' + dopt + ' < ' + tf.name
        if verbose:
            print(cmd)
        self.ss_log = subprocess.check_output(cmd, shell=True).decode()
        # [u'>rna_seq\nGGCAGGGGCGCUUCGGCCCCCUAUGCC\n((((((((.((....)).)))).))))', u'(-13.50)']
        return float(self.ss_log.strip().split(' ')[-1].replace('(','').replace(')', ''))


    def get_foldability(self, ss='', verbose=False):
        """Calculate foldability based on EntRNA.

        Steps:

        - parse SS into basepairs,
        - calculate foldabilty

        Configuration:

        - Set ENTRNA_PATH to the folder where ENTRNA_predict.py is.

        Cmd wrapper in here::

            python ENTRNA_predict.py --seq_file pseudoknotted_seq.txt --str_file pseudoknotted_str.txt

        Su, C., Weir, J. D., Zhang, F., Yan, H., & Wu, T. (2019).
        ENTRNA: a framework to predict RNA foldability. BMC Bioinformatics, 20(1), 1–11.
        http://doi.org/10.1186/s12859-019-2948-5
        """
        if ss:
            self.ss = ss

        # parse SS into base-pairs
        def dp_to_bp(dp):
            import numpy as np
            a_list = []
            bp_array = np.zeros(len(dp),dtype = int)
            for i in range(len(dp)):
                if dp[i] == "(":
                    a_list.append(i)
                if dp[i] == ")":
                    bp_array[i] = a_list[-1] + 1
                    bp_array[a_list[-1]] = i + 1
                    a_list.pop()
            return list(bp_array)

        bp = dp_to_bp(self.ss)
        if verbose: print(bp)

        tempstr = tempfile.NamedTemporaryFile(delete=False)
        with open(tempstr.name, 'w') as f:
            f.write(str(bp))

        tempseq = tempfile.NamedTemporaryFile(delete=False)
        with open(tempseq.name, 'w') as f:
            f.write(self.seq)

        # -W to silent warnings See [1]
        cmd = "cd " + ENTRNA_PATH + " && python -W ignore ENTRNA_predict.py --seq_file " + tempseq.name + " --str_file " + tempstr.name
        log = subprocess.check_output(cmd, shell=True).decode()
        if verbose:
            print(cmd)
            print(log)
        for l in log.split('\n'):
            if l.startswith('Foldability: '):
                return round(float(l.replace('Foldability: ', '')), 2)
        return -1
        ## [1]:
        ##     /Users/magnus/work/evoClustRNA/rna-foldability/ENTRNA/util/pseudoknot_free.py:22: SettingWithCopyWarning:
        ## A value is trying to be set on a copy of a slice from a DataFrame.
        ## Try using .loc[row_indexer,col_indexer] = value instead

        ## See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
        ##   df_v1['length'] = df_v1['seq'].apply(lambda x:len(x))
        ## /home/magnus/miniconda2/lib/python2.7/site-packages/sklearn/linear_model/logistic.py:433: FutureWarning: Default solver will be changed to 'lbfgs' in 0.22. Specify a solver to silence this warning.
        ##   FutureWarning)
        ## cd /Users/magnus/work/evoClustRNA/rna-foldability/ENTRNA/ && python ENTRNA_predict.py --seq_file /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpUORegp --str_file /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmp1ERCcD

    def predict_ss(self, method="RNAfold", constraints='', enforce_constraint=False, shapefn='', explore='', verbose=0):
        """Predict secondary structure of the seq.

        Args:
          method:
          onstraints:
          shapefn (str): path to a file with shape reactivites
          verbose (boolean)

        It creates a seq fasta file and runs various methods for secondary structure
        prediction. You can provide also a constraints file for RNAfold and RNAsubopt.

        Methods that can be used with contraints: RNAsubopt, RNAfold, mcfold.

        Methods that can be used with SHAPE contraints: RNAfold.

        **ContextFold**

        Example::

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


        **RNAstructure**

        Example::

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


        **MC-Fold**

        MC-Fold uses the online version of the tool, this is very powerful with constraints::

            rna_seq
            acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg
            ((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))
            curl -Y 0 -y 300 -F "pass=lucy" -F mask="((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))" -F sequence="acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
            mcfold::energy best dynamics programming: -53.91
            (-53.91, '((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))')

            curl -Y 0 -y 300 -F "pass=lucy" -F mask="((((........)))).......((((..............((((((((((..............))))))))))..))))" -F sequence="acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
            mcfold::energy best dynamics programming: -34.77
            (-34.77, '((((........)))).......((((..............((((((((((..............))))))))))..))))')

            acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg
            ((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))
            curl -Y 0 -y 300 -F "pass=lucy" -F mask="((((xxxxxxxx))))xxxxxxx((((xxxxxxxxxxxxxx((((((((((xxxxxxxxxxxxxx))))))))))xx))))" -F sequence="acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
            mcfold::energy best dynamics programming: -34.77
            (-34.77, '((((........)))).......((((..............((((((((((..............))))))))))..))))')


            acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg
            ((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))
            curl -Y 0 -y 300 -F "pass=lucy" -F mask="((((********))))*******((((**************((((((((((**************))))))))))**))))" -F sequence="acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
            mcfold::energy best dynamics programming: -77.30
            (-71.12, '(((((((..))))))).......((((((.(((...)))..(((((((((((((((....)))))))))))))))))))))')

            acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg
            ((((........)))).......((((..............(((((((((((((((....)))))))))))))))..))))
            curl -Y 0 -y 300 -F "pass=lucy" -F mask="((((**[[[[[**))))*******((((****]]]]]****(((((((((((((((****)))))))))))))))**))))" -F sequence="acucggcuaggcgaguauaaauagccgucaggccuagcgcguccaagccuagccccuucuggggcugggcgaagggucggg" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
            mcfold::energy best dynamics programming: -77.30
            ('-77.30', '((((**[[[[[**))))*******((((****]]]]]****(((((((((((((((****)))))))))))))))**))))')

        **explore**

         The sub-optimal search space can be constrained within a percentage of the minimum free energy structure, as MC-fold makes use of the Waterman-Byers algorithm [18, 19]. Because the exploration has an exponential time complexity, increasing this value can have a dramatic effect on MC-Fold’s run time.

        Parisien, M., & Major, F. (2009). RNA Modeling Using the MC-Fold and MC-Sym Pipeline [Manual] (pp. 1–84).


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

        if method == "RNAfoldX" and constraints:
            if enforce_constraint:
                cmd = 'RNAfold -p -d2 --noLP -C --enforceConstraint < ' + tf.name
            else:
                cmd = 'RNAfold -p -d2 --noLP -C < ' + tf.name
            if verbose:
                print(cmd)

            try:
                self.ss_log = subprocess.check_output(cmd, shell=True).decode()
            except subprocess.CalledProcessError:
                 print('Error')
                 return 0, 'error', 0, '', 0, '', 0, 0

            if verbose:
                print(self.ss_log)
            # parse the results
            lines = self.ss_log.split('\n')
            if 'Supplied structure constraints create empty solution set for sequence' in self.ss_log:
                return 0, 'Supplied structure constraints create empty solution set for sequence', 0, '', 0, '', 0, 0
            #if not 'frequency of mfe structure' in self.ss_log:

            # RNAfold -p -d2 --noLP -C < /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpGiUoo7.fa
            # >rna_seq
            # AAAUUAAGGGGAAGCGUUGAGCCGCUACCCAUAUGUGGUUCACUCGGAUAGCGGGGAGCUAAUAGUGAAACCGGCCCUUUAGGGG
            # ...((((((((.(((......((((((.((....(((...)))..)).))))))...)))..............))))))))... (-19.80)
            # ...{(((((((.(((......((((((.((....(((...)))..)).))))))...)))..............)))))))}... [-21.05]
            #...((((((((.(((......((((((.((....(((...)))..)).))))))...)))..............))))))))... {-19.80 d=2.34}
            # frequency of mfe structure in ensemble 0.131644; ensemble diversity 3.68

            mfess = lines[2].split()[0]
            mfe = ' '.join(lines[2].split()[-1:])
            mfe = float(mfe.replace('(', '').replace(')', '')) # (-19.80) ->-19.80

            efess = lines[3].split()[0]  # ensamble free energy
            efe = ' '.join(lines[3].split()[-1:])
            efe = float(efe.replace('[', '').replace(']', '')) # (-19.80) ->-19.80

            cfess = lines[4].split()[0]  # ensamble free energy
            cfe, d = ' '.join(lines[4].split()[1:]).split('d')
            cfe = float(cfe.replace('{', '').replace('}', '')) # (-19.80) ->-19.80

            words = lines[5].split()  # ensamble free energy
            freq = round(float(words[6].replace(';', '')), 2) # frequency of mfe structure in ensemble
            diversity = float(words[9])     # ensemble diversity

            if verbose:
                print(mfe, mfess, efe, efess, cfe, cfess, freq, diversity)
            return mfe, mfess, efe, efess, cfe, cfess, freq, diversity

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
            #  -F tope=1
            if explore:
                explore_str = " -F explore=" + str(explore)
            else:
                explore_str = ''

            #if constraints:
            #cmd = "curl -Y 0 -y 300 -F \"pass=lucy\" -F mask=\"" + constraints + "\" " + explore_str + \
            #" -F sequence=\"" + self.seq + "\" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi"
            cmd = "curl https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi\?pass\=lucy\&sequence\=" + self.seq + "\&top\=20\&explore\=15\&name\=\&mask\='" + constraints + "'\&singlehigh\=\&singlemed\=\&singlelow\="
            # cmd = "curl -Y 0 -y 300 -F \"pass=lucy\" -F sequence=\"" + self.seq + "\" https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi"
            if verbose:
                print(cmd)
                
            o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out = o.stdout.read().decode(errors='ignore').strip()
            err = o.stderr.read().decode(errors='ignore').strip()

            if verbose:
                print(out)

            # If the structure can't be find, detect this statement and finish this routine.
            if 'Explored 0 structures' in out:
                return 0.00, '', 'Explored 0 structures'

            comment = ''
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
                    energy_bdp = l.split(':')[1].strip()
                    if verbose:
                        print ('mcfold::energy best dynamics programming: ' + energy_bdp)
                        comment = 'energy best dynamics programming'
                    ss = constraints
                    # return float(energy), constraints # I'm not sure if this is good

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

            # if evenn this will not find ss, then set ss to null not to crash
            # and it's possible, like in here
            # curl -Y 0 -y 300 -F "pass=lucy" -F mask="((******)))" -F sequence="CCUgcgcaAGG" \
            #            http://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi
            ss = ''
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
                        energy = energy_bdp
                        comment = 'BP energy'
                        return energy_bdp, constraints, comment
                    # break

            # prepare outputs, return and self-s
            self.log = out
            self.ss = ss
            return float(energy), ss, comment

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

def load_fasta_ss_into_RNAseqs(fn, debug=True):
    seqs = []
    with open(fn) as f:
        for line in f:
            if debug: print(line)
            name = line.replace('>', '').strip()
            seq = next(f).strip()
            ss = next(f).strip()
            rs = RNASequence(seq, ss, name)
            seqs.append(rs)
    return seqs

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

    #
    # test of MethodNotChose
    # print(seq.predict_ss("test"))
