#!/usr/bin/env python2
"""RNAalignment - a module to work with RNA sequence alignments.

To see a full demo what you can do with this util, please take a look at the jupiter notebook (https://github.com/mmagnus/rna-pdb-tools/blob/master/rna_tools/tools/rna_alignment/rna_alignment.ipynb)

Load an alignment in the Stockholm::

    alignment = ra.RNAalignment('test_data/RF00167.stockholm.sto')

 or fasta format::

     import rna_alignment as ra
     alignment = ra.fasta2stokholm(alignment.fasta)
     alignment = ra.RNAalignment

Parameters of the aligmnent::

     print(alignment.describe())

Consensus SS::

    print(alignment.ss_cons_with_pk)

Get sequnce/s from teh aligment::

    >>> seq = a.io[0]

"""

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator
from rna_tools import SecondaryStructure
from rna_tools.rna_tools_config import RCHIE_PATH
from collections import OrderedDict

import Levenshtein
import tempfile
import subprocess
import os
import shutil
import re
import gzip
import copy


class RNAalignmentError(Exception):
    pass


class RChieError(Exception):
    pass


class RFAMFetchError(Exception):
    pass


class RChie:
    """RChie - plotting arc diagrams of RNA secondary structures.

    .. image:: ../pngs/rchie.png

    http://www.e-rna.org/r-chie/

    The offline version of R-chie, which requires first installing R4RNA is available here, or clone our git repository here
    How to install it:

    - Ensure R is installed already, or download it freely from http://www.r-project.org/
    - Download the R4RNA (https://github.com/jujubix/r-chie), open R and install the package::

            install.packages("<path_to_file>/R4RNA", repos = NULL, type="source")
            # Install the optparse and RColorBrewer
            install.packages('optparse')
            install.packages('RColorBrewer')

    - Go to rna_tools/rna_tools_config_local.py and set RCHIE_PATH to the folder with RChie, e.g. ``"/home/magnus/work/opt/r-chie/"``.

    To test if Rchie works on your machine (from rna_align folder)::

      <path to your rchie>/rchie.R --msafile test_data/rchie_test_files/fasta.txt test_data/rchie_test_files/helix.txt

    you should have rchie.png file in the folder.

    More at http://www.e-rna.org/r-chie/download.cgi

    Cite: Daniel Lai, Jeff R. Proctor, Jing Yun A. Zhu, and Irmtraud M. Meyer (2012) R-chie: a web server and R package for visualizing RNA secondary structures. Nucleic Acids Research, first published online March 19, 2012. doi:10.1093/nar/gks241
    """

    def __init__(self):
        pass

    def plot_cov(self, seqs, ss_cons, plot_fn='rchie.png', verbose=False):
        """Plot an RChie plot_conv.

        :param seqs: seqs in the fasta format
        :param ss_cons: a string of secondary structure consensus, use only ``().``. Works with pseuoknots.
        """
        fasta_alignment = tempfile.NamedTemporaryFile(delete=False)
        with open(fasta_alignment.name, 'w') as f:
            f.write(seqs)

        plot = tempfile.NamedTemporaryFile(delete=False)
        plot.name += '.png'
        ss = tempfile.NamedTemporaryFile(delete=False)
        with open(ss.name, 'w') as f:
            f.write(ss_cons)
        if not RCHIE_PATH:
            raise RChieError('RChie path not set up!')
        cmd = RCHIE_PATH + \
            "rchie.R --msafile='%s' --format1 vienna '%s' --output '%s'" % (
                fasta_alignment.name, ss.name, plot.name)
        if verbose:
            print(cmd)
        o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = o.stdout.read().strip()
        err = o.stderr.read().strip()
        # *****PROCESSING MSA FILE*****
        # Error in readFasta(opt$msafile, filter = TRUE) : no FASTA sequences found
        # Error: ERROR: Invalid FASTA file
        # Execution halted
        if "error" in str(err).lower():
            raise Exception('\n'.join([cmd, err]))
        if verbose:
            print('\n'.join([cmd, err]))
        self.plotfn = plot.name
        if verbose:
            print(self.plotfn)
        if plot_fn:
            shutil.move(plot.name, plot_fn)
            print('Rchie: plot saved to %s' % plot_fn)
        from IPython.display import Image
        return Image(filename=plot_fn)

    def show(self):
        from IPython.display import Image
        return Image(filename=self.plotfn)

    def write(self, outfn):
        shutil.copyfile(self.plotfn, outfn)
        print('Write to %s' % outfn)


class RNASeq(object):
    """RNASeq.

    Args:

       id (str)    : id of a sequence
       seq (str)   : seq, it be uppercased.
       ss (str)    : secondary structure, default None

    Attributes:

       seq_no_gaps(str) : seq.replace('-', '')
       ss_no_gaps(str)  : ss.replace('-', '')

    """

    def __init__(self, id, seq, ss=None):
        self.id = id
        self.seq = seq.upper()

        # self.ss_raw = ss  # this will not be changed after remove_gaps.
        # so maybe don't use ss_raw at call
        self.ss = ss
        self.ss = self.get_ss_std()

        self.seq_no_gaps = seq.replace('-', '')
        self.ss_no_gaps = ss.replace('-', '')

    #@property
    def get_ss_std(self):
        nss = ''
        for s in self.ss:
            nss += get_rfam_ss_notat_to_dot_bracket_notat(s)
        return nss

    def __repr__(self):
        return self.id

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, i):
        if self.ss:
            return RNASeq(self.id + '_slice', self.seq[i], self.ss[i])
        else:
            return RNASeq(self.id + '_slice', self.seq[i])

    def remove_columns(self, to_remove):
        """indexing from 0"""
        nseq = ''
        for i, s in enumerate(self.seq):
            if i not in to_remove:
                nseq += s

        nss = ''
        if self.ss:
            for i, s in enumerate(self.ss):
                if i not in to_remove:
                    nss += s
        self.seq = nseq
        self.ss = nss

    def draw_ss(self, title='', verbose=False, resolution=1.5):
        """Draw secondary structure of RNA with VARNA.

        VARNA: Visualization Applet for RNA
        A Java lightweight component and applet for drawing the RNA secondary structure

        .. image :: ../pngs/varna.png

        Cite: VARNA: Interactive drawing and editing of the RNA secondary structure Kevin Darty, Alain Denise and Yann Ponty Bioinformatics, pp. 1974-197,, Vol. 25, no. 15, 2009

        http://varna.lri.fr/"""
        drawfn = tempfile.NamedTemporaryFile(delete=False).name + '.png'
        SecondaryStructure.draw_ss(title, self.seq, self.ss, drawfn, resolution, verbose=verbose)
        from IPython.display import Image
        return Image(filename=drawfn)

    def remove_gaps(self, check_bps=True, only_canonical=True, allow_gu=True):
        """Remove gaps from seq and secondary structure of the seq.

        Args:

            check_bps (bool)      : fix mistakes as
            only_canonical (bool) : keep in ss only pairs GC, AU
            allow_gu (bool)       : keep in ss also GU pair

        .. image :: ../pngs/ss_misgap.png

        A residue "paired" with a gap.

        .. image :: ../pngs/ss_misgap_wrong.png

        .. when check_bps (by default), then when removing gaps, the functions check is the gap is
        paired with any residues (in the blue circle). If yes, then this residues is unpair (in this case ``)`` -> ``.``).

        .. image :: ../pngs/ss_misgap_ok.png

        if ``only_canonical`` (by default) is True then only GC, AU can be paired.


        .. image :: ../pngs/ss_only_canonical.png


        If ``allow_gu`` is False (be default is True) then GU pair is also possible.

        .. image :: ../pngs/ss_no_gu.png

        If you provide seq and secondary structure such as::

             GgCcGGggG.GcggG.cc.u.aAUACAAuACCC.GaAA.GGGGAAUAaggCc.gGCc.gu......CU.......uugugcgGUuUUcaAgCccCCgGcCaCCcuuuu
             (((((((((....((.((............(((......)))......))))..(((.(.....................)))).......)))))))))........

        gaps will be remove as well.
        """
        GAPS = ['-', '.']

        nseq = ''
        nss = ''
        for (nt, nt_ss) in zip(self.seq, self.ss):
            if nt in GAPS and nt_ss in GAPS:
                pass
            else:
                nseq += nt
                nss += nt_ss
        self.seq = nseq
        self.ss = nss

        nss = list()
        bps = self.ss_to_bps()

        if check_bps:
            nss = list(self.ss)
            for bp in bps:
                nt_left = self.seq[bp[0]]
                nt_right = self.seq[bp[1]]
                if nt_left == '-' or nt_right == '-':
                    nss[bp[0]] = '.'
                    nss[bp[1]] = '.'
            self.ss = ''.join(nss)

        if only_canonical:
            nseq = list(self.seq)
            nss = list(self.ss)
            for bp in bps:
                nt_left = nseq[bp[0]]
                nt_right = nseq[bp[1]]
                if (nt_left == 'A' and nt_right == 'U') or (nt_left == 'U' and nt_right == 'A'):
                    pass
                elif (nt_left == 'G' and nt_right == 'C') or (nt_left == 'C' and nt_right == 'G'):
                    pass
                elif (nt_left == 'G' and nt_right == 'U') or (nt_left == 'U' and nt_right == 'G'):
                    if allow_gu:
                        pass
                    else:
                        nss[bp[0]] = '.'
                        nss[bp[1]] = '.'
                else:
                    nss[bp[0]] = '.'
                    nss[bp[1]] = '.'
            self.ss = ''.join(nss)

        # two?????? what is "two"?
        nss = []
        nseq = ''
        for i, (c, s) in enumerate(zip(self.seq, self.ss)):
            if c != '-':
                nseq += c
                nss.append(s)

        self.seq = nseq
        self.ss = ''.join(nss)

    def ss_to_bps(self):
        """Convert secondary structure into a list of basepairs.

        Returns:

            bps (list): a list of base pairs, e.g. [[0, 80], [1, 79], [2, 78], [4, 77], [6, 75], [7, 74], ...]

        """
        j = []
        bps = []
        pair_types = ['()', '[]', '<>', '{}']
        for pair_type in pair_types:
            for i, s in enumerate(self.ss):
                if s == pair_type[0]:
                    j.append(i)
                if s == pair_type[1]:
                    bps.append([j.pop(), i])
            if len(j):
                # if something left, this is a problem (!!!)
                raise Exception('Mis-paired secondary structure')
        bps.sort()
        return bps

    def get_conserved(self, consensus, start=0, to_pymol=True, offset=0):
        """Start
        UCGGGGUGCCCUUCUGCGUG--------------------------------------------------AAGGC-UGAGAAAUACCCGU-------------------------------------------------AUCACCUG-AUCUGGAU-AAUGC
        XXXXXXXXXXXXGXGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX----------------------------XXXXX-XCUGAGAXXXXXXXXXXXXXXXXXXXXXX----------------------------------XXXXXXXX-XXXXXXXX-ACXUG
        """
        c = start + offset
        index = []
        print(self.seq)
        print(consensus)
        for nt_seq, nt_consensus in zip(self.seq, consensus):
            if nt_consensus in ['-', 'X']:
                pass
            else:
                index.append(c)
            if nt_seq != '-':
                c += 1
        if to_pymol:
            return("color red, " + str(index).replace(', ', '+').replace('[','').replace(']',''))
        else:
            return index

    def get_distance_to(self, nseq):
        """Get distance of self.seq to nseq."""
        return round(Levenshtein.ratio(self.seq, nseq), 2)

class RNAalignment(object):
    """RNA alignment - adapter class around BioPython to do RNA alignment stuff

    Usage (for more see IPython notebook https://github.com/mmagnus/rna-tools/blob/master/rna_tools/tools/rna_alignment/rna_alignment.ipynb)

    >>> a = RNAalignment('test_data/RF00167.stockholm.sto')
    >>> print(a.tail())
    >>> print(a.ss_cons)

    Args:

      fn (str): Filename
      io (Bio.AlignIO): AlignIO.read(fn, "stockholm")
      lines (list): List of all lines of fn
      seqs (list): List of all sequences as class:`RNASeq` objects
      rf (str): ReFerence annotation, the consensus RNA sequence

    Read more:

    - http://biopython.org/DIST/docs/api/Bio.AlignIO.StockholmIO-module.html

    and on the format itself

    - https://en.wikipedia.org/wiki/Stockholm_format
    - http://sonnhammer.sbc.su.se/Stockholm.html

    .. warning:: fetch requires urllib3
    """

    def __init__(self, fn='', fetch=''):
        if fetch:
            import urllib3
            http = urllib3.PoolManager()
            response = http.request('GET', 'http://rfam.xfam.org/family/' +
                                    fetch + '/alignment/stockholm?gzip=1&download=1')
            if not response.status == 200:
                raise RFAMFetchError(
                    "The alignment could not be downloaded. Please check the RFAM id that you requested! (don't put .stk etc in the id)")
            with open(fetch + '.stk.gz', 'wb') as f:
                f.write(response.data)
            with gzip.open(fetch + '.stk.gz', 'rb') as f:
                file_content = f.read()
            with open(fetch + '.stk', 'wb') as f:
                f.write(file_content)
            fn = fetch + '.stk'

        self.fn = fn
        self.lines = open(fn).read().split('\n')
        self.io = AlignIO.read(fn, "stockholm")
        self.ss_cons = self.get_ss_cons()
        self.ss_cons_pk = self.get_ss_cons_pk()
        self.copy_ss_cons_to_all()
        self._ss_cons_std = self.ss_cons
        self.rf = self.get_gc_rf()
        self.shift = self.get_shift_seq_in_align()

        # get all lines not # nor //
        # fix for blocked alignment
        seq_lines = [l for l in self.lines if (
            not l.startswith('#')) and (not l.startswith('//')) and (l)]
        seqs_dict = OrderedDict()
        for seq in seq_lines:
            seq_id, seq_seq = seq.split()
            if seq_id not in seqs_dict:
                seqs_dict[seq_id] = seq_seq
            else:
                seqs_dict[seq_id] += seq_seq
        self.seqs = []
        for seq in seqs_dict:
            self.seqs.append(RNASeq(seq, seqs_dict[seq], ss=self.ss_cons_with_pk_std))

        # this is sick!
        # I create a Cols object to be able to slice alignments
        class Cols:
            def __init__(self, alignment):
                self.alignment = alignment

            def __getitem__(self, i):
                """Return new alignment"""
                if type(i) is list:
                    pass
                else:
                    # collect "new" sequences
                    n_seqs = []
                    for s in self.alignment:
                        new_seq = RNASeq(s.id, s.seq[i], s.ss[i])
                        n_seqs.append(new_seq)

                    # this is not very smart :(
                    # save new seqs to a file
                    # and load it as RNAalignment
                    tf = tempfile.NamedTemporaryFile(delete=False)
                    tf.name += '.stk'
                    with open(tf.name, 'w') as f:
                        f.write('# STOCKHOLM 1.0\n')
                        for s in n_seqs:
                            f.write(' '.join([s.id, s.seq, '\n']))
                        # add ss_cons & //
                        f.write('#=GC SS_cons ' + self.alignment.ss_cons[i] + '\n')
                        if self.alignment.ss_cons_pk:
                            f.write('#=GC SS_cons_pk' + self.alignment.ss_cons_pk[i] + '\n')
                        f.write('#=GC RF ' + self.alignment.rf[i] + '\n')
                        f.write('//\n')
                    return RNAalignment(tf.name)
        self.cols = Cols(self)
        # ^^^^ sick ^^^^^^^^^^^

    def reload_alignment(self):
        tmpfn = tempfile.NamedTemporaryFile(delete=False).name
        self.write(tmpfn)
        self.io = AlignIO.read(tmpfn, "stockholm")

    def __len__(self):
        """Return length of all sequenes."""
        return len(self.seqs)

    def __getitem__(self, i):
        if type(i) is str:
            for s in self:
                if s.id == i:
                    return s
        elif type(i) is list:
            seqs = []
            for j in i:
                seqs.append(self.seqs[j])
            return seqs
        else:
            return self.seqs[i]

    @property
    def ss_cons_std(self):
        return get_rfam_ss_notat_to_dot_bracket_notat(self.ss_cons)

    @ss_cons_std.setter
    def ss_cons_std(self, ss):
        self._ss_cons_std = ss
        print(self._ss_cons_std)

    def subset(self, ids, verbose=False):
        """Get subset for ids::

            # STOCKHOLM 1.0
            #=GF WK Tetrahydrofolate_riboswitch
            ..
            AAQK01002704.1/947-1059              -U-GC-AAAAUAGGUUUCCAUGC..
            #=GC SS_cons                         .(.((.((----((((((((((...
            #=GC RF                              .g.gc.aGAGUAGggugccgugc..
            //

        """
        nalign = ''
        for l in self.lines:
            if l.startswith('//'):
                nalign += l + '\n'
            if l.startswith('#'):
                nalign += l + '\n'
            else:
                for i in ids:
                    if l.startswith(i):
                        nalign += l + '\n'

        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.name += '.stk'
        print('Saved to ', tf.name)
        if verbose:
            print(nalign)
        f = open(tf.name, 'w')
        f.write(nalign)
        #f.write(self.rf + '\n')
        f.close()
        return RNAalignment(tf.name)

    def __add__(self, rna_seq):
        self.seqs.append(rna_seq)
        self.reload_alignment()

    def write(self, fn, verbose=False):
        """Write the alignment to a file"""
        if verbose:
            print('Save to ', fn)
        with open(fn, 'w') as f:
            f.write('# STOCKHOLM 1.0\n')
            shift = max([len(x) for x in [s.id for s in self.seqs] + ['#=GC=SS_cons']])
            for s in self:
                f.write(s.id.ljust(shift + 2, ' ') + s.seq + '\n')
            f.write('#=GC SS_cons'.ljust(shift + 2, ' ') + self.ss_cons + '\n')
            f.write('//')

    def copy_ss_cons_to_all(self, verbose=False):
        for s in self.io:
            if verbose:
                self.ss_cons
                self.io[0].seq
            try:
                s.letter_annotations['secondary_structure'] = self.ss_cons
            except TypeError:
                raise Exception(
                    'Please check if all your sequences and ss lines are of the same length!')
            s.ss = self.ss_cons
            s.ss_clean = self.get_clean_ss(s.ss)
            s.seq_nogaps = str(s.seq).replace('-', '')
            s.ss_nogaps = self.get_ss_remove_gaps(s.seq, s.ss_clean)

    def copy_ss_cons_to_all_editing_sequence(self, seq_id, before, after):
        """Change a sequence's sec structure.

        :param seq_id: string, sequence id to change, eg: ``AE009948.1/1094322-1094400``
        :param before: string, character to change from, eg: ``,``
        :param after: string, character to change to, eg: ``.``

        .. warning:: before and after has to be one character long
        """
        for s in self.io:
            if s.id == seq_id:
                s.letter_annotations['secondary_structure'] = self.ss_cons.replace(before, after)
            else:
                s.letter_annotations['secondary_structure'] = self.ss_cons

    def get_ss_remove_gaps(self, seq, ss):
        """
        :param seq: string, sequence
        :param ss: string, ss

        UAU-AACAUAUAAUUUUGACAAUAUGG-GUCAUAA-GUUUCUACCGGAAUACC--GUAAAUAUUCU---GACUAUG-UAUA-
        (((.(.((((,,,(((((((_______.))))))).,,,,,,,,(((((((__.._____))))))...),,)))).)))).
        """
        nss = ''
        for i, j in zip(seq, ss):
            if i != '-':
                nss += j
        return nss

    def plot(self, plot_fn='rchie.png'):
        return RChie().plot_cov(self.io.format("fasta"), self.ss_cons_std, plot_fn)

    def get_ss_cons_pk(self):
        """
        :return: SS_cons_pk line or None if there is now SS_cons_pk:"""
        ss_cons_pk = ''
        for l in self.lines:
            if l.startswith('#=GC SS_cons_pk'):
                ss_cons_pk += l.replace('#=GC SS_cons_pk', '').strip()
        return ss_cons_pk

    def get_ss_cons(self):
        """
        :return: SS_cons_pk line or None if there is now SS_cons_pk.
        """
        ss_cons = ''
        for l in self.lines:
            if '#=GC SS_cons' in l and '#=GC SS_cons_pk' not in l:
                ss_cons += l.replace('#=GC SS_cons', '').strip()
        return ss_cons

    @property
    def ss_cons_with_pk(self):
        """go over ss_cons and overwrite bp is there is pk (ss_cons_pk)

        ss_cons:         (((.(.((((,,,(((((((_______.))))))).,,,,,,,,(((((((__.._____))))))...),,)))).)))).
        ss_cons_pk:      .......................[[...............................]]........................
        ss_cons_with_pk: (((.(.((((,,,(((((((___[[__.))))))).,,,,,,,,(((((((__.._]]__))))))...),,)))).)))).

        "return ss_cons_with_pk: string, e.g. (((.(.((((,,,(((((((___[[__.))))"""
        if self.ss_cons_pk:
            ss_cons_with_pk = ''
            for i, (s, p) in enumerate(zip(self.ss_cons, self.ss_cons_pk)):
                if p != '.':
                    ss_cons_with_pk += p
                else:
                    ss_cons_with_pk += s
            return ss_cons_with_pk
        else:
            return self.ss_cons

    @property
    def ss_cons_with_pk_std(self):
        if self.ss_cons_pk:
            ss_cons_with_pk_std = ''
            for i, (s, p) in enumerate(zip(get_rfam_ss_notat_to_dot_bracket_notat(self.ss_cons), self.ss_cons_pk)):
                if p != '.':
                    ss_cons_with_pk_std += p
                else:
                    ss_cons_with_pk_std += s
            return ss_cons_with_pk_std
        else:
            return self.ss_cons

    def get_gc_rf(self):
        """Return (str) ``#=GC RF`` or '' if this line is not in the alignment.
        """
        for l in self.lines:
            if l.startswith('#=GC RF'):
                return l.replace('#=GC RF', '').replace('_cons', '').strip()
        else:
            return ''
            # raise RNAalignmentError('There is on #=GC RF in the alignment!')

    def get_shift_seq_in_align(self):
        """RF_cons vs '#=GC RF' ???"""
        for l in self.lines:
            if l.startswith('#=GC RF'):
                # #=GC RF                        .g.gc.a
                l = l.replace('#=GC RF', '')
                c = 7  # 12 # len of '#=GC RF'
                #                         .g.gc.a
                for i in l:
                    if i == ' ':
                        c += 1
                self.shift = c
                return c

    def map_seq_on_seq(self, seq_id, seq_id_target, resis, v=True):
        """
        :param seq_id: seq_id, 'AAML04000013.1/228868-228953'
        :param seq_id_target: seq_id of target, 'CP000721.1/2204691-2204778'
        :param resis: list resis, [5,6]

        map::

            [4, 5, 6]
            UAU-A
            UAU-AA
            UAU-AAC
            [5, 6, 7]
            CAC-U
            CAC-U-
            CAC-U-U
            [4, None, 5]

        """
        # get number for first seq
        print(resis)
        nresis = []
        for s in self.io:
            if s.id.strip() == seq_id.strip():
                for i in resis:
                    print(s.seq[:i + 1])  # UAU-A
                    nresis.append(i + s.seq[:i].count('-'))
        print(nresis)

        print(self.map_seq_on_align(seq_id_target, nresis))
        return

        resis_target = []
        for s in self.io:
            if s.id.strip() == seq_id_target.strip():
                for i in nresis:
                    if v:
                        print(s.seq[:i])
                    if s.seq[i - 1] == '-':
                        resis_target.append(None)
                    else:
                        resis_target.append(i - s.seq[:i].count('-'))
        return resis_target

    def map_seq_on_align(self, seq_id, resis, v=True):
        """
        :param seqid: seq_id, 'CP000721.1/2204691-2204775'
        :param resis: list resis, [5,6]

        maps::

            [5, 6, 8]
            CAC-U
            CAC-U-
            CAC-U-UA
            [4, None, 6]

        """
        if v:
            print(resis)
        nresis = []
        for s in self.io:
            if s.id.strip() == seq_id.strip():
                for i in resis:
                    if v:
                        print(s.seq[:i])
                    if s.seq[i - 1] == '-':
                        nresis.append(None)
                    else:
                        nresis.append(i - s.seq[:i].count('-'))
        return nresis

    def head(self):
        return '\n'.join(self.lines[:5])

    def tail(self):
        return '\n'.join(self.lines[-5:])

    def describe(self):
        """Describe the alignment.

           > print(a.describe())
           SingleLetterAlphabet() alignment with 13 rows and 82 columns

        """
        return str(self.io).split('\n')[0]

    def remove_empty_columns(self, verbose=False):
        """Remove empty columns in place.

        Example::

            >>> a = RNAalignment("test_data/zmp.stk")
            >>> print(a)
            SingleLetterAlphabet() alignment with 6 rows and 319 columns
            ---ACCUUGCGCGACUGGCGAAUCC-------------------...AAU CP001644.1/756294-756165
            --GCUCUCGCGCGACUGGCGACUUUG------------------...GAA CU234118.1/352539-352459
            UGAGUUUUCUGCGACUGACGGAUUAU------------------...CUG BAAV01000055.1/2897-2982
            GCCCGUUCGCGUGACUGGCGCUAGU-------------------...CGA CP000927.1/5164264-5164343
            -----GGGUCGUGACUGGCGAACA--------------------...--- zmp
            UCACCCCUGCGUGACUGGCGAUA---------------------...GUU AP009385.1/718103-718202
            >>> a.remove_empty_columns()
            >>> print(a)
            SingleLetterAlphabet() alignment with 6 rows and 138 columns
            ---ACCUUGCGCGACUGGCGAAUCC-UGAAGCUGCUUUG-AGCG...AAU CP001644.1/756294-756165
            --GCUCUCGCGCGACUGGCGACUUUG------------------...GAA CU234118.1/352539-352459
            UGAGUUUUCUGCGACUGACGGAUUAU------------------...CUG BAAV01000055.1/2897-2982
            GCCCGUUCGCGUGACUGGCGCUAGU-------------------...CGA CP000927.1/5164264-5164343
            -----GGGUCGUGACUGGCGAACA--------G-----------...--- zmp
            UCACCCCUGCGUGACUGGCGAUA--------GAACCCUCGGGUU...GUU AP009385.1/718103-718202

        go over all seq
        modifes self.nss_cons"""
        cols_to_rm = []

        # get only seqs
        for i in range(len(self[0])):
            gap = True
            for seq in self:
                if seq.seq[i] != '-':
                    gap = False
            if gap:
                cols_to_rm.append(i)

        # remove from sequences
        for s in self:
            s.remove_columns(cols_to_rm)

        # update io # hack #
        tmpfn = tempfile.NamedTemporaryFile(delete=False).name
        self.write(tmpfn, verbose=False)
        self.io = AlignIO.read(tmpfn, "stockholm")

        # nss_cons update
        nss_cons = ''
        for i, s in enumerate(self.ss_cons):
            if i not in cols_to_rm:
                nss_cons += s
        self.ss_cons = nss_cons

    def format_annotation(self, t):
        return self.shift * ' ' + t

    def find_core(self, ids=None):
        """Find common core for ids.

        .. image:: ../pngs/find_core.png
        Fig. By core, we understand columns that have all homologous residues. The core is here marked by `x`.

        :param id: list, ids of seq in the alignment to use
        """
        if not ids:
            ids = []
            for s in self.io:
                ids.append(s.id)

        xx = list(range(0, len(self.io[0])))
        for i in range(0, len(self.io[0])):  # if . don't use it
            for s in self.io:
                # print s.id
                if s.id in ids:
                    if s.seq[i] == '-':
                        xx[i] = '-'
                        break
                    else:
                        xx[i] = 'x'

        return ''.join(xx)
        shift = self.get_shift_seq_in_align()

        fnlist = open(self.fn).read().strip().split('\n')
        fnlist.insert(-2, 'x' + ' ' * (shift - 1) + ''.join(xx))
        # print fnlist
        for l in fnlist:
            print(l)

    def find_seq(self, seq, verbose=False):
        """Find seq (also subsequences) and reverse in the alignment.

        Args:

          seq (str): seq is upper()
          verbose (bool): be verbose

        ::

            seq = "ggaucgcugaacccgaaaggggcgggggacccagaaauggggcgaaucucuuccgaaaggaagaguaggguuacuccuucgacccgagcccgucagcuaaccucgcaagcguccgaaggagaauc"
            hit = a.find_seq(seq, verbose=False)
            ggaucgcugaacccgaaaggggcgggggacccagaaauggggcgaaucucuuccgaaaggaagaguaggguuacuccuucgacccgagcccgucagcuaaccucgcaagcguccgaaggagaauc
            Match: AL939120.1/174742-174619
            ID: AL939120.1/174742-174619
            Name: AL939120.1
            Description: AL939120.1/174742-174619
            Number of features: 0
            /start=174742
            /end=174619
            /accession=AL939120.1
            Per letter annotation for: secondary_structure
            Seq('CCAGGUAAGUCGCC-G-C--ACCG---------------GUCA-----------...GGA', SingleLetterAlphabet())
            GGAUCGCUGAACCCGAAAGGGGCGGGGGACCCAGAAAUGGGGCGAAUCUCUUCCGAAAGGAAGAGUAGGGUUACUCCUUCGACCCGAGCCCGUCAGCUAACCUCGCAAGCGUCCGAAGGAGAAUC

        """
        seq = seq.replace('-', '').upper()
        for s in self.io:
            seq_str = str(s.seq).replace('-', '').upper()
            if verbose:
                print(seq_str)
            if seq_str.find(seq) > -1 or seq.find(seq_str) > -1:
                print('Match:', s.id)
                print(s)
                print(seq)
        return s
        print('Not found')

    def find_seq_exact(self, seq, verbose=False):
        """Find seq (also subsequences) and reverse in the alignment.

        :param seq: string, seq, seq is upper()
        :param verbose: boolean, be verbose or not
        """
        seq = seq.replace('-', '').upper()
        for s in self.io:
            seq_str = str(s.seq).replace('-', '').upper()
            if verbose:
                print(seq_str)
            if seq_str == seq:
                print('Match:', s.id)
                print(s)
                print(seq)
                return s
        print('Not found')

    def get_clean_ss(self, ss):
        nss = ''
        for s in ss:
            nss += get_rfam_ss_notat_to_dot_bracket_notat(s)
        return nss

    def get_seq_ss(self, seq_id):  # seq,ss):
        s = self.get_seq(seq_id).seq
        # print seq, ss
        # new
        nseq = ''
        nss = ''
        for i, j in zip(seq, ss):
            if i != '-':  # gap
                # print i,j
                nseq += i
                nss += get_rfam_ss_notat_to_dot_bracket_notat(j)
        return nseq.strip(), nss.strip()

    def get_seq(self, seq_id):
        for s in self.io:
            if s.id == seq_id:
                return s
        raise Exception('Seq not found')

    def get_seq_with_name(self, seq_name):
        for s in self.io:
            if s.name == seq_name:
                return s
        raise Exception('Seq not found')

    def align_seq(self, seq):
        """Align seq to the alignment.

        Using self.rf.

        Args:

            seq (str): sequence, e.g. ``-GGAGAGUA-GAUGAUUCGCGUUAAGUGUGUGUGA-AUGGGAUGUC...``

        Returns:

           str: seq that can be inserted into alignemnt, ``.-.GG.AGAGUA-GAUGAUUCGCGUUA`` ! . -> -
        """
        seq = list(seq)
        seq.reverse()
        nseq = ''
        for n in self.rf:  # n nuclotide
            if n != '.':
                try:
                    j = seq.pop()
                except:
                    j = '.'
                nseq += j
            if n == '.':
                nseq += '.'  # + j
        return nseq.replace('.', '-')

    def __repr__(self):
        return (str(self.io))

    def trimmed_rf_and_ss(self):
        """Remove from RF and SS gaps.

        Returns:

          (str,str): trf, tss - new RF and SS

        """
        trf = ''
        tss = ''
        for r, s in zip(self.rf, self.ss_cons_std):
            if r not in ['-', '.']:
                trf += r
                tss += s
        return trf, tss

    def get_distances(self):
        """Get distances (seq identity) all-vs-all.

        With BioPython.

        blastn: ``Bad alphabet 'U' in sequence 'AE008922.1/409481-409568' at position '7'`` only for DNA?

        read more (also about matrix at <http://biopython.org/wiki/Phylo> and
        HTTP://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceCalculator-class.html
        """
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(self.io)
        return dm

    def get_the_closest_seq_to_ref_seq(self, verbose=False):
        """

        Example::

            >>> a = RNAalignment("test_data/RF02221.stockholm.sto")
            >>> a.get_the_closest_seq_to_ref_seq()
            AF421314.1/431-344

        """
        self + RNASeq('ConSeq', self.rf, '')
        dist = self.get_distances()
        distConSeq = dist['ConSeq'][:-1]  # to remove ending 0, bc of distance to itself
        minimal = min(distConSeq)

        index = dist['ConSeq'].index(minimal)
        #id = dist.names[index]
        if verbose:
            print('dist:\n', str(dist))
            print('distConSeq:', dist['ConSeq'])
        return self[index]


class CMAlign():
    """CMAalign class around cmalign (of Inferal).

    cmalign - aligns  the RNA sequences in <seqfile> to the covariance model
    (CM) in <cmfile>.  The new alignment is output to stdout  in  Stockholm
    format.

    Example::

        cma = ra.CMAlign()
        cma.run_cmalign("ade_seq.fa", "RF00167.cm")
        seq = cma.get_seq()
        print 'cma hit  ', seq
        print 'seq      ', a.align_seq(seq)
        print 'a.rf     ', a.rf

        cmd cmalign -g RF00167.cm ade_seq.fa

        # STOCKHOLM 1.0
        #=GF AU Infernal 1.1.2

        ade          ----------------CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAG-CCUUAAA-CUCUUGAUUAUGAAGUGA------------
        #=GR ade PP  ................99*********************************************.*******.***************999............
        #=GC SS_cons :::::::::::::::::((((((((,,,<<<<<<<_______>>>>>>>,,,,,,,,<<<<<<<_______>>>>>>>,,))))))))::::::::::::::
        #=GC RF      aaaaaauaaaaaaaauucccuCgUAUAAucccgggAAUAUGGcccgggaGUUUCUACCaggcagCCGUAAAcugccuGACUAcGagggaaauuuuuuuuuuu
        //
        cma hit   ----------------CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAG-CCUUAAA-CUCUUGAUUAUGAAGUGA------------
        seq       ----------------CGCU-U-CAUAUAAUCCUAAUGAUAUGG-UUUGGGA-GUUUCUACCAAGAG-CC--UUAAA-CUCUU---GAUUAUG-AAGUGA-------------
        a.rf      aaaaaauaaaaaaaauuccc.u.CgUAUAAucccgggAAUAUGG.cccggga.GUUUCUACCaggcagCC..GUAAAcugccu...GACUAcG.agggaaauuuuuuuuuuu.

    Install http://eddylab.org/infernal/

    Cite: Nawrocki and S. R. Eddy, Infernal 1.1: 100-fold faster RNA homology searches, Bioinformatics 29:2933-2935 (2013). """

    def __init__(self, outputfn=None):
        """Use run_cmalign or load cmalign output from a file"""
        if outputfn:
            self.output = open(outputfn).read().strip().split('\n')

    def run_cmalign(self, seq, cm, verbose=True):
        """Run cmalign and process the result.

        :param seq: seq string
        :param cm: cm fn

        Run::

            $ cmalign RF01831.cm 4lvv.seq
            # STOCKHOLM 1.0
            #=GF AU Infernal 1.1.2

            4lvv         -GGAGAGUA-GAUGAUUCGCGUUAAGUGUGUGUGA-AUGGGAUGUCG-UCACACAACGAAGC---GAGA---GCGCGGUGAAUCAUU-GCAUCCGCUCCA
            #=GR 4lvv PP .********.******************9999998.***********.8999999******8...5555...8**************.************
            #=GC SS_cons (((((----(((((((((((,,,,,<<-<<<<<<<<___________>>>>>>>>>>,,,<<<<______>>>>,,,)))))))))))-------)))))
            #=GC RF      ggcaGAGUAGggugccgugcGUuAAGUGccggcgggAcGGGgaGUUGcccgccggACGAAgggcaaaauugcccGCGguacggcaccCGCAUcCgCugcc
            //

        .. warning :: requires cmalign to be set in your shell
        """
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.name += '.seq'

        with open(tf.name, 'w') as f:
            f.write('>target\n')
            f.write(seq + '\n')

        cmd = 'cmalign -g ' + cm + ' ' + tf.name  # global
        if verbose:
            print('cmd' + cmd)
        o = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout = o.stdout.read().strip()
        stderr = o.stderr.read().strip()
        if verbose:
            print(stdout)
        self.output = stdout.split('\n')

    def get_gc_rf(self):
        """Get ``#=GC RF``.

        :var self.output: string
        """
        for l in self.output:
            if l.startswith('#=GC RF'):
                return l.replace('#=GC RF', '').strip()

    def get_seq(self):
        """
        :var self.output: output of cmalign, string
        """
        for l in self.output:
            if l.strip():
                if not l.startswith('#'):
                    #  4lvv         -GGAGAGUA-GAUGAU
                    return l.split()[1].strip()


def clean_seq_and_ss(seq, ss):
    nseq = ''
    nss = ''
    for i, j in zip(seq, ss):
        if i != '-':  # gap
            # print i,j
            nseq += i
            nss += get_rfam_ss_notat_to_dot_bracket_notat_per_char(j)
    return nseq.strip(), nss.strip()


def get_rfam_ss_notat_to_dot_bracket_notat_per_char(c):
    """Take (c)haracter and standardize ss (including pks in letter (AAaa) notation).

    .. warning:: DD DD will be treated as BB BB (<< >>) (it might be wrong!)"""
    if c in [',', '_', ':', '-']:
        return '.'
    if c == '<':
        return '('
    if c == '{':
        return '('
    if c == '}':
        return ')'
    if c == ']':
        return ')'
    if c == '[':
        return '('
    if c == '>':
        return ')'
    if c == '>':
        return ')'
    if c == 'A':
        return '['
    if c == 'a':
        return ']'
    if c == 'B':
        return '<'
    if c == 'b':
        return '>'
    if c == 'C':
        return '{'
    if c == 'c':
        return '}'
    # !!!!!!!!!!!!!!!!!!!!!!!!!!
    if c == 'D':
        return '<'
    if c == 'd':
        return '>'
    # !!!!!!!!!!!!!!!!!!!!!!!!!!
    return c


def get_rfam_ss_notat_to_dot_bracket_notat(ss):
    """Change all <<>> to "standard" dot bracket notation.

    Works also with pseudknots AA BB CC etc."""
    nss = ''
    for s in ss:
        ns = get_rfam_ss_notat_to_dot_bracket_notat_per_char(s)
        nss += ns
    return nss


def fasta2stokholm(fn):
    """Take a gapped fasta file and return an RNAalignment object.

    A fasta file should look like this::

        >AE009948.1/1094322-1094400
        UAC-U-UAUUUAUGCUGAGGAU--UGG--CUUAGC-GUCUCUACAAGACA-CC--GU-AA-UGUCU---AACAAUA-AGUA-
        ...
        >CP000721.1/2204691-2204778
        CAC-U-UAUAUAAAUCUGAUAAUAAGG-GUCGGAU-GUUUCUACCAAGCUACC--GUAAAUUGCUAUAGGACUAUA-AGUG-
        >SS_cons
        (((.(.((((...(((((((........))))))).........(((((((.........))))))...)..)))).)))).
        >SS_cons_pk
        .........................[[........................]].............................

    ``SS_cons_pk`` in optionally and is an extra line used to define a pseudoknot. You
    can also define second psedoknot as ``<<<...>>>`` and the third one with ``{{{ }}}``.

    :param fn: file
    :return: RNAalignment object
    """
    seqs = []
    s = None
    for l in open(fn):
        if l.startswith('>'):
            if s:
                seqs.append(s)
            s = RNASeq(l.replace('>', '').strip(), '')
        else:
            s.seq += l.strip()

    txt = ''
    for l in open(fn):
        if l.startswith('>'):
            id = '\n' + l.replace('>', '\n').strip() + ' '
            id = re.sub('ss_cons_pk', '#=GC SS_cons_pk', id, flags=re.IGNORECASE)
            id = re.sub('ss_cons', '#=GC SS_cons', id, flags=re.IGNORECASE)
            txt += id
        else:
            txt += l.strip()
    txt = txt.strip() + '\n'  # clean upfront \n and add tailing \n

    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.name += '.stk'
    with open(tf.name, 'w') as f:
        f.write('# STOCKHOLM 1.0\n')
        f.write(txt)
        f.write('//\n')
    return RNAalignment(tf.name)


def fetch_stokholm(rfam_acc, dpath=None):
    """Fetch Stokholm file from Rfam.

    :param rfam_acc: str, Rfam accession number, eg. RF00028
    :param dpath: str or None, if None saves to current location, otherwise save to dpath folder

    .. warning :: urllib3 is required, pip install urllib3
    """
    import urllib3
    http = urllib3.PoolManager()
    # try:
    print(dpath)
    response = http.request('GET', url='http://rfam.xfam.org/family/' + rfam_acc.lower() +
                            '/alignment?acc=' + rfam_acc.lower() + '&format=stockholm&download=1')
    # except urllib3.HTTPError:
    #    raise Exception('The PDB does not exists: ' + pdb_id)
    txt = response.data
    if dpath is None:
        npath = rfam_acc + '.stk'
    else:
        npath = dpath + os.sep + rfam_acc + '.stk'
    print('downloading...' + npath)
    with open(npath, 'wb') as f:
        f.write(txt)
    print('ok')
    return rfam_acc + '.stk'


# def test_seq_multilines_alignment()
def test_alignment_with_pk():
    a = RNAalignment('test_data/RF00167.stockholm.sto')
    print(a.tail())
    print(a.ss_cons)
    print(a.ss_cons_pk)
    print(a.ss_cons_with_pk)

    print(a[[1, 2]])


# main
if __name__ == '__main__':
    ## a = RNAalignment('test_data/RF00167.stockholm.sto')
    # print a.get_shift_seq_in_align()
    # print a.map_seq_on_align('CP000721.1/2204691-2204778', [5,6,8])
    # print a.map_seq_on_seq('AAML04000013.1/228868-228953', 'CP000721.1/2204691-2204778', [4,5,6])

    # print(record.letter_annotations['secondary_structure'])
    ## seq = a.get_seq('AAML04000013.1/228868-228953')
    # print seq.seq
    # print seq.ss
    # print a.ss_cons

    # print 'x'
    # for s in a.io:
    # print s.seq
    # print s.ss_clean #letter_annotations['secondary_structure']

    # for s in a.io:
    # print s.seq_nogaps
    # print s.ss_nogaps

    # a.write('test_output/out.sto')

    # a = RNAalignment("/home/magnus/work/rna-evo/rp12/seq/RF00379.stockholm.stk")##_rm85.stk')
    # print a.get_seq_ss("AL939120.1/174742-174619")
    # subset = a.subset(["AL939120.1/174742-174619",
    # "rp12bsubtilis",
    # "AAWL01000001.1/57092-56953",
    # "CP000612.1/87012-87130",
    # "BA000028.3/1706087-1706245"])
    # % cat /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpTzPenx
    # for s in subset:
    # print s.seq
    # subset.remove_empty_columns()
    # subset.write('/home/magnus/out2.stk')
    # for s in subset:
    # print s.seq

    # print 'subset.ss_cons_std:', subset.ss_cons_std

    # a = RNAalignment("/home/magnus/work/rna-evo/rp12/seq/RF00379.stockholm.stk")##_rm85.stk')
    # a.ss_cons = "::::::::::{{{{,,,<<<_..."
    # print a.ss_cons
    # print a.ss_cons_std
    # print get_rfam_ss_notat_to_dot_bracket_notat(subset.ss_cons)

    ## a.ss_cons_std = 'test of setter'
    # print a._ss_cons_std

    # pass
    # slice

    ## slices = a.cols[0:10]
    # for s in a:
    # print s.seq
    # add pk

    # for s in slices:
    # print s, s.seq, s.ss

    #a = fasta2stokholm('test_output/ade_gapped.fa')
    # print a

    #a = RNAalignment('test_data/RF00001.blocked.stk')
    # print a
    # print a.ss_cons
    # a.plot('rchie.png')

    # a = RNAalignment("test_data/RF02221.stockholm.sto")
    ## a = RNAalignment('test_data/test_data/RF02221.stockholm.sto')
    ## a + RNASeq('ConSeq', '-A-GU-AGAGUA-GGUCUUAUACGUAA-----------------AGUG-UCAUCGGA-U-GGGGAGACUUCCGGUGAACGAA-G-G-----------------------------GUUA---------------------------CCGCGUUAUAUGAC-C-GCUUCCG-CUA-C-U-','')
    ## dist = a.get_distances()
    # distConSeq = dist['ConSeq'][:-1]  # to remove ending 0, bc of distance to itself
    ## minimal = min(distConSeq)
    ## index = dist['ConSeq'].index(minimal)
    # print(dist.names[index])

    a = RNAalignment("test_data/dist_test2.stk")
    rep = a.get_the_closest_seq_to_ref_seq()
    rep.remove_gaps()
    print(rep)
    print(rep.ss)
    print(rep.seq)

    # a.write('tmp.stk')
    # s.remove_gaps()
    # print(s.seq)
    # print(s.ss)

    import doctest
    doctest.testmod()
