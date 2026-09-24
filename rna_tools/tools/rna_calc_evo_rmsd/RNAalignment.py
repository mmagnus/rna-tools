#!/usr/bin/env python

"""RNAalignment

Example::

    # STOCKHOLM 1.0

    AACY023581040                --CGUUGGACU------AAA--------AGUCGGAAGUAAGC-----AAU-C------GCUGAAGCAACGC---
    AJ630128                     AUCGUUCAUUCGCUAUUCGCA-AAUAGCGAACGCAA--AAG------CCG-A-------CUGAAGGAACGGGAC
    target                       --CGUUGACCCAG----GAAA-----CUGGGCGGAAGUAAGGCCCAUUGCACUCCGGGCCUGAAGCAACGCG--
    #=GC SS_cons                 ::(((((,<<<<<<<.._____..>>>>>>>,,,,,,,,<<<<...._____.....>>>>,,,,)))))::::
    x                            --xxxxxxxxx-----------------xxxxxxxx--xxx------------------xxxxxxxxxxx----
    #=GC RF                      AUCGUUCAuCucccc..uuuuu..ggggaGaCGGAAGUAGGca....auaaa.....ugCCGAAGGAACGCguu
//

x line is used to pick resides to calculate RMSD.

x line could be renamed to EvoClust.  When no selector sequence is available
the ``#=GC RF`` reference annotation is used instead, and if neither is
present the selector is derived from columns without gaps."""
from __future__ import print_function
import warnings
warnings.filterwarnings("ignore")

from Bio import AlignIO

class RNAalignment:
    """RNAalignemnt"""

    def __init__(self, fn, verbose=False):
        """Load the alignment in the Stockholm format using biopython"""
        self.verbose = verbose
        self.alignment = AlignIO.read(open(fn), "stockholm")
        if self.verbose:
            print('Alignment records:')
            for record in self.alignment:
                print(' ', record.id, str(record.seq))
        self._selection_track, self._selection_source = self._get_selection_track()
        if self.verbose:
            print('Selector source:', self._selection_source)
            print('Selector (x-line):', self._selection_track)

    def _get_selection_track(self):
        """Return a string that marks positions to be used for RMSD calculation.

        Historically an "x" (or EvoClust) sequence was appended to the alignment
        and the x characters were used as selectors.  Modern Infernal-generated
        Stockholm files tend to provide a ``#=GC RF`` line instead.  Support both
        by checking for an explicit sequence first and falling back to the
        reference annotation.
        """

        def is_selector_record(record_id):
            rid = record_id.strip().lower()
            return rid in ("x", "evoclust")

        for record in reversed(self.alignment):
            if is_selector_record(record.id):
                return str(record.seq), "sequence"

        column_annotations = getattr(self.alignment, "column_annotations", {})
        reference = column_annotations.get("reference_annotation")
        if reference:
            return reference, "reference_annotation"

        auto_selector = self._build_gapless_selector()
        if auto_selector:
            return auto_selector, "gapless_columns"

        raise Exception(
            "Alignment must contain an EvoClust/x line, a #=GC RF reference annotation, or columns without gaps"
        )

    def get_column_positions(self, seqid, offset=0):
        """Get the residue position (1-based, + offset) for every selected column.

        Returns a list with one entry per selected column of the alignment; the
        entry is ``None`` when the sequence has a gap in that column. Lists
        returned for different sequences are therefore column-wise comparable,
        which is what is needed to pair residues between two structures.

        If seqid is not found in the alignment, raise an exception."""
        x = self._selection_track  # ---(((((((----xxxxx-- or RF annotation
        for record in self.alignment:
            if record.id.strip().lower() in ("x", "evoclust"):
                continue
            if record.id == seqid.strip():
                positions = []
                spos = 0
                for xi, si in zip(x, str(record.seq)):
                    is_gap = si in ('-', '.')
                    if not is_gap:
                        spos += 1
                    if self._include_column(xi):
                        positions.append(None if is_gap else spos + offset)
                return positions
        raise Exception('Seq not found in the alignment: %s' % seqid)

    def get_range(self, seqid, offset=0, verbose=None):
        """Get a list of positions for selected residues based on the selector line.

        Columns where the sequence has a gap are skipped.

        If seqis not found in the alignment, raise an exception, like ::

            Exception: Seq not found in the alignment: 'CP000879.1/21644622164546
        """
        if verbose is None:
            verbose = self.verbose
        x_range = [p for p in self.get_column_positions(seqid, offset) if p is not None]
        if verbose:
            print(' Selected residues for %s: %s' % (seqid, x_range))
        if not x_range:
            raise Exception('Seq not found or selector line is malformed')
        return x_range

    def get_paired_positions(self, seqid1, seqid2, offset=0, verbose=None):
        """Get residue positions of two sequences paired by selected alignment columns.

        Only columns where both sequences have a residue are kept.

        :returns: (positions1, positions2), two lists of the same length"""
        if verbose is None:
            verbose = self.verbose
        pairs = [(a, b) for a, b in zip(self.get_column_positions(seqid1, offset),
                                         self.get_column_positions(seqid2, offset))
                 if a is not None and b is not None]
        if not pairs:
            raise Exception('No selected columns shared by %s and %s' % (seqid1, seqid2))
        positions1 = [a for a, b in pairs]
        positions2 = [b for a, b in pairs]
        if verbose:
            print(' Paired residues %s <-> %s: %s' % (seqid1, seqid2,
                                                    ' '.join('%s:%s' % p for p in pairs)))
        return positions1, positions2

    def _include_column(self, selector_char):
        """Return True when the selector character marks a column to use."""
        if not selector_char:
            return False
        if selector_char in ('-', '.', ' '):
            return False
        return True

    def _build_gapless_selector(self):
        """Build a selector string that keeps columns without gaps in any sequence."""
        sequences = [str(record.seq) for record in self.alignment
                     if record.id.strip().lower() not in ("x", "evoclust")]
        if not sequences:
            return None

        selector_chars = []
        for column in zip(*sequences):
            if all(base not in ('-', '.') for base in column):
                selector_chars.append('x')
            else:
                selector_chars.append('-')
        selector = ''.join(selector_chars)
        if 'x' not in selector:
            return None
        if self.verbose:
            print(' Inferred selector from gapless columns.')
        return selector


if __name__ == '__main__':
    ra = RNAalignment('test_data/rp13finalX.stk')
    print(len(ra.get_range('rp13target', offset=0)))
    print(len(ra.get_range('cp0016', offset=0)))
    print(len(ra.get_range('NZ_ABBD01000528', offset=0)))

    print(ra.get_range('rp13target', offset=0))
    print(ra.get_range('cp0016', offset=0))
    print(ra.get_range('NZ_ABBD01000528', offset=0))
