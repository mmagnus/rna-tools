
## History

170608 Add `--get_ss` (secondary structure) using x3dna.

170518 Edit in place [experimental, only for `get_rnapuzzle_ready`] `rna_pdb_tools.py --rpr 7_Das_7_rpr.pdb --inplace`. (2) get a structure in org-mode format <sick!>

170517 Fix #37 mis-align atom names after rpr-ing bug

170515 Fix fixing missing O2'

170404	`rna_simrna_extract.py -t template.pdb -f *05.trafl -c -n 1 # extract only the first model`

170331 rna-pdb-tools meets Emacs!

![](docs/pngs/rpt_emacs.png)

170325 Seq: secondary structure prediction with constraints

    >>> seq = Seq("CCCCUUUUGGGG")
    >>> seq.name = 'RNA03'
    >>> print(seq.predict_ss("RNAfold", constraints="((((....))))"))
    >RNA03
    CCCCUUUUGGGG
    ((((....)))) ( -6.40)

170324 Starting converting to Python3, fetch_align by Pietro

170320 `rna_cartoon` in PyMOL

![](docs/pngs/rna_cartoon_small.png)

170319 Add clanstix (move it from its own GitHub repository).

170315 SimRNA_trajectory:
  - get len of frame, and trajectory
  - warn about broken frame
  - `only_first_frame` to get only the first frame

170311 Get seq (v2) gets segments of chains with correct numbering

	> 6_solution_0 A:1-19 26-113 117-172
	GGCGGCAGGUGCUCCCGACGUCGGGAGUUAAAAGGGA

170308 Add fixing missing atoms of bases, and O2'

... many things! :-)

~2011 Prelimiary version as rnastruc, yapdb_parser etc.
