PyMOL: color by conservation
--------------------------------------

Get an alignment for your protein sequence from <https://toolkit.tuebingen.mpg.de/tools/psiblast>

Download the alignment and open the alignment in JalView and save Export Annotation -> To CSV

![](docs/jalview.png)

The file should like this:

![](docs/consensus-file.png)
	 
Run `pymol_color_by_conserv.py` on the alignment and the consensus file:

    (base) [mm] pymol_color_by_conserv$ git:(master) ✗ ./pymol_color_by_conserv.py -h
    usage: pymol_color_by_conserv.py [-h] [-v] conserv alignment chain

    positional arguments:
      conserv        consigns file
      alignment
      chain

    optional arguments:
      -h, --help     show this help message and exit
      -v, --verbose  be verbose

example:

    python pymol_color_by_conserv.py data/cwc15/cwc15_consensus.csv data/cwc15/psiblast_8183496.aln A

Load your protein in PyMOL and load the script (pml file) to color residues according to conservation:

    # in PyMOL
    @ <path to a file> cwc15.pml

![](docs/cwc15.png)

Enjoy!