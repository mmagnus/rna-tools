rna-pdb-tools
===============================================

rna_pdb_tools - a swiss army knife to manipulation of RNA pdb structures

Remove atoms with XYZ equals 0::

    (base) ➜  MiniRoseTTA git:(main) ✗ cat  4GXY_min_at0_chemicals.pdb
    ATOM      1  OP1   G A   1      50.150  76.113  39.198  1.00  0.00
    ATOM      2  P     G A   1      50.001  77.254  40.137  1.00  0.00
    ATOM      3  OP2   G A   1      48.880  77.258  41.111  1.00  0.00
    ATOM      4  O5'   G A   1      51.362  77.417  40.948  1.00  0.00
    ATOM      5  C5'   G A   1       0.000   0.000   0.000  1.00  0.00
    ATOM      6  C4'   G A   1       0.000   0.000   0.000  1.00  0.00
    ATOM      7  O4'   G A   1       0.000   0.000   0.000  1.00  0.00
    ATOM      8  C3'   G A   1       0.000   0.000   0.000  1.00  0.00

to get::

    (base) ➜  MiniRoseTTA git:(main) ✗ rna_pdb_tools.py --remove0  4GXY_min_at0_chemicals.pdb
    ATOM      1  OP1   G A   1      50.150  76.113  39.198  1.00  0.00
    ATOM      2  P     G A   1      50.001  77.254  40.137  1.00  0.00
    ATOM      3  OP2   G A   1      48.880  77.258  41.111  1.00  0.00
    ATOM      4  O5'   G A   1      51.362  77.417  40.948  1.00  0.00
    ATOM     35  OP1   C A   2      54.648  73.216  44.394  1.00  0.00
    ATOM     36  P     C A   2      53.712  74.058  43.607  1.00  0.00
    ATOM     37  OP2   C A   2      53.842  74.111  42.128  1.00  0.00
    ATOM     38  O5'   C A   2      52.223  73.613  43.957  1.00  0.00


.. autoprogram:: rna_tools.rna_pdb_tools:get_parser()
   :prog: rna_pdb_tools.py

get RNAPuzzle ready
-----------------------------------------

.. autoclass:: rna_tools.rna_tools_lib.RNAStructure
   :members: get_rnapuzzle_ready

get sequence
-----------------------------------------

Example::

      $ rna_pdb_tools.py --get-seq 5_solution_1.pdb
      > 5_solution_1.pdb A:1-576
      CAUCCGGUAUCCCAAGACAAUCUCGGGUUGGGUUGGGAAGUAUCAUGGCUAAUCACCAUGAUGCAAUCGGGUUGAACACUUAAUUGGGUUAAAACGGUGGGGGACGAUCCCGUAACAUCCGUCCUAACGGCGACAGACUGCACGGCCCUGCCUCAGGUGUGUCCAAUGAACAGUCGUUCCGAAAGGAAG

.. autoclass:: rna_tools.rna_tools_lib.RNAStructure
   :members: get_seq

fetch
-----------------------------------------

Example::

  $ rna_pdb_tools.py --fetch 1xjr
  downloading...1xjr ok

.. autofunction:: rna_tools.rna_tools_lib.fetch

fetch Biological Assembly
-----------------------------------------

Example::

  $ rna_pdb_tools.py --fetch-ba 1xjr
  downloading...1xjr_ba.pdb ok

or over a list of pdb ids in a text file::

  $ cat data/pdb_ids.txt
  1y26
  1fir

  $ while read p; do rna_pdb_tools.py --fetch-ba $p; done < data/pdb_ids.txt
  downloading...1y26_ba.pdb ok
  downloading...1fir_ba.pdb ok

  $ ls *.pdb
  1fir_ba.pdb 1y26_ba.pdb

.. autofunction:: rna_tools.rna_tools_lib.fetch_ba

delete
-----------------------------------------

Examples::

    $ for i in *pdb; do rna_pdb_tools.py --delete A:46-56 $i > ../rpr_rm_loop/$i ; done

go over all files in the current directory, remove a fragment of chain A, residues between 46-56 (including them) and save outputs to in the folder `rpr_rm_loops`.

edit
-----------------------------------------

.. autofunction:: rna_tools.rna_tools_lib.edit_pdb

the library
-----------------------------------------

  .. automodule:: rna_tools.rna_tools_lib
   :members:


PDB Edit Bfactor/Occupancy
------------------------------------------

.. autoprogram:: rna_tools.tools.rna_pdb_edit_occupancy_bfactor.rna_pdb_edit_occupancy_bfactor:get_parser()
   :prog: rna_pdb_edit_occupancy_bfactor.py

.. autofunction:: rna_tools.tools.rna_pdb_edit_occupancy_bfactor.rna_pdb_edit_occupancy_bfactor.edit_occupancy_of_pdb

Add chain to a file
------------------------------------------

.. automodule:: rna_tools.tools.misc.rna_add_chain
   :members:
   :undoc-members:

.. autoprogram:: rna_tools.tools.misc.rna_add_chain:get_parser()
   :prog: rna_add_chain.py

Measure distance between atoms
------------------------------------------

.. autoprogram:: rna_tools.tools.pdbs_measure_atom_dists.pdbs_measure_atom_dists:get_parser()
   :prog: pdbs_measure_atom_dists.py

.. automodule:: rna_tools.tools.pdbs_measure_atom_dists.pdbs_measure_atom_dists
   :members:
   :undoc-members:
