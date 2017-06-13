rna-pdb-tools
========================================

rna-pdb-tools is the main code of this library.

.. argparse::
   :ref: rna_pdb_tools.rna_pdb_tools.get_parser
   :prog: rna_pdb_tools.py

get RNAPuzzle ready
-----------------------------------------

.. autoclass:: rna_pdb_tools.pdb_parser_lib.StrucFile
   :members: get_rnapuzzle_ready

get sequence
-----------------------------------------

Example::

      $ rna_pdb_tools.py --get_seq 5_solution_1.pdb
      > 5_solution_1.pdb A:1-576
      CAUCCGGUAUCCCAAGACAAUCUCGGGUUGGGUUGGGAAGUAUCAUGGCUAAUCACCAUGAUGCAAUCGGGUUGAACACUUAAUUGGGUUAAAACGGUGGGGGACGAUCCCGUAACAUCCGUCCUAACGGCGACAGACUGCACGGCCCUGCCUCAGGUGUGUCCAAUGAACAGUCGUUCCGAAAGGAAG

.. autoclass:: rna_pdb_tools.pdb_parser_lib.StrucFile
   :members: get_seq

fetch
-----------------------------------------

Example::

  $ rna_pdb_tools.py --fetch 1xjr
  downloading...1xjr ok

.. autofunction:: rna_pdb_tools.pdb_parser_lib.fetch

fetch Biological Assembly
-----------------------------------------

Example::

  $ rna_pdb_tools.py --fetch_ba 1xjr
  downloading...1xjr_ba.pdb ok

or over a list of pdb ids in a text file::

  $ cat data/pdb_ids.txt
  1y26
  1fir

  $ while read p; do rna_pdb_tools.py --fetch_ba $p; done < data/pdb_ids.txt
  downloading...1y26_ba.pdb ok
  downloading...1fir_ba.pdb ok

  $ ls *.pdb
  1fir_ba.pdb 1y26_ba.pdb

.. autofunction:: rna_pdb_tools.pdb_parser_lib.fetch_ba

delete
-----------------------------------------

Examples::

    $ for i in `ls *pdb`; do rna_pdb_tools.py --delete A:46-56 $i > ../rpr_rm_loop/$i ; done

go over all files in the current directory, remove a fragment of chain A, residues between 46-56 (including them) and save outputs to in the folder `rpr_rm_loops`.

edit
-----------------------------------------

.. autofunction:: rna_pdb_tools.pdb_parser_lib.edit_pdb

the library
-----------------------------------------

  .. automodule:: rna_pdb_tools.pdb_parser_lib
   :members:

