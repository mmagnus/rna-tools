rna-pdb-tools
========================================

Fetch
-----------------------------------------

Example::

  $ ./rna-pdb-tools.py --fetch 1xjr
  downloading...1xjr ok

.. autofunction:: rna_pdb_tools.pdb_parser_lib.fetch

Fetch Biological Assembly
-----------------------------------------

Example::

  $ rna-pdb-tools.py --fetch_ba 1xjr
  downloading...1xjr_ba.pdb ok

or over a list of pdb ids in a text file::

  $ cat data/pdb_ids.txt
  1y26
  1fir

  $ while read p; do rna-pdb-tools.py --fetch_ba $p; done < data/pdb_ids.txt
  downloading...1y26_ba.pdb ok
  downloading...1fir_ba.pdb ok

  $ ls *.pdb
  1fir_ba.pdb 1y26_ba.pdb

.. autofunction:: rna_pdb_tools.pdb_parser_lib.fetch_ba

Selection
-----------------------------------------

.. automodule:: rna_pdb_tools.utils.extra_functions.select_fragment
		:members:
Delete
-----------------------------------------

Examples::

	for i in `ls *pdb`; do rna-pdb-tools.py --delete A:46-56 $i > ../rpr_rm_loop/$i ; done

go over all files in the current directory, remove a fragment of chain A, residues between 46-56 (including them) and save outputs to in the folder `rpr_rm_loops`.

Edit
-----------------------------------------

.. autofunction:: rna_pdb_tools.pdb_parser_lib.edit_pdb

Examples::

  rna-pdb-tools.py --edit 'A:3-21>A:1-19' 1f27_clean.pdb > 1f27_clean_A1-19.pdb

or even::

  $ md_1f27_clx rna-pdb-tools.py --edit 'A:3-21>A:1-19,B:22-32>B:20-30' 1f27_clean.pdb > 1f27_clean_renumb.pdb

The library
-----------------------------------------

  .. automodule:: rna_pdb_tools.pdb_parser_lib
   :members:

