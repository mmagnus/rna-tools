Standardize your PDB files
===============================================

.. autoprogram:: rna_tools.rna_standardize:get_parser()
   :prog: rna_standardize.py

.. autoclass:: rna_tools.rna_tools_lib.RNAStructure
   :members: get_rnapuzzle_ready

OP3
------------------------------------------------

The first residue will have only OP1 and OP2 (OP3 will be removed)::

	ATOM      1  OP3   G A   1      50.193  51.190  50.534  1.00 99.85           O  
	ATOM      2  P     G A   1      50.626  49.730  50.573  1.00100.19           P  
	ATOM      3  OP1   G A   1      49.854  48.893  49.562  1.00100.19           O  
	ATOM      4  OP2   G A   1      52.137  49.542  50.511  1.00 99.21           O  
	ATOM      5  O5'   G A   1      50.161  49.136  52.023  1.00 99.82           O  
	ATOM      6  C5'   G A   1      50.216  49.948  53.210  1.00 98.63           C  
	ATOM      7  C4'   G A   1      50.968  49.231  54.309  1.00 97.84           C  

Listing. An example: `1ehz.pdb`.
