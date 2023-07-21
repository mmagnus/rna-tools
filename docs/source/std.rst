Standardize your PDB files
===============================================

.. autoprogram:: rna_tools.rna_standardize:get_parser()
   :prog: rna_standardize.py

.. autoclass:: rna_tools.rna_tools_lib.RNAStructure
   :members: get_rnapuzzle_ready

Atoms order
-----------------------------------------------

Atoms order, A as an example::

	ATOM      1  P     G A   1      50.626  49.730  50.573  1.00100.19           P  
	ATOM      2  OP1   G A   1      49.854  48.893  49.562  1.00100.19           O  
	ATOM      3  OP2   G A   1      52.137  49.542  50.511  1.00 99.21           O  
	ATOM      4  O5'   G A   1      50.161  49.136  52.023  1.00 99.82           O  
	ATOM      5  C5'   G A   1      50.216  49.948  53.210  1.00 98.63           C  
	ATOM      6  C4'   G A   1      50.968  49.231  54.309  1.00 97.84           C  
	ATOM      7  O4'   G A   1      50.450  47.888  54.472  1.00 97.10           O  
	ATOM      8  C3'   G A   1      52.454  49.030  54.074  1.00 98.07           C  
	ATOM      9  O3'   G A   1      53.203  50.177  54.425  1.00 99.39           O  
	ATOM     10  C2'   G A   1      52.781  47.831  54.957  1.00 96.96           C  
	ATOM     11  O2'   G A   1      53.018  48.156  56.313  1.00 96.77           O  
	ATOM     12  C1'   G A   1      51.502  47.007  54.836  1.00 95.70           C  
	ATOM     13  N9    G A   1      51.628  45.992  53.798  1.00 93.67           N  
	ATOM     14  C8    G A   1      51.064  46.007  52.547  1.00 92.60           C  
	ATOM     15  N7    G A   1      51.379  44.966  51.831  1.00 91.19           N  
	ATOM     16  C5    G A   1      52.197  44.218  52.658  1.00 91.47           C  
	ATOM     17  C6    G A   1      52.848  42.992  52.425  1.00 90.68           C  
	ATOM     18  O6    G A   1      52.826  42.291  51.404  1.00 90.38           O  
	ATOM     19  N1    G A   1      53.588  42.588  53.534  1.00 90.71           N  
	ATOM     20  C2    G A   1      53.685  43.282  54.716  1.00 91.21           C  
	ATOM     21  N2    G A   1      54.452  42.733  55.671  1.00 91.23           N  
	ATOM     22  N3    G A   1      53.077  44.429  54.946  1.00 91.92           N  
	ATOM     23  C4    G A   1      52.356  44.836  53.879  1.00 92.62           C  

	ATOM     24  P     C A   2      54.635  50.420  53.741  1.00100.19           P  
	ATOM     25  OP1   C A   2      55.145  51.726  54.238  1.00100.19           O  


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

.. image :: ../pngs/op3.jpg # op3.pse 
