===========
 Workflows
===========

Example #1
-----------------------------------------------------------------------

The native::

	[mq] md_1msy_clx cat 1msy_clean.pdb.outCR
	Classifier: Clarna
	chains:  A 2647 2673
	A 2648   A 2672          bp G U                  WW_cis   0.8732
	A 2649   A 2671          bp C G                  WW_cis   0.9160
	A 2650   A 2670          bp U A                  WW_cis   0.9289
	A 2651   A 2669          bp C G                  WW_cis   0.9439
	A 2652   A 2668          bp C G                  WW_cis   0.9281
	A 2655   A 2656          bp G U                  SH_cis   0.9227
	A 2656   A 2665          bp U A                 WH_tran   0.8526
	A 2657   A 2664          bp A G                 HS_tran   0.8513
	A 2658   A 2663          bp C G                  WW_cis   0.9421
	A 2659   A 2662          bp G A                 SH_tran   0.7619

but analyzed structres are like::

	[mq] md_1msy_clx cat struc/1msy_rnakbmd_decoy1478_clx.pdb.outCR
	Classifier: Clarna
	chains: A 1 27
	 2       26          bp G U                  WW_cis   0.7196
	 3       25          bp C G                  WW_cis   0.6702
	 4       24          bp U A                  WW_cis   0.8911
	 5       23          bp C G                  WW_cis   0.8925
	 6       22          bp C G                  WW_cis   0.9026
	 9       10          bp G U                  SH_cis   0.8714
	 10       19          bp U A                 WH_tran   0.7279
	 11       18          bp A G                 HS_tran   0.8810
	 12       17          bp C G                  WW_cis   0.9115
	
You have to renumber 1msy_clean.pdb to 1:27::

    rna-pdb-tools.py --edit 'A:2647-2673>A:1:17' 1msy_clean.pdb > 1msy_clean_renumb.pdb

.. image:: ../pngs/edit.png

Example #2
-----------------------------------------------------------------------

Listing::
  
    $ rna-pdb-tools.py --get_seq 1nuj_rnakbmd_decoy1000_clx.pdb
    > 1nuj_rnakbmd_decoy1000_clx.pdb A:1-13
    CGGACCGAGCCAG
    > 1nuj_rnakbmd_decoy1000_clx.pdb B:14-24
    GCUGGGAGUCC

    $ rna-pdb-tools.py --get_seq 1nuj_clean.pdb
    > 1nuj_clean.pdb A:18-30
    CGGACCGAGCCAG
    > 1nuj_clean.pdb B:39-49
    GCUGGGAGUCC

    $ rna-pdb-tools.py --edit 'A:18-30>A:1-13,B:39-49>B:14-24' 1nuj_clean.pdb > 1nuj_clean_renumber.pdb

    $ rna-pdb-tools.py --get_seq 1nuj_clean_renumber.pdb
    > 1nuj_clean_renumber.pdb A:1-13
    CGGACCGAGCCAG
    > 1nuj_clean_renumber.pdb B:14-24
    GCUGGGAGUCC

Examples #3
------------------------------------------------------------------------

Starting structure don't have chain id:
 
.. code-block:: console
		  
  # add chain A
  $ parallel "rna_add_chain.py -c A {} > ../struc_with_chain/{}" ::: *.pdb
  # edit the second part of the new chain A as B
  $ parallel "rna-pdb-tools.py --edit 'A:14-27>B:14-27' {} > out/{}" ::: *.pdb

.. image:: ../pngs/1duq.png
