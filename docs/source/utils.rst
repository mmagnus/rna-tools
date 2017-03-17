Utils
=========================================

Apps
-----------------------------------------

SimRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select low energy frames
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.simrna_trajectory.rna_simrna_lowest.get_parser
   :prog: rna_simrna_lowest.py
	  
Extract
``````````````````````````````````````````

.. argparse::
   :ref: rna_pdb_tools.utils.simrna_trajectory.rna_simrna_extract.get_parser
   :prog: rna_simrna_extract.py
	  
ClaRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. automodule:: rna_pdb_tools.utils.clarna_app.clarna_app
   :members:

ROSETTA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run (modeling)
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_run.get_parser
   :prog: rna_rosetta_run.py

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_run
   :members:

Get a number of structures
```````````````````````````````````````````

.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_n.get_parser
   :prog: rna_rosetta_n.py

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_n
   :members:

Minimize
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_min.get_parser
   :prog: rna_rosetta_min.py

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_min
   :members:

Cluster
```````````````````````````````````````````

.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_cluster.get_parser
   :prog: rna_rosetta_cluster.py

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_cluster
   :members:

RNA Alignment
------------------------------------------

RNASeq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_pdb_tools.utils.rna_alignment.rna_alignment.RNASeq
   :members:

RNAalignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_pdb_tools.utils.rna_alignment.rna_alignment.RNAalignment
   :members:


CMAlign
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_pdb_tools.utils.rna_alignment.rna_alignment.CMAlign
   :members:

RChie
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_pdb_tools.utils.rna_alignment.rna_alignment.RChie
   :members:

Shell utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rna_align_get_ss_from_alignment.py
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.rna_alignment.rna_align_get_ss_from_alignment.get_parser
   :prog: rna_align_get_ss_from_alignment.py

rna_align_seq_to_alignment.py
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.rna_alignment.rna_align_seq_to_alignment.get_parser
   :prog: rna_align_seq_to_alignment.py

rna_align_find_seq_in_alignment.py
```````````````````````````````````````````          
.. argparse::
   :ref: rna_pdb_tools.utils.rna_alignment.rna_align_find_seq_in_alignment.get_parser
   :prog: rna_align_find_seq_in_alignment.py

rna_align_find_core.py
`````````````````````````````````````````

.. argparse::
   :ref: rna_pdb_tools.utils.rna_alignment.rna_align_find_core.get_parser
   :prog: rna_align_find_core.py



RNA Helix Vis
-----------------------------------------

.. argparse::
   :ref: rna_pdb_tools.utils.rna_helix_vis.rna_helix_vis.get_parser
   :prog: rna_helix_vis

.. automodule:: rna_pdb_tools.utils.rna_helix_vis.rna_helix_vis
   :members:

Calculate RMSD
-----------------------------------------

rna_calc_rmsd
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_pdb_tools.utils.rna_calc_rmsd.rna_calc_rmsd
   :members:

rna_calc_rmsd_all_vs_all
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_pdb_tools.utils.rna_calc_rmsd.rna_calc_rmsd_all_vs_all
   :members:

Calculate Interaction Network Fidelity (INF) and not only
----------------------------------------------------------

rna_calc_inf
~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_pdb_tools.utils.rna_calc_inf.rna_calc_inf.get_parser
   :prog: rna_calc_inf

.. automodule:: rna_pdb_tools.utils.rna_calc_inf.rna_calc_inf
   :members:

diffpdb
-----------------------------------------

.. image:: ../../rna_pdb_tools/utils/diffpdb/doc/diffpdb_osx_diffmerge.png

rna_filter
-----------------------------------------

.. automodule:: rna_pdb_tools.utils.rna_filter.rna_filter
   :members:

SimRNATrajectory
-----------------------------------------

.. automodule:: rna_pdb_tools.utils.simrna_trajectory.simrna_trajectory
   :members:

PyMOL drawing
-----------------------------------------

.. automodule:: rna_pdb_tools.utils.pymol_drawing.pymol_drawing
   :members:
     


Secondary structure format conversion
-----------------------------------------

rna\_convert\_pseudoknot\_formats

Run this as::

	python rna-pk-simrna-to-one-line.py test_data/simrna.ss

Convert::

	> a
	....((.(..(((....)))..((((.(.........).)))....).)).......((((......))))
	..............................((((...................))))..............

to::

	> a
	....((.(..(((....)))..((((.(..[[[[...).)))....).))...]]]]((((......))))

and::

	>2 chains
	(((((......)))))........(.((....(.......)..(((. .)))...)).)
	.....((((((......................))))))........ ...........

to::

	>2 chains
	((((([[[[[[)))))........(.((....(]]]]]].)..(((. .)))...)).)

and::

	> b
	..(.......(((....)))..((((.(.........).))))).............((((......))))
	....((.(......................................).)).....................
	..............................((((...................))))..............

to::

	> b
	..(.[[.[..(((....)))..((((.(..{{{{...).)))))..].]]...}}}}((((......))))
	
and it works with VARNA:

.. image:: ../../rna_pdb_tools/utils/rna_convert_pseudoknot_formats/doc/varna_2pk.png

.. automodule:: rna_pdb_tools.utils.rna_convert_pseudoknot_formats.rna_ss_pk_to_simrna
   :members:


Secondary structure (secstruc)
-----------------------------------------

.. automodule:: rna_pdb_tools.utils.secstruc.secstruc
   :members:

Cluster load
-----------------------------------------

A very simple tool to see your cluster load per user::

  MAX_JOBS: 1000
  #jobs cluster 917 load:  0.917  to use: 83
  #jobs you     749 load:  0.749  to use: 251
  {'deepak': 160, 'azyla': 8, 'magnus': 749}
  1 azyla        r 8
  20 magnus       r 10
  16 deepak       r 10
  329 magnus       r 1
  22 magnus       qw 10

.. automodule:: rna_pdb_tools.utils.cluster_load.cluster_load
   :members:

Misc
------------------------------------------

rna_sali2dotbracket
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_pdb_tools.utils.rna_sali2dotbracket.rna_sali2dotbracket.get_parser
   :prog: rna_sali2dotbracket
   
.. automodule:: rna_pdb_tools.utils.rna_sali2dotbracket.rna_sali2dotbracket
   :members:

rna_add_chain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_pdb_tools.utils.misc.rna_add_chain
   :members:

.. argparse::
   :ref: rna_pdb_tools.utils.misc.rna_add_chain.get_parser
   :prog: rna_add_chain
