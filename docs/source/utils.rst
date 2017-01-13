Utils
=========================================

Apps
-----------------------------------------

ClaRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. automodule:: rna_pdb_tools.utils.clarna_app.clarna_app
   :members:

ROSETTA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run all modeling steps
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_run.get_parser
   :prog: rna_rosetta_run

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_run
   :members:

Minimize
``````````````````````````````````````````
.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_min.get_parser
   :prog: rna_rosetta_min

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_min
   :members:

Cluster
```````````````````````````````````````````

.. argparse::
   :ref: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_cluster.get_parser
   :prog: rna_rosetta_min

.. automodule:: rna_pdb_tools.utils.rna_rosetta.rna_rosetta_cluster
   :members:

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
     

Secondary Structure (secstruc)
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
