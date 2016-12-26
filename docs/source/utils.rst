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

.. argparse::
   :ref: rna_pdb_tools.utils.rna_run_rosetta.rna_run_rosetta.get_parser
   :prog: rna_run_rosetta

.. automodule:: rna_pdb_tools.utils.rna_run_rosetta.rna_run_rosetta
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

  magnus@peyote2:~/src/cluster_load$ ./load
  jobs: 743 load (1k max):  0.743  to use: 257
  jobs: 200 load (1k max):  0.2  -- magnus
  1 magnus       r 200        
  4 gchojnowski  r 10        
  17 wdawson      r 10        
  4 bharat       r 10        
  29 mboni        r 10        
  2 mmaciejczyk  r 1        

<https://github.com/mmagnus/rna-pdb-tools/tree/master/rna_pdb_tools/utils/cluster_load>

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
