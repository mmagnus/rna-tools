RNA Tools
=========================================

RNA Sequence
------------------------------------------
.. automodule:: rna_tools.Seq
   :members:
   :undoc-members:

RNA Secondary Structure
------------------------------------------
.. automodule:: rna_tools.SecondaryStructure
   :members:
   :undoc-members:

rna_dot2ct.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. argparse::
   :ref: rna_tools.rna_dot2ct.get_parser
   :prog: rna_dot2ct.py

Secondary structure format conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. image:: ../../rna_tools/tools/rna_convert_pseudoknot_formats/doc/varna_2pk.png

.. automodule:: rna_tools.tools.rna_convert_pseudoknot_formats.rna_ss_pk_to_simrna
   :members:
   :undoc-members:

Search
------------------------------------------
Blast PDB
~~~~~~~~~~~~~~~~~~~~

A super-simple wrapper around Blast on the PDB db (online).

.. automodule:: rna_tools.BlastPDB
   :members:
   :undoc-members:

Rfam Search
~~~~~~~~~~~~~~~~~~~~

A super-simple wrapper around cmscan (Infernal) on local Rfam database.

.. automodule:: rna_tools.RfamSearch
   :members:
   :undoc-members:

RNA Alignment
------------------------------------------

.. automodule:: rna_tools.tools.rna_alignment.rna_alignment

RNASeq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_tools.tools.rna_alignment.rna_alignment.RNASeq
   :members:
   :undoc-members:

RNAalignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_tools.tools.rna_alignment.rna_alignment.RNAalignment
   :members:
   :undoc-members:


rna_alignment_get_species.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.utils.rna_alignment_get_species.get_parser
   :prog: rna_alignment_get_species

rna_alignment_calc_energy.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.utils.rna_alignment_calc_energy.get_parser
   :prog: rna_alignment_calc_energy

rna_align_get_ss_from_fasta.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.rna_align_get_ss_from_fasta.get_parser
   :prog: rna_align_get_ss_from_fasta

rna_align_get_ss_from_stk.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.rna_align_get_ss_from_stk.get_parser
   :prog: rna_align_get_ss_from_stk

rna_align_distance_to_seq.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.rna_align_distance_to_seq.get_parser
   :prog: rna_align_distance_to_seq.py

rna_align_foldability.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.rna_align_foldability.get_parser
   :prog: rna_align_foldability.py


Random assignment of nucleotides
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_alignment.random_assignment_of_nucleotides.get_parser
   :prog: random_assignment_of_nucleotides

.. automodule:: rna_tools.tools.rna_alignment.random_assignment_of_nucleotides
   :members:
   :undoc-members:


CMAlign
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_tools.tools.rna_alignment.rna_alignment.CMAlign
   :members:
   :undoc-members:


RChie
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: rna_tools.tools.rna_alignment.rna_alignment.RChie
   :members:
   :undoc-members:

Renumber a pdb file according to alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.renum_pdb_to_aln.renum_pdb_to_aln.get_parser
   :prog: renum_to_aln

.. automodule:: rna_tools.tools.renum_pdb_to_aln.renum_pdb_to_aln
   :members:
   :undoc-members:

Root Mean Square Deviation (RMSD)
---------------------------------------------

rna_calc_rmsd
~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_calc_rmsd.rna_calc_rmsd.get_parser
   :prog: rna_calc_rmsd

.. automodule:: rna_tools.tools.rna_calc_rmsd.rna_calc_rmsd
   :members:
   :undoc-members:

rna_calc_rmsd_multi_targets
~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_calc_rmsd.rna_calc_rmsd_multi_targets.get_parser
   :prog: rna_calc_evo_rmsd

rna_calc_rmsd_trafl
~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_calc_rmsd_trafl.rna_calc_rmsd_trafl.get_parser
   :prog: rna_calc_evo_rmsd

.. argparse::
   :ref: rna_tools.tools.rna_calc_rmsd_trafl.rna_cal_rmsd_trafl_plot.get_parser
   :prog: rna_cal_rmsd_trafl_plot


rna_calc_rmsd_all_vs_all
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_calc_rmsd.rna_calc_rmsd_all_vs_all.get_parser
   :prog: rna_calc_rmsd_all_vs_all

.. automodule:: rna_tools.tools.rna_calc_rmsd.rna_calc_rmsd_all_vs_all
   :members:
   :undoc-members:

Interaction Network Fidelity (INF)
------------------------------------------------------------

rna_calc_inf
~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_calc_inf.rna_calc_inf.get_parser
   :prog: rna_calc_inf

.. automodule:: rna_tools.tools.rna_calc_inf.rna_calc_inf
   :members:
   :undoc-members:

rna_calc_dinf
~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_calc_inf.rna_calc_dinf.get_parser
   :prog: rna_calc_inf

.. automodule:: rna_tools.tools.rna_calc_dinf.rna_calc_dinf
   :members:
   :undoc-members:

RNA filter (DCA)
-----------------------------------------

rna_filter.py - calculate distances based on given restrants on PDB files or SimRNA trajectories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_filter.rna_filter.get_parser
   :prog: rna_filter.py

rna_dca_mapping.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_filter.rna_dca_mapping.get_parser
   :prog: rna_dca_mapping.py


show_dists - show distances in PyMOL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_tools.tools.rna_filter.pymol_dists
   :undoc-members:

rna_ex2x.py - analyze an evolutionary coupling file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

.. argparse::
   :ref: rna_tools.tools.rna_filter.rna_ec2x.get_parser
   :prog: rna_ec2x.py

rna_pairs2SimRNArestrs.py - convert pairs to SimRNA restraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_filter.rna_pairs2SimRNArestrs.get_parser
   :prog: rna_pairs2SimRNArestrs.py

rna_ss_get_bps.py - get a list of base pairs for a given "fasta ss" file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_filter.rna_ss_get_bps.get_parser
   :prog: rna_ss_get_bps

rna_pairs_diff.py - get a diff of pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_filter.rna_pairs_diff.get_parser
   :prog: rna_pairs_diff.py

Contacts classification & secondary structure detection
------------------------------------------------------------------

3DNA (contacts classification & secondary structure detection)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_tools.tools.rna_x3dna.rna_x3dna
   :members:
   :undoc-members:

ClaRNA (contacts classification)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to calculate "Interaction Network Fidelity (INF) and not only" see rna_calc_inf in the Utils.

.. automodule:: rna_tools.tools.clarna_app.clarna_app
   :members:
   :undoc-members:


RNA 3D structure prediction
-----------------------------------------

ROSETTA
~~~~~~~~

A set of wrappers around Rosetta (https://www.rosettacommons.org/), mostly based on C. Y. Cheng, F. C. Chou, and R. Das, Modeling complex RNA tertiary folds with Rosetta, 1st ed., vol. 553. Elsevier Inc., 2015. http://www.sciencedirect.com/science/article/pii/S0076687914000524

Run (modeling)
^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_run.get_parser
   :prog: rna_rosetta_run.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_run.py
   :members:
   :undoc-members:

Get a number of structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_n.get_parser
   :prog: rna_rosetta_n.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_n.py
   :members:
   :undoc-members:

Get a head of a Rosetta silent file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_head.get_parser
   :prog: rna_rosetta_n.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_head.py
   :members:
   :undoc-members:

Cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :nodescription:
   :nosubcommands:
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_cluster.get_parser
   :prog: rna_rosetta_cluster.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_cluster
   :members:
   :undoc-members:

Minimize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_min.get_parser
   :prog: rna_rosetta_min.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_min.py
   :members:
   :undoc-members:

Extract lowscore decoy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_extract_lowscore_decoys.get_parser
   :prog: rna_rosetta_extract_lowscore_decoys.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_extract_lowscore_decoys.py
   :members:
   :undoc-members:

Check progress
^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.rna_rosetta.rna_rosetta_check_progress.get_parser
   :prog: rna_rosetta_cluster.py

.. automodule:: rna_tools.tools.rna_rosetta.rna_rosetta_check_progress.py
   :members:
   :undoc-members:

SimRNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Select low energy frames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.simrna_trajectory.rna_simrna_lowest.get_parser
   :prog: rna_simrna_lowest.py

Extract
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.tools.simrna_trajectory.rna_simrna_extract.get_parser
   :prog: rna_simrna_extract.py

SimRNAweb
~~~~~~~~~~~~~

Download files of a SimRNAweb run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. argparse::
   :ref: rna_tools.rna_simrnaweb_download_job.get_parser
   :prog: rna_simrnaweb_download_job.py

SimRNATrajectory
~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_tools.tools.simrna_trajectory.simrna_trajectory
   :members:
   :undoc-members:

RNA Refinement (QRNAS)
-----------------------------------------

.. argparse::
   :ref: rna_tools.tools.rna_refinement.rna_refinement.get_parser
   :prog: rna_refinement.py


diffpdb
-----------------------------------------

.. automodule:: rna_tools.tools.diffpdb.diffpdb
   :undoc-members:

RNA clustering with CLANS (clanstix)
------------------------------------------

.. automodule:: rna_tools.tools.clanstix.rna_clanstix
   :members:
   :undoc-members:

Misc
-----------

Plotting
~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.plotting.rna_plot_hist.get_parser
   :prog: rna_plot_hist

.. argparse::
   :ref: rna_tools.tools.plotting.rna_plot_density.get_parser
   :prog: rna_plot_density


rna_sali2dotbracket
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. argparse::
   :ref: rna_tools.tools.rna_sali2dotbracket.rna_sali2dotbracket.get_parser
   :prog: rna_sali2dotbracket

.. automodule:: rna_tools.tools.rna_sali2dotbracket.rna_sali2dotbracket
   :members:
   :undoc-members:

Cluster load
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. automodule:: rna_tools.tools.cluster_load.cluster_load
   :members:
   :undoc-members:

RNAkb
~~~~~~~~~~~~~~~~~~~~

.. automodule:: rna_tools.tools.rnakb_utils.rnakb_utils
   :members:
   :undoc-members:
