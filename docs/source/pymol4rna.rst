PyMOL4RNA
=========================================

.. automodule:: rna_pdb_tools.utils.PyMOL4RNA.PyMOL4RNA
   :members:
   :undoc-members:

PyMOL Drawing
-----------------------------------------

.. automodule:: rna_pdb_tools.utils.pymol_drawing.pymol_drawing
   :members:
   :undoc-members:

Install PyMOL plugin to view the interactions with PyMOL::

    run <path>rna-pdb-tools/utils/pymol_drawing/pymol_dists.py

and type::

    draw_dists([[29, 41], [7, 66], [28, 42], [51, 63], [50, 64], [2, 71], [5, 68], [3, 70], [31, 39], [4, 69], [6, 67], [12, 23], [52, 62], [30, 40], [49, 65], [27, 43], [11, 24], [1, 72], [10, 25], [15, 48], [53, 61], [19, 56], [13, 22], [36, 37], [18, 19], [22, 46], [35, 73], [32, 38], [9, 13], [19, 20], [18, 20], [54, 60], [9, 23], [34, 35], [36, 38], [53, 54], [20, 56], [9, 12], [26, 44], [18, 55], [54, 61], [32, 36]])

.. image:: ../pngs/pymol_dists.png

RNA Helix Vis (draw helices using PyMOL)
-----------------------------------------

.. argparse::
   :ref: rna_pdb_tools.utils.rna_helix_vis.rna_helix_vis.get_parser
   :prog: rna_helix_vis

.. automodule:: rna_pdb_tools.utils.rna_helix_vis.rna_helix_vis
   :members:
   :undoc-members:

Install
-----------------------------------------

Open your `~/.pymolrc` and set up following variables as you need::

   # rna-pdb-tools
   RNA_PDB_TOOLS="/Users/magnus/work-src/rna-pdb-tools"
   EXECUTABLE="/bin/zsh" # set up your shell, usually /bin/bash or /bin/zsh
   SOURCE="source ~/.zshrc" # set up the path to the file where you keep your shell variables
   CLARNA_RUN="/Users/magnus/work-src/clarna_play/clarna_run.py"  # if you want to run clarna_run.py set up the path
   sys.path.append('/Users/magnus/work-src/rna-pdb-tools')
   run ~/work-src/rna-pdb-tools/rna_pdb_tools/utils/PyMOL4RNA/PyMOL4RNA.py
   run ~/work-src/rna-pdb-tools/rna_pdb_tools/utils/pymol_drawing/pymol_drawing.py
   run ~/work-src/rna-pdb-tools/rna_pdb_tools/utils/rna_filter/pymol_dists.py
