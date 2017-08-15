Install & Configure
=============================================

Download rna-pdb-tools by clicking here https://github.com/mmagnus/rna-pdb-tools/archive/master.zip or using ``git``::

   $ git clone https://github.com/mmagnus/rna-pdb-tools.git
   $ cd rna-pdb-tools

(``git`` is better if you want to contribute to the package and if you want to get pretty frequent updates).

The first step is "zero" because not all requirments are needed to start working with rna-pdb-tools.

To install the full set of requirements, use ``pip``:

0. ``pip3 install -r docs/requirements.txt``

and install the package itself in three steps:

1. add the path to the package to your PYTHONPATH (in ~/.bashrc), e.g. ``export PYTHONPATH=$PYTHONPATH:/home/magnus/src/rna-pdb-tools/``
   
2. add the path to the bin folder of the package to your PATH (in ~/.bashrc), e.g.  ``export PATH=$PATH:/home/magnus/src/rna-pdb-tools/bin/``
   
3. add the path to the bin folder of the package to your PATH (in ~/.bashrc), e.g.  ``export RNA_PDB_TOOLS=/home/magnus/src/rna-pdb-tools/``

4. and run the install script::

    ➜  rna-pdb-tools git:(master) ✗ ./install_links_bin.sh
    Installed in ./bin
    rmsd_calc_to_target.py

should be OK now :-)

To set you own configuration, please first::

    cp rpt_config_local.py_sample rpt_config_local.py # in rna-pdb-tools/rna_pdb_tools

and then edit ``rpt_config_local.py`` as you need. In my case it is::

    rna_pdb_tools git:(master) ✗ cat rpt_config_local.py
    VARNA_PATH  = '/Users/magnus/skills/rnax/varna_tut/'
    VARNA_JAR_NAME = 'VARNA.jar'


