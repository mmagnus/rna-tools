Install & Configure
=============================================

Download rna-pdb-tools by clicking here https://github.com/mmagnus/rna-tools/archive/master.zip or using ``git``::

   $ git clone https://github.com/mmagnus/rna-tools.git
   $ cd rna-pdb-tools

``git`` is better if you want to contribute to the package and if you want to get pretty frequent updates.

The first step is "zero" because not all requirements are needed to start working with rna-tools.

To install the full set of requirements, use ``pip``:

0. ``pip install -r docs/requirements.txt``

and setup the package itself in three steps:

1. Add the path to the package to your PYTHONPATH (in ~/.bashrc or ~/.zshrc), e.g., ``export PYTHONPATH=$PYTHONPATH:/home/magnus/src/rna-tools/``
   
2. Add the path to the bin folder of the package to your PATH (in  ~/.bashrc or ~/.zshrc), e.g.,  ``export PATH=$PATH:/home/magnus/src/rna-tools/bin/``
   
3. Add the path to the bin folder of the package to your PATH (in  ~/.bashrc or ~/.zshrc), e.g.,  ``export RNA_TOOLS_PATH=/home/magnus/src/rna-tools/``

4. And run the installation script::

    ➜  rna-tools git:(master) ✗ ./install_links_bin.sh
    Installed in ./bin
    rmsd_calc_to_target.py

should be OK now :-)

To set you own configuration, please first::

    cp rna_tools_config_local.py_sample rna_tools_config_local.py # in rna-tools/rna_tools

and then edit ``rpt_config_local.py`` as you need. In my case it is::

    rna_tools git:(master) ✗ cat rna_tools_config_local.py
    VARNA_PATH  = '/Users/magnus/skills/rnax/varna_tut/'
    VARNA_JAR_NAME = 'VARNA.jar'


