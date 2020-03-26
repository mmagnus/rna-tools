Install & Configure
=============================================

Download rna-tools by clicking here https://github.com/mmagnus/rna-tools/archive/master.zip (~30 MB) and unpack the zip file, enter the folder::

   $ cd rna-tools-master
   
OR use ``git``::

   $ git clone https://github.com/mmagnus/rna-tools.git
   $ cd rna-tools

``git`` is better if you want to contribute to the package or/and you want to get pretty frequent updates.

The first step is "zero" because not all requirements are needed to start working with rna-tools. If anything is missing you can install it later.

To install the full set of requirements, use ``pip``:

0. ``pip install -r docs/requirements.txt``

1. Setup the package paths by adding to your  ~/.bashrc or ~/.zshrc following code::

      export RNA_TOOLS_PATH=<PATH TO YOUR RNA_TOOLS>
      export PYTHONPATH=$PYTHONPATH:$RNA_TOOLS_PATH
      export PATH=$PATH:$RNA_TOOLS_PATH'/bin/'

for example in my case it looks as::

   export RNA_TOOLS_PATH=/home/magnus/work-src/rna-tools/
   export PYTHONPATH=$PYTHONPATH:$RNA_TOOLS_PATH
   export PATH=$PATH:$RNA_TOOLS_PATH'/bin/'

2. And run the installation script::

    ➜  rna-tools git:(master) ✗ ./install_links_bin.sh
    Installed in ./bin
    rmsd_calc_to_target.py

3. Run ``rna_tools_test_all.py`` to see if you got any errors, this should look like::

      (py37) [mm] rna-tools$ git:(master) rna_tools_test_all.py
      BlastPDB requires urllib3
      - Python: 3.7.4 (default, Aug 13 2019, 15:17:50) [Clang 4.0.1 (tags/RELEASE_401/final)]
      - rna-tools: b'py2-78-g3b3dd5f'
      - RNA_TOOLS_PATH set to  /home/magnus/work-src/rna-tools/
      - See full list of tools <https://github.com/mmagnus/rna-tools/blob/master/rna-tools-index.csv
      Seems OK

or for Python 2::

   (base) [mm] rna-tools$ git:(master) rna_tools_test_all.py
   - Python: 2.7.16 |Anaconda, Inc.| (default, Mar 14 2019, 16:24:02) [GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)]
   - rna-tools: py2-78-g3b3dd5f
   - RNA_TOOLS_PATH set to  /home/magnus/work-src/rna-tools/
   - See full list of tools <https://github.com/mmagnus/rna-tools/blob/master/rna-tools-index.csv
   Seems OK

To set you own configuration, please first::

    cp rna_tools_config_local.py_sample rna_tools_config_local.py # in rna-tools/rna_tools

and then edit ``rna_tools_config_local.py`` as you need. In my case it is::

    rna_tools git:(master) ✗ cat rna_tools_config_local.py
    VARNA_PATH  = '/Users/magnus/skills/rnax/varna_tut/'
    VARNA_JAR_NAME = 'VARNA.jar'

For git based pip run::

   pip install git+http://github.com/mmagnus/rna-tools.git
