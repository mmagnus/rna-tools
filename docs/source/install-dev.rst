Install & Configure
=============================================

Pip as a developer
-------------------------------------------
You can get rna-tools and install them from the current directory with this pip::

    pip install -e git+http://github.com/mmagnus/rna-tools.git

This was is better (than ``pip install rna-tools``) if you're going to do some coding in the tools.

Configuration
------------------------------------------
To set up your own configuration, create ~/.rna_tools.py in your HOME directory and redefine variables, e.g.::

    (py37) [mx] rna-tools$ git:(master) ✗ cat ~/.rna_tools.py
    VARNA_PATH = '/Users/magnus/work/opt/varna'
    VARNA_JAR_NAME = 'VARNA.jar'
    SIMRNA_DATA_PATH = '/Users/magnus/work/opt/simRNA/SimRNA_64bitIntel_MacOSX_staticLibs/data'
    QRNAS_PATH = "/Users/magnus/work/opt/qrnas/"
    RCHIE_PATH = "/Users/magnus/work/opt/r-chie/"
    RFAM_DB_PATH = "/Users/magnus/work/db/rfam/Rfam.cm"
    CONTEXTFOLD_PATH = "/Users/magnus/work/opt/ContextFold_1_00/"
    DIFF_TOOL = "open -a diffmerge"
    CPUS_CLUSTER = 630
    RNASTRUCTURE_PATH = "/Users/magnus/work/opt/RNAstructure/6.1/"
    ENTRNA_PATH = "/Users/magnus/work/opt/ENTRNA"

All requirements
-------------------------------------------
To get ALL requirements, use ``pip``::

     pip install -r docs/requirements.txt

Be default rna-tools will not install all requirements, because some of them are heavy or might cause various problems, so you will be asked to install them when needed. For example, installing `matplotlib` is not essential for many other tools, `python-Levenshtein` is only used in one function.

Test all
-------------------------------------------
This is still under active development.

To test (almost) all rna-tools functionality, you can run ``rna_tools_test_all.py`` to see if you got any errors, this should look like::

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

For crude testing you can also use ``./test.sh`` script and then see for errors in output and also check output/ folder to see if there are differences between your output and output committed to GitHub by me.
