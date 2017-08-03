Configuration
------------------------------------

Keep configuration sytanx like::

    from rna_pdb_tools.rpt_config import CPUS_CLUSTER
    # since we use export PYTHONPATH=$PYTHONPATH:/home/magnus/src/rna-pdb-tools/

vs::

    try:
        RNA_ROSETTA_RUN_ROOT_DIR_MODELING = os.environ['RNA_ROSETTA_RUN_ROOT_DIR_MODELING']
    except:
        print ('Set up RNA_ROSETTA_RUN_ROOT_DIR_MODELING in .bashrc')

Documentation
------------------------------------

We are using(at least we are moving towards) the Google style docstrings via Napoleon. Napoleon is a Sphinx Extensions that enables Sphinx to parse both NumPy and Google style docstrings - the style recommended by Khan Academy. http: // www.sphinx - doc.org / en / stable / ext / napoleon.html  # type-annotations

Add a new util
------------------------------------

1. Create a new folder in ``rna-pdb-tools/rna_pdb_tools/utils`` with your util. The folder will be seen online after your push at https://github.com/mmagnus/rna-pdb-tools/tree/master/rna_pdb_tools/utils. We will walk you through this simple example https://github.com/mmagnus/rna-pdb-tools/tree/master/rna_pdb_tools/utils/renum_pdb_to_aln .

2. Make sure that there is a simple test as ``test.sh``::

    #!/bin/bash
    python renum_pdb_to_aln.py - -residue_index_start 1 obj1 test_data / ALN_OBJ1_OBJ2.fa test_data / obj01.pdb

and there is a ``test_data`` folder with some test inputs and outputs. See the example.

3. Add your util to ``install_links_bin.sh`` at the top folder of ``rna-pdb-tools``::

    ln -s $curr_dir/rna_pdb_tools /utils/<util folder>/<util script name with .py> $curr_dir/bin/<util script name with .py>

    e.g.

    ln -s $curr_dir/rna_pdb_tools/utils/renum_pdb_to_aln/renum_pdb_to_aln.py $curr_dir/bin/rna_renum_pdb_to_aln.py

This will "install" you script in bin directory of the project so it can be used system-wide.

Run this script to see if there is any error, ``./install_links_bin.sh``.

4.  Add your util to the documentation. Go to ``rna-pdb-tools/docs/source`` and edit ``utils.rst``. Add - wherever you think your tool will fit - lines like::

     Renumber a pdb file according to alignment
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     .. argparse::
        :ref: rna_pdb_tools.utils.<util folder>.<utl script name>.get_parser
        :prog: <utl script name>

     .. automodule:: rna_pdb_tools.utils.<util folder>.<utl script name>
        : members:

     e.g:

      Renumber a pdb file according to alignment
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      .. argparse::
         :ref: rna_pdb_tools.utils.renum_pdb_to_aln.renum_pdb_to_aln.get_parser
         :prog: renum_pdb_to_aln

      .. automodule:: rna_pdb_tools.utils.renum_pdb_to_aln.renum_pdb_to_aln
        :members:

and run ``make html`` in the folder to check if the documenatation is compiled without any errors.

If you are using any external library such us ``scipy``, please make sure that they are listed in ``rna-pdb-tools/docs/requirements.txt``. If the library is not there, please add it. This file is read by the Read The Docs to compile the documentation and also by Travis for continuous testing.

You can open the documentation compiled locally under a link link file://<path to rna-pdb>/rna-pdb-tools/docs/build/html/index.html, e.g. file:///Users/magnus/work/src/rna-pdb-tools/docs/build/html/index.html.

5. The very last step is to add your util ``test.sh`` to the main testing script. Edit ``rna-pdb-tools/test.sh`` and add ::

       cd ./utils/<util folder>/
       ./test.sh
       cd ../..

      e.g.

      cd ./utils/renum_pdb_to_aln/
      ./test.sh
      cd ../..

6. Run this main test (``./test.sh``) and see if the util works as expected.

7. Now we are ready to push the changes. In the terminal, type::

     $ git pull
     $ git add <files> # or use git gui
     $ git commit -m <desc the util>
     $ git push

to commit all your changes and push it to the Github repository!

.. warning:: This testing is very, very rough and we are moving to have more test in py.test.
