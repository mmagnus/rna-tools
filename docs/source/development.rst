Configuration
------------------------------------

Keep configuration sytanx like::

    from rna_tools.rna_tools_config import CPUS_CLUSTER
    # since we use export PYTHONPATH=$PYTHONPATH:/home/magnus/src/rna-tools/

vs::

    try:
        RNA_ROSETTA_RUN_ROOT_DIR_MODELING = os.environ['RNA_ROSETTA_RUN_ROOT_DIR_MODELING']
    except:
        print ('Set up RNA_ROSETTA_RUN_ROOT_DIR_MODELING in .bashrc')

Documentation
------------------------------------

We are using (at least we are moving towards)the Google style docstrings via Napoleon. Napoleon is a Sphinx Extensions that enables Sphinx to parse both NumPy and Google style docstrings - the style recommended by Khan Academy. http://www.sphinx-doc.org/en/stable/ext/napoleon.html#type-annotations

Add a new tool to the package
------------------------------------

1. Create a new folder in ``rna-tools/rna_tools/tools`` with your tool. The folder will be seen online after your push at https://github.com/mmagnus/rna-tools/tree/master/rna_tools/tools. We will walk you through this simple example https://github.com/mmagnus/rna-tools/tree/master/rna_tools/tools/renum_pdb_to_aln .

2. Make sure that there is a simple test as ``test.sh``::

    #!/bin/bash
    python renum_pdb_to_aln.py --residue_index_start 1 obj1 test_data/ALN_OBJ1_OBJ2.fa test_data/obj01.pdb

and there is a ``test_data`` folder with some test inputs and outputs. See the example.

3. Add your tool to ``install_links_bin.sh`` at the top folder of ``rna-tools``::

    ln -s $curr_dir/rna_tools /tools/<tool folder>/<util script name with .py> $curr_dir/bin/<util script name with .py>

    e.g.

    ln -s $curr_dir/rnatools/utils/renum_pdb_to_aln/renum_pdb_to_aln.py $curr_dir/bin/rna_renum_pdb_to_aln.py

This will "install" you script in bin directory of the project so it can be used system-wide.

Run this script to see if there is any error, ``./install_links_bin.sh``.

4.  Add your tool to the documentation. The tool has to be "importable", so don't forget to create ``__init__.py`` inside your tool directory. Next, go to ``rna-tools/docs/source`` and edit ``tools.rst``. Add, wherever you think your tool will fit, lines like::

     Renumber a pdb file according to alignment
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     .. argparse::
        :ref: rna_tools.tools.<tool folder>.<tool script name>.get_parser
        :prog: <utl script name>

     .. automodule:: rna_tools.tools.<tool folder>.<tool script name>
        : members:

     e.g.:

      Renumber a pdb file according to alignment
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      .. argparse::
         :ref: rna_tools.tools.renum_pdb_to_aln.renum_pdb_to_aln.get_parser
         :prog: renum_pdb_to_aln

      .. automodule:: rna_tools.tools.renum_pdb_to_aln.renum_pdb_to_aln
        :members:

and run ``make html`` in the folder to check if the documentation is compiled without any errors.

If you are using any external library such us ``scipy``, please make sure that they are listed in ``rna-tools/docs/requirements.txt``. If the library is not there, please add it. This file is read by the Read The Docs to compile the documentation online and also by Travis for continuous testing.

You can open the documentation compiled locally under a link link file://<path to rna-tools>/rna-tools/docs/build/html/index.html, e.g. file:///Users/magnus/work/src/rna-tools/docs/build/html/index.html.

5. The very last step is to add your tool ``test.sh`` to the main testing script. Edit ``rna-tools/test.sh`` and add ::

       cd ./tools/<tool folder>/
       ./test.sh
       cd ../..

      e.g.

      cd ./tools/renum_pdb_to_aln/
      ./test.sh
      cd ../..

6. Run this main test (``./test.sh``) and see if the tool works as expected.

7. Now we are ready to push the changes. In the terminal, type::

     $ git pull
     $ git add <files> # or use git gui
     $ git commit -m <desc the tool>
     $ git push

to commit all your changes and push it to the Github repository!

.. warning:: This testing is very, very rough and we are moving to have more test in py.test at some point.
