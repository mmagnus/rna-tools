Configuration
------------------------------------

Keep configuration sytanx like ::

     from rna_pdb_tools.rpt_config import CPUS_CLUSTER
     # since we use export PYTHONPATH=$PYTHONPATH:/home/magnus/src/rna-pdb-tools/

vs ::

  try:
    RNA_ROSETTA_RUN_ROOT_DIR_MODELING = os.environ['RNA_ROSETTA_RUN_ROOT_DIR_MODELING']
  except:
    print ('Set up RNA_ROSETTA_RUN_ROOT_DIR_MODELING in .bashrc')

Documentation
------------------------------------

We are using (at least we are moving towards) the Google style docstrings via Napoleon. Napoleon is a Sphinx Extensions that enables Sphinx to parse both NumPy and Google style docstrings - the style recommended by Khan Academy. http://www.sphinx-doc.org/en/stable/ext/napoleon.html#type-annotations
