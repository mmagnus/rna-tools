from distutils.core import setup

setup(
    name='rna-pdb-tools',
    version='',
    packages=['', 'rna_tools., 'rna_tools.tools', 'rna_tools.tools.misc', 'rna_tools.tools.rna_bp',
              'rna_tools.tools.diffpdb', 'rna_tools.tools.clanstix', 'rna_tools.tools.plotting',
              'rna_tools.tools.ClashCalc', 'rna_tools.tools.PyMOL4RNA', 'rna_tools.tools.rna_x3dna',
              'rna_tools.tools.clarna_app', 'rna_tools.tools.rna_filter', 'rna_tools.tools.rmsd_signif',
              'rna_tools.tools.rna_rosetta', 'rna_tools.tools.rnakb_tools', 'rna_tools.tools.cluster_load',
              'rna_tools.tools.pdb_formatix', 'rna_tools.tools.pdb_formatix.test',
              'rna_tools.tools.rna_calc_inf', 'rna_tools.tools.pymol_drawing',
              'rna_tools.tools.rna_alignment', 'rna_tools.tools.rna_calc_rmsd',
              'rna_tools.tools.rna_calc_rmsd.lib', 'rna_tools.tools.rna_calc_rmsd.lib.rmsd',
              'rna_tools.tools.rna_helix_vis', 'rna_tools.tools.rna_refinement',
              'rna_tools.tools.extra_functions', 'rna_tools.tools.renum_pdb_to_aln',
              'rna_tools.tools.rna_calc_evo_rmsd', 'rna_tools.tools.simrna_trajectory',
              'rna_tools.tools.rna_calc_rmsd_trafl', 'rna_tools.tools.rna_sali2dotbracket',
              'rna_tools.tools.pdbs_measure_atom_dists', 'rna_tools.tools.rna_convert_pseudoknot_formats',
              'rna_tools.tools.rna_pdb_edit_occupancy_bfactor',
              'rna_tools.tools.rna_pdb_merge_structure_with_fragments'],
    url='https://github.com/mmagnus/rna-pdb-tools',
    license='GPLv3',
    author='Marcin Magnus',
    author_email='',
    description='rna-pdb-tools: a toolbox to analyze structures and simulations of RNA'
)
