from distutils.core import setup

setup(
    name='rna-pdb-tools',
    version='',
    packages=['', 'rna_pdb_tools', 'rna_pdb_tools.utils', 'rna_pdb_tools.utils.misc', 'rna_pdb_tools.utils.rna_bp',
              'rna_pdb_tools.utils.diffpdb', 'rna_pdb_tools.utils.clanstix', 'rna_pdb_tools.utils.plotting',
              'rna_pdb_tools.utils.ClashCalc', 'rna_pdb_tools.utils.PyMOL4RNA', 'rna_pdb_tools.utils.rna_x3dna',
              'rna_pdb_tools.utils.clarna_app', 'rna_pdb_tools.utils.rna_filter', 'rna_pdb_tools.utils.rmsd_signif',
              'rna_pdb_tools.utils.rna_rosetta', 'rna_pdb_tools.utils.rnakb_utils', 'rna_pdb_tools.utils.cluster_load',
              'rna_pdb_tools.utils.pdb_formatix', 'rna_pdb_tools.utils.pdb_formatix.test',
              'rna_pdb_tools.utils.rna_calc_inf', 'rna_pdb_tools.utils.pymol_drawing',
              'rna_pdb_tools.utils.rna_alignment', 'rna_pdb_tools.utils.rna_calc_rmsd',
              'rna_pdb_tools.utils.rna_calc_rmsd.lib', 'rna_pdb_tools.utils.rna_calc_rmsd.lib.rmsd',
              'rna_pdb_tools.utils.rna_helix_vis', 'rna_pdb_tools.utils.rna_refinement',
              'rna_pdb_tools.utils.extra_functions', 'rna_pdb_tools.utils.renum_pdb_to_aln',
              'rna_pdb_tools.utils.rna_calc_evo_rmsd', 'rna_pdb_tools.utils.simrna_trajectory',
              'rna_pdb_tools.utils.rna_calc_rmsd_trafl', 'rna_pdb_tools.utils.rna_sali2dotbracket',
              'rna_pdb_tools.utils.pdbs_measure_atom_dists', 'rna_pdb_tools.utils.rna_convert_pseudoknot_formats',
              'rna_pdb_tools.utils.rna_pdb_edit_occupancy_bfactor',
              'rna_pdb_tools.utils.rna_pdb_merge_structure_with_fragments'],
    url='https://github.com/mmagnus/rna-pdb-tools',
    license='GPLv3',
    author='Marcin Magnus',
    author_email='',
    description='rna-pdb-tools: a toolbox to analyze structures and simulations of RNA'
)
