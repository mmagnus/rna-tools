from setuptools import setup, find_packages
import versioneer

# read the contents of your README file
with open('README.md') as f:
    long_description = f.read()
    
setup(
    name='rna_tools',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    test_suite="tests",
    url='https://github.com/mmagnus/rna-tools',
    scripts=['rna_tools/rna_pdb_tools.py',
             
             'rna_tools/rna_calc_fenergy.py',
             'rna_tools/rna_pdb_replace.py',
             'rna_tools/rna_standardize.py',
             'rna_tools/rna_pdb_inspect.py',
             'rna_tools/rna_cif2pdb.py',
             'rna_tools/rna_pdb_fetch_header.py',             

             'rna_tools/tools/misc/rna_tools_which.py',
             'rna_tools/tools/misc/rna_tools_demo.py',

             'rna_tools/tools/workflow/rna_workflow.py',

             'rna_tools/tools/rna_calc_rmsd/rna_calc_rmsd.py',
             'rna_tools/tools/rna_calc_rmsd/rna_calc_rmsd_all_vs_all.py',
             'rna_tools/tools/rna_calc_rmsd/rna_calc_rmsd_biopython.py',
             
             'rna_tools/rna_torsions.py',

             'rna_tools/tools/diffpdb/diffpdb.py',
             'rna_tools/tools/clanstix/rna_clanstix.py',
             'rna_tools/rna_dot2ct.py',
             'rna_tools/rna_secondary_structure_prediction.py',
             'rna_tools/tools/rna_prediction_significance/rna_prediction_significance.py',
             'rna_tools/tools/rna_multimodels/rna_pdb_merge_into_one.py',
             'rna_tools/tools/rna_calc_inf/rna_calc_inf.py',
             'rna_tools/tools/rna_convert_pseudoknot_formats/rna_pk_simrna_to_one_line.py',
             'rna_tools/tools/clarna_app/rna_clarna_app.py',
             'rna_tools/tools/rna_helix_vis/rna_helix_vis.py',
             'rna_tools/tools/misc/rna_add_chain.py',
             'rna_tools/tools/misc/rna_dict2fasta.py',
             'rna_tools/tools/rna_sali2dotbracket/rna_sali2dotbracket.py',

             'rna_tools/tools/rna_rosetta/rna_rosetta_run.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_cluster.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_min.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_n.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_silent_file_sort_and_select.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_extract_lowscore_decoys.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_check_progress.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_head.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_extract.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_get_score.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_silent_file_sort_and_select.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_silent_random.py',
             'rna_tools/tools/rna_rosetta/rna_rosetta_silent_split.py',
             
             'rna_tools/tools/simrna_trajectory/rna_simrna_lowest.py',
             'rna_tools/tools/simrna_trajectory/rna_simrna_extract.py',
             'rna_tools/tools/simrna_trajectory/rna_simrna_cluster.py', 

             'rna_tools/tools/plotting/rna_plot_hist.py',
             'rna_tools/tools/plotting/rna_plot_density.py',
             'rna_tools/tools/plotting/rna_plot_heatmap.py',
             'rna_tools/tools/plotting/rna_plot_dendogram.py',             
             'rna_tools/tools/plotting/rna_plot_boxplotlike.py',             

             'rna_tools/tools/simrna_trajectory/rna_simrna_lowest.py',
             'rna_tools/tools/simrna_trajectory/rna_simrna_extract.py',
             'rna_tools/tools/simrna_trajectory/rna_simrna_cluster.py',

             'rna_tools/tools/rna_filter/rna_filter.py',
             'rna_tools/tools/rna_filter/rna_ec2x.py',
             'rna_tools/tools/rna_filter/rna_pairs2SimRNArestrs.py',
             'rna_tools/tools/rna_filter/rna_ss_get_bps.py',
             'rna_tools/tools/rna_filter/rna_pairs_diff.py',

             'rna_tools/tools/rna_rfam/rna_download_rfam.py',
             
             'rna_tools/tools/rna_alignment/rna_align_coverage.py',
             'rna_tools/tools/rna_alignment/rna_align_seq_to_alignment.py',
             'rna_tools/tools/rna_alignment/rna_align_find_core.py',
             'rna_tools/tools/rna_alignment/rna_align_get_ss_from_stk.py',
             'rna_tools/tools/rna_alignment/rna_align_seq_to_alignment.py',
             'rna_tools/tools/rna_alignment/rna_align_fetch_seed.py',
             'rna_tools/tools/rna_alignment/rna_align_fetch_cm.py',             
             'rna_tools/tools/rna_alignment/utils/rna_alignment_get_species.py',
             'rna_tools/tools/rna_alignment/utils/rna_alignment_process_id.py',
             'rna_tools/tools/rna_alignment/utils/rna_alignment_r2r.py',
             'rna_tools/tools/rna_alignment/rna_align_find_seq_in_alignment.py',
             'rna_tools/tools/rna_alignment/rna_align_foldability.py',
             'rna_tools/tools/rna_pdb_edit_occupancy_bfactor/rna_pdb_edit_occupancy_bfactor.py',
             'rna_tools/tools/rna_refinement/rna_refinement.py',
             'rna_tools/rna_simrnaweb_download_job.py',
             'rna_tools/tools/rna_calc_evo_rmsd/rna_calc_evo_rmsd.py',
             'rna_tools/tools/rna_calc_rmsd/rna_calc_rmsd_multi_targets.py',
             'rna_tools/tools/rna_calc_rmsd_trafl/rna_calc_rmsd_trafl.py',
             'rna_tools/tools/rna_calc_rmsd_trafl/rna_cal_rmsd_trafl_plot.py',
             'rna_tools/tools/rna_x3dna/rna_x3dna.py',
             'rna_tools/tools/renum_pdb_to_aln/renum_pdb_to_aln.py',
             'rna_tools/tools/pdbs_measure_atom_dists/pdbs_measure_atom_dists.py',
             'rna_tools/tools/rna_alignment/random_assignment_of_nucleotides.py',
             'rna_tools/tools/rna_seq_search_BLASTn_outfmt-6/select_seq_fromBLAStn_6outfm.py',

             'rna_tools/tools/misc/translate.py',
             'rna_tools/tools/misc/reverse.py',
             'rna_tools/tools/misc/reverse.py',

             'rna_tools/tools/rna_alignment/utils/rna_align_strip_stk.py',
             'rna_tools/tools/pymol_color_by_conserv/pymol_color_by_conserv.py',
             'rna_tools/tools/pymol_preview_generator/pymol_preview_generator.py',
             'rna_tools/tools/pymol_preview_generator/pymol_preview_install.py',

             'rna_tools/tools/rna_binding_affinity/rna_dG2kd.py',
             'rna_tools/tools/rna_binding_affinity/rna_dG2kd.py',
             'rna_tools/tools/rna_binding_affinity/rna_kd2dG.py',
             'rna_tools/tools/rna_binding_affinity/rna_kd2pkd.py',
             'rna_tools/tools/rna_binding_affinity/rna_pkd2dG.py',
             'rna_tools/tools/rna_binding_affinity/rna_pkd2kd.py',

             'rna_tools/tools/clarna_play/rna_clarna_compare.py',
             'rna_tools/tools/clarna_play/rna_clarna_run.py',
             'rna_tools/tools/clarna_play/ClaRNAlib/clarna.py',

             'rna_tools/tools/mq/rna_mq_score.py',
             'rna_tools/tools/mq/rna_mq_collect.py',             
             'rna_tools/tools/mq/rna_mq_exp.py',
             'rna_tools/tools/mq/rna_mq_plot.py',
             'rna_tools/tools/mq/rna_mq_remote.sh',
             'rna_tools/tools/mq/rna_mq_wx_collect.py',

             'rna_tools/tools/mq/RNA3DCNN/rna_mq_rna3dcnn.py',
             'rna_tools/tools/mq/Dfire/rna_mq_dfire.py',
             'rna_tools/tools/mq/ClashScore/rna_mq_clashscore.py',
             'rna_tools/tools/mq/AnalyzeGeometry/rna_mq_analyzegeometry.py',
             'rna_tools/tools/mq/RNAkb/rna_mq_rnakb.py',
             'rna_tools/tools/mq/rsRNASP/rna_mq_rsRNASP.py',
             'rna_tools/tools/mq/eSCORE/rna_mq_eSCORE.py',
             'rna_tools/tools/mq/FARNA/rna_mq_farna.py',
             'rna_tools/tools/mq/FARFAR2/rna_mq_farfar2.py',
             'rna_tools/tools/mq/SimRNA/rna_mq_simrna.py',
             'rna_tools/tools/mq/RASP/rna_mq_rasp.py',
             
             'rna_tools/tools/md/rna_gromacs_ready.py',
             'rna_tools/tools/md/rna_md.py',
             'rna_tools/tools/md/rna_minimize.py',
             'rna_tools/tools/md/rna_traj_n.py',
             'rna_tools/tools/md/rna_traj_range.py',
             'rna_tools/tools/md/rna_traj_nohetm.py',
             
             'rna_tools/tools/triplexibility/triplexibility.py',
             'rna_tools/tools/triplexibility/triplexibility2.py',

             'rna_tools/tools/triplexibility/trx_mutagen_allvsall.py',
             'rna_tools/tools/triplexibility/trx_mutagen.py',
             'rna_tools/tools/triplexibility/trx_score_allvsall.py',
             'rna_tools/tools/triplexibility/trx_score.py',
             'rna_tools/tools/triplexibility/trx_score_onCol.py',
             
             'rna_tools/tools/df/rna_df_csv_sort.py',
             'rna_tools/tools/df/rna_df_concat.py',
             'rna_tools/tools/df/rna_df_merge.py',
             'rna_tools/tools/df/rna_df_merge_on.py',
             'rna_tools/tools/df/rna_df_merge_two.py',
             
             'rna_tools/tools/misc/rna_ok_mail.py',
             'rna_tools/tools/misc/rna_ok',
             'rna_tools/tools/misc/rna_test_cuda.py',
             
             'rna_tools/tools/PyMOL4RNA/rna_draw_edges.py',

             'rna_tools/tools/rna_calc_hbonds/rna_calc_hbonds.py',

             'rna_tools/tools/rna_pdb_merge_structure_with_fragments/rna_pdb_merge_structure_with_fragments.py',
             'rna_tools/tools/clarna_app/rna_clarna2graph.py',
             ],

    license='GPLv3',
    author='Marcin Magnus',
    author_email='mag_dex@o2.pl',
    description='a toolbox to analyze structures and simulations of RNA',
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    install_requires=[
        'numpy',
        'biopython',
        'progressbar2',
        'tqdm',
        'icecream',
        'csvsort',
        # keep pip minimal
        #'sphinx==1.6.7',
        #'sphinx-argparse==0.1.15',
        #'sphinx-rtd-theme',
        #'sphinxcontrib-napoleon',
        #'forgi',
        #'sphinxcontrib-autoprogram'
        'simplejson',
        'urllib3',
        'wget',

        #'pandas', # rna_calc
        #'matplotlib',
        #'scipy',
        #'python-Levenshtein',

        'versioneer',
        'configparser',
       
      ],
)

