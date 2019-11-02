#!/bin/bash
set -v
export PYTHONPATH=$PYTHONPATH:$(pwd)

cd rna_tools

python update_readme.py

./rna_tools_lib.py

./rna_pdb_toolsx.py -h | tee rna_pdb_tools.out

./rna_pdb_toolsx.py -r input/1xjr.pdb

./rna_pdb_toolsx.py --no-hr -c input/1xjr.pdb > output/1xjr_clx.pdb

##  --rosetta2generic
./rna_pdb_toolsx.py --no-hr --rosetta2generic input/farna.pdb > output/farna_clx.pdb

## --get_chain
./rna_pdb_toolsx.py --get-chain A input/1xjr.pdb > output/1xjr_A_clx.pdb

## extract
./rna_pdb_toolsx.py --no-hr --extract A:1-2 input/pistol_thrs0.50A_clust99-000001_AA.pdb > output/pistol_thrs0.50A_clust99-000001_AA_A-2.pdb

## -rename-chain
./rna_pdb_toolsx.py --no-hr --rename-chain 'B>A' input/8_Chen_10_rpr_ChainB.pdb output/8_Chen_10_rpr_ChainA.pdb

## --clean
./rna_pdb_toolsx.py --no-hr --clean input/1a9l_NMR_1_2_models.pdb > output/1a9l_NMR_1_2_models_tool.pdb
./rna_pdb_toolsx.py --no-hr --clean input/1xjr_GTP.pdb > output/1xjr_GTP.pdb
./rna_pdb_toolsx.py --no-hr --clean input/1osw.pdb > output/1osw_NMR_1.pdb

## --get-rnapuzzle-ready
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/1xjr_onlyGTP.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/1xjr_onlyGTP.pdb > output/1xjr_onlyGTP_rnapuzzle_ready.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/1_das_1_rpr_fixed.pdb > output/1_das_1_rpr_fixed.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/4GXY_3firstNt.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/gtp.pdb  > output/gtp.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle input/377D.pdb # should finish with error
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/rp13_Dokholyan_1_URI_CYT_ADE_GUA_hydrogens.pdb > output/rp13_Dokholyan_1_URI_CYT_ADE_GUA_hydrogens_rpr.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/7_Chen_2_rpr.pdb > output/7_Chen_2_rpr.pdb
./rna_pdb_toolsx.py --no-hr --rpr input/7_Chen_7_rpr.pdb > output/7_Chen_7_rpr.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready input/1I9V_YG_HETATM_ATOM.pdb > output/1I9V_YG_HETATM_ATOM_rpr.pdb
./rna_pdb_toolsx.py --no-hr --get-rnapuzzle-ready --replace-hetatm input/1I9V_A.pdb > output/1I9V_A_rpr.pdb

./rna_pdb_toolsx.py --no-hr --rpr input/A_YG_A.pdb --renumber-residues > output/A_YG_A_renumbered.pdb

./rna_pdb_toolsx.py --no-hr --rpr input/A_YG_A.pdb > output/A_YG_A_norenumbered.pdb

## CCC.pdb
./rna_pdb_toolsx.py --no-hr --rpr --replace-hetatm input/CCC.pdb > output/CCC.pdb
./rna_pdb_toolsx.py --no-hr --rpr input/CCC.pdb > output/CCC_empty.pdb

#./rna_pdb_toolsx.py --no-hr --rpr input/4W90_protein_rna.pdb > output/4W90_protein_rna_onlyRNA.pdb

## --get-rnapuzzle-ready and --dont_rename_chains
./rna_pdb_toolsx.py --no-hr --rpr input/1osw_nt1-4_ChainRenamedToBC.pdb > output/1osw_nt1-4_ChainRenamedToBC_toAB.pdb
./rna_pdb_toolsx.py --no-hr --rpr --dont-rename-chains input/1osw_nt1-4_ChainRenamedToBC.pdb > output/1osw_nt1-4_ChainRenamedToBC_dontRenameChains.pdb
./rna_pdb_toolsx.py --no-hr --swap-chain 'B>A' input/20_Bujnicki_1_mini_2chains_on_TER.pdb > output/20_Bujnicki_1_mini_2chains_on_TER_swapped_chains.pdb
./rna_pdb_toolsx.py --no-hr --swap-chain 'B>B' input/20_Bujnicki_1_mini_2chains_on_TER.pdb > output/20_Bujnicki_1_mini_2chains_on_TER_swapped_chains_B2B.pdb

# --rpr inplace fix
cp input/7_Chen_7_rpr.pdb output/7_Chen_7_rpr_inplacefix.pdb
./rna_pdb_toolsx.py --no-hr --rpr output/7_Chen_7_rpr_inplacefix.pdb --inplace

## --delete
./rna_pdb_toolsx.py --no-hr --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb
./rna_pdb_toolsx.py --no-hr --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb

## --edit
./rna_pdb_toolsx.py --no-hr --edit 'A:6>B:200' input/tetraloop.pdb > output/tetraloop_a6_b200.pdb
./rna_pdb_toolsx.py --no-hr --edit 'A:1-5>B:200-204' input/tetraloop.pdb > output/tetraloop_a1-b200-204.pdb
./rna_pdb_toolsx.py --no-hr --edit 'A:2672>A:1' input/1msy_A2672.pdb > output/1msy_A1.pdb

## --get_seq
./rna_pdb_toolsx.py --get-seq input/5k7c.pdb > output/get_seq.txt
./rna_pdb_toolsx.py --get-seq input/tetraloop.pdb >> output/get_seq.txt
./rna_pdb_toolsx.py --get-seq input/1xjr.pdb > output/1xjr.seq

./rna_pdb_toolsx.py --get-seq input/2_bujnicki_1_rpr.pdb > output/2_bujnicki_1_rpr.txt
./rna_pdb_toolsx.py --get-seq input/2_bujnicki_1_rpr_BA_chain_swap.pdb > output/2_bujnicki_1_rpr_BA_chain_swap.txt

./rna_pdb_toolsx.py --get-seq --uniq '[:5]' --compact input/pistol* > output/pistol_compact.txt
# pistol_thrs0.50A_clust99-000001_AA
#CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGGUGCAA # A:1-62

./rna_pdb_toolsx.py --get-seq --uniq '[:5]' input/pistol* > output/pistol_uniq.txt
# pistol_thrs0.50A_clust99-000001_AA
#CGUGGUUAGGGCCACGUUAAAUAGUUGCUUAAGCCCUAAGCGUUGAUAAAUAUCAGGUGCAA # A:1-62

## --replace_hetatm
./rna_pdb_toolsx.py --no-hr --rpr input/A_YG_A.pdb --replace_hetatm  > output/A_YG_A.pdb

## --get_ss
./rna_pdb_toolsx.py --get-ss input/1xjr*.pdb > output/secondary_structures.txt
# off for now
#./rna_pdb_seq.py input/1ykq_clx.pdb > output/1ykq_clx.seq
#./rna_pdb_seq.py input/1xjr.pdb > output/1xjr2.seq
#./rna_pdb_seq.py input/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb > output/5k7c_clean_onechain_renumber_as_puzzle_srr.seq
#./rna_pdb_seq.py input/6_solution_0.pdb > output/6_solution_0.seq
./rna_pdb_toolsx.py --get-seq input/1ykq_clx.pdb  > output/1ykq_clx.seq
./rna_pdb_toolsx.py --get-seq input/1xjr.pdb > output/1xjr2.seq
./rna_pdb_toolsx.py --get-seq input/5k7c_clean_onechain_renumber_as_puzzle_srr.pdb > output/5k7c_clean_onechain_renumber_as_puzzle_srr.seq
./rna_pdb_toolsx.py --get-seq input/6_solution_0.pdb > output/6_solution_0.seq

## --orgmode
./rna_pdb_toolsx.py --orgmode input/2_das_1_rpr.pdb > output/2_das_1_rpr.org

## --replace-chain
./rna_pdb_toolsx.py --replace-chain output/205d_rmH2o_mutant_A.pdb input/205d_rmH2o.pdb

## --mutate
./rna_pdb_toolsx.py --mutate 'A:1A+2A+3A+4A' input/205d_rmH2o.pdb > output/205d_rmH2o_mutA1234.pdb
cp input/205d_rmH2o.pdb output/205d_rmH2o_mutA1234-B1_inplace.pdb && ./rna_pdb_toolsx.py --mutate 'A:1A+2A+3A+4A,B:13A' --inplace output/205d_rmH2o_mutA1234-B1_inplace.pdb
./rna_pdb_toolsx.py --mutate 'A:1A+2A+3A+4A,B:13A' input/205d_rmH2o.pdb > output/205d_rmH2o_mutA1234-B1.pdb

## --is_pdb
./rna_pdb_toolsx.py --is-pdb input/1I9V_A.pdb
./rna_pdb_toolsx.py --is-pdb input/image.png
./rna_pdb_toolsx.py --is-pdb input/image.png.zip

./rna_pdb_toolsx.py --no-hr --delete-anisou input/3xmg-ANISOU.pdb > output/3xmg-ANISOU.pdb
./rna_pdb_toolsx.py --no-hr --split-alt-locations input/3xmg-polim_ATL_LOC.pdb > output/3xmg-polim_ATL_LOC.pdb

## --renumber_residues
./rna_pdb_toolsx.py --no-hr --renumber-residues input/rp03_solution.pdb > output/rp03_solution_renumber.pdb

./rna_calc_fenergy.py --file input/u2/* --cstinfile --method all  > output/rna_calc_fenergy.txt

./BlastPDB.py

if [ "$1" == "--full" ]; then
    ./RfamSearch.py
fi

./Seq.py

if [ "$1" == "--full" ]; then
    ./SecondaryStructure.py
fi

# ClashCalc
cd ./tools/ClashCalc/
./ClashCalc.py
cd ../..

cd ./tools/rna_calc_rmsd/
./test.sh
cd ../..

cd ./tools/rnashape2ascii/
./test.sh
cd ../..

cd ./tools/rna_pdb_edit_occupancy_bfactor
./test.sh
cd ../..

cd ./tools/rna_filter/
./test.sh
cd ../..

cd ./tools/renum_pdb_to_aln/
./test.sh
cd ../..

if [ "$1" == "--full" ]; then
    cd ./tools/rna_refinement/
    ./test.sh
    cd ../..
fi

if [ "$1" == "--full" ]; then
   cd ..
   codecov --token=e78310dd-7a28-4837-98ef-c93533a84c5b
fi
