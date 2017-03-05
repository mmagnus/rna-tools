#!/bin/bash

export PYTHONPATH=$PYTHONPATH:$(pwd)

cd rna_pdb_tools

./pdb_parser_lib.py

./rna_pdb_tools.py -h | tee rna_pdb_tools.out

./rna_pdb_tools.py -r input/1xjr.pdb 

./rna_pdb_tools.py --no_hr -c input/1xjr.pdb > output/1xjr_clx.pdb

##  --rosetta2generic
./rna_pdb_tools.py --no_hr --rosetta2generic input/farna.pdb > output/farna_clx.pdb

## --get_chain
./rna_pdb_tools.py --get_chain A input/1xjr.pdb > output/1xjr_A_clx.pdb

## --clean
./rna_pdb_tools.py --no_hr --clean input/1a9l_NMR_1_2_models.pdb > output/1a9l_NMR_1_2_models_tool.pdb
./rna_pdb_tools.py --no_hr --clean input/1xjr_GTP.pdb > output/1xjr_GTP.pdb
./rna_pdb_tools.py --no_hr --clean input/1osw.pdb > output/1osw_NMR_1.pdb

## --get_rnapuzzle_ready
./rna_pdb_tools.py --no_hr --get_rnapuzzle_ready input/1xjr_onlyGTP.pdb
./rna_pdb_tools.py --no_hr --get_rnapuzzle_ready input/1xjr_onlyGTP.pdb > output/1xjr_onlyGTP_rnapuzzle_ready.pdb
python rna_pdb_tools.py --get_rnapuzzle_ready input/1_das_1_rpr_fixed.pdb > output/1_das_1_rpr_fixed.pdb
./rna_pdb_tools.py --no_hr --get_rnapuzzle_ready input/4GXY_3firstNt.pdb
./rna_pdb_tools.py --no_hr --get_rnapuzzle_ready input/gtp.pdb  > output/gtp.pdb
./rna_pdb_tools.py --no_hr --get_rnapuzzle input/377D.pdb # should finish with error
./rna_pdb_tools.py --no_hr  --get_rnapuzzle_ready input/rp13_Dokholyan_1_URI_CYT_ADE_GUA_hydrogens.pdb > output/rp13_Dokholyan_1_URI_CYT_ADE_GUA_hydrogens_rpr.pdb

## --delete
./rna_pdb_tools.py --no_hr --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb
./rna_pdb_tools.py --no_hr --delete A:10-60 input/rp17.out.1.pdb > output/rp17_rmA10-60.pdb

## --edit
./rna_pdb_tools.py --no_hr --edit 'A:6>B:200' input/tetraloop.pdb > output/tetraloop_a6_b200.pdb
./rna_pdb_tools.py --no_hr --edit 'A:1-5>B:200-204' input/tetraloop.pdb > output/tetraloop_a1-b200-204.pdb
./rna_pdb_tools.py --no_hr --edit 'A:2672>A:1' input/1msy_A2672.pdb > output/1msy_A1.pdb

## --get_seq
./rna_pdb_tools.py --no_hr --get_seq input/5k7c.pdb > output/get_seq.txt
./rna_pdb_tools.py --no_hr --get_seq input/tetraloop.pdb >> output/get_seq.txt
./rna_pdb_tools.py --get_seq input/1xjr.pdb > output/1xjr.seq

## --renumber_residues
./rna_pdb_tools.py --no_hr --renumber_residues input/rp03_solution.pdb > output/rp03_solution_renumber.pdb

# ClashCalc
cd ./utils/ClashCalc/
./ClashCalc.py
cd ../..

cd ./utils/rna_calc_rmsd/
./test.sh
cd ../..

cd ./utils/rnashape2ascii/
./test.sh
cd ../..

cd ./utils/simrna_trajectory
./test.sh
cd ../..

cd ./utils/rna_filter/
./test.sh
cd ../..

cd ..
codecov --token=e78310dd-7a28-4837-98ef-c93533a84c5b
