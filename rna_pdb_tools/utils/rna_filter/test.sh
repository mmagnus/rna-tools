#!/bin/bash

cmd="python rna_filter.py -r restraints.txt -s test_data/CG.pdb -v"
echo $cmd
$cmd
echo
cmd="python rna_filter.py -r restraints.txt -t test_data/CG.trafl -v"
echo $cmd
$cmd

