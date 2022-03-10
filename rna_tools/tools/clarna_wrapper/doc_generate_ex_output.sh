#!/bin/bash

exe='./clarna_wrapper.py test_data/1mme.pdb'
echo "exe: $exe" > examplary_output.txt

$exe | tee -a examplary_output.txt
