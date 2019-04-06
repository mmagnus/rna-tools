#!/bin/bash

cmd="python select_seq_fromBLAStn_6outfm.py  test_data/db_nt.fasta test_data/output_55 --outfn  fasta_file"
echo $cmd
$cmd
