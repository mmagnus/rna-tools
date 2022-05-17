# add this to you .bashrc or .zshrc to get this functions
# e.g., source ~/work/src/rna-tools/rna_tools/shell_rc.sh

rna_rosetta_progress(){
    parallel 'echo {} && cd {} && rna_rosetta_n.py *.out' ::: `ls -d */`
}

rna_rosetta_run_for_folders(){
    parallel 'cd {} && ./*.sh ' ::: *
}


trx(){
    rna_rosetta_extract_lowscore_decoys.py 10 *out
    parallel rna_mq_farfar2.py -qr {} ::: *out*.pdb | tee ff2.csv
}

alias ta="tmux attach"
