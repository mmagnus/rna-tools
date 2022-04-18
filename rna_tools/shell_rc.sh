# add this to you .bashrc or .zshrc to get this functions
# e.g., source ~/work/src/rna-tools/rna_tools/shell_rc.sh

rna_rosetta_progress(){
    parallel 'echo {} && cd {} && rna_rosetta_n.py default.out' ::: `ls -d */`
}

rna_rosetta_run_for_folders(){
 parallel 'cd {} && ./*.sh ' ::: *
}
