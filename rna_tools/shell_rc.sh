# add this to you .bashrc or .zshrc to get this functions
# e.g., source ~/work/src/rna-tools/rna_tools/shell_rc.sh

rna_rosetta_progress(){
    parallel 'echo {} && cd {} && rna_rosetta_n.py *.out' ::: `ls -d */`
}

rna_rosetta_run_for_folders(){
 parallel 'cd {} && ./*.sh ' ::: *
}

#
#  shell-aliases
#  install in .bashrc/.zshrc
#
#  source "/home/magnus/src/rna-pdb-tools/rna_pdb_tools/shell-aliases.sh"
#

alias rna_rosetta_n_loop="find . -maxdepth 2 -iname '*out' -type f -not -name 'helix*' -exec rna_rosetta_n.py --verbose {} \; | tee looplog.txt && echo 'sorted' && cat looplog.txt | sort"

if [ "$(uname)" == "Darwin" ]; then
    alias duu="du -h -d 1"
elif
    alias duu='du -s -B G * 2> /dev/null | sort -nr | tee duu' #du -h -d 1'
fi
