#
#  shell-aliases
#  install in .bashrc/.zshrc
#
#  source "/home/magnus/src/rna-pdb-tools/rna_pdb_tools/shell-aliases.sh"
#
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

alias rna_rosetta_n_loop="find . -maxdepth 2 -iname '*out' -type f -not -name 'helix*' -exec rna_rosetta_n.py --verbose {} \; | tee looplog.txt && echo 'sorted' && cat looplog.txt | sort"

setopt NO_EQUALS
if [ "$(uname)" == "Darwin" ]; then
    alias duu="du -h -d 1"
else
    alias duu='du -s -B G * 2> /dev/null | sort -nr | tee duu' #du -h -d 1'
fi

xx(){
    set -x
    for 1 in $@; do
       time 7zz a -m0=lzma2 -mmt=on ${1}.lzma2.7z ${1} && trash ${1} # mv -v ${1} ${1}_7zed
    done;
}

x9(){
    set -x
    for 1 in $@; do
	     time 7zz a  -mx=9 -m0=lzma2 -mmt=on ${1}.`date +%y%m%d`.lzma2.7z ${1} && trash ${1} # mv -v ${1} ${1}_7zed
    done;

    #
    #
    ## | tee ${1}.log
    # zipinfo ${1}.7z >> ${1}.log
    # trash ${1}
}
