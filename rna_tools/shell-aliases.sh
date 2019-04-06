#
#  shell-aliases
#  install in .bashrc/.zshrc
#
#  source "/home/magnus/src/rna-pdb-tools/rna_pdb_tools/shell-aliases.sh"
#

alias rna_rosetta_n_loop="find . -maxdepth 2 -iname '*out' -type f -not -name 'helix*' -exec rna_rosetta_n.py --verbose {} \; | tee looplog.txt && echo 'sorted' && cat looplog.txt | sort"
