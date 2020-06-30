#!/bin/zsh

# created with git clone rna-tools rna-tools-pip
cd ..
git clone rna-tools rna-tools-pip
cd rna-tools-pip # run this script from the root package level
git pull --tags

# strip some data for pip package to keep it under 60 MB
trash rna_tools/input
trash dist
trash build
trash rna_tools/output
trash U6MolCell
trash notes
trash docs
trash rna_tools/tools/spotifier/imgs/

python setup.py bdist_wheel


# py3
#unset PYTHONPATH
#conda activate base
#python setup.py bdist_wheel

twine upload dist/* --verbose
