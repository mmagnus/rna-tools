#!/bin/zsh

git push --tags
# created with git clone rna-tools rna-tools-pip
cd /Users/magnus/work/src # ugly
pwd
trash rna-tools-pip
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

# build wheel and sdist via the pinned build backend (see pyproject.toml)
python -m build --wheel --sdist


# py3
#unset PYTHONPATH
#conda activate base
#python setup.py bdist_wheel

twine upload dist/* --verbose

# clean for mdfind
cd ..
trash rna-tools-pip
trash rna-tools/build
