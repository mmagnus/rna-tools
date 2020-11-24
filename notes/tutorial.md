# rna-tools

    3.6.x

# rna-tools & PyMOL

```
fetch 4ts2
inspect 4ts2
rpr 4ts2
clarna 4ts2
clarna (sele) # clarna sele

    PyMOL>clarna (sele)
    rna_clarna_run.py -ipdb /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpfs3hbhcw.pdb -bp+stack
    chains:  X 15 16
    X   15   X   16          bp A A                      ><   0.9427 
    # stacking interaction
    
```

# Docs

https://rna-tools.readthedocs.io/en/latest/install.html

# Install

    pip install rna-tools
    python -c 'import rna_tools'

or

    pip install -e git+http://github.com/mmagnus/rna-tools.git#egg=rna-tools

# RNAseq

    $ rna_secondary_structure_prediction.py --method cyclefold --seq ggCUUccGUAAggAUGcc
    ggCUUccGUAAggAUGc

    $ rna_secondary_structure_prediction.py --method mcfold --seq ggCUUccGUAAggAUGcc
    ggCUUccGUAAggAUGcc
    ((((((((..)))))))) -20.54

--verbose:

    na_secondary_structure_prediction.py --method cyclefold --seq ggCUUccGUAAggAUGcc --verbose
    ggCUUccGUAAggAUGcc
    /Users/magnus/work/opt/RNAstructure/6.1//exe/CycleFold /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmp3xqx9ami.fa > /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmp3xqx9ami.fa.ct
    /Users/magnus/work/opt/RNAstructure/6.1//exe/ct2dot /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmp3xqx9ami.fa.ct 1 /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmp3xqx9ami.fa.dot
    ((((((((..)))))))) -20.5493


# Calc RMSD

# Calc INFs

# RNAalignment

# Misc

    
    $ rna_pdb_toolsx.py --fetch 5lj3
    /usr/local/lib/python2.7/site-packages/urllib3/connectionpool.py:858: InsecureRequestWarning: Unverified HTTPS request is being made. Adding certificate verification is strongly advised. See: https://urllib3.readthedocs.io/en/latest/advanced-usage.html#ssl-warnings
      InsecureRequestWarning)
    downloading..../5lj3.pdb
    ok

# Notes

```
PyMOL>clarna
chains:  A 4 7
A    4   A    7   bp G A SH_tran 0.72
```
