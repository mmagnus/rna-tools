set -x
splix_trx.py GCU --edge cWW_cHS #; echo 0 # bad 
# H
#splix_trx.py UAU
#
splix_trx.py GCA --edge cWW_cHS #; echo 0 # bad
splix_trx.py UAU --edge cWW_cHS #; echo 1 # good
splix_trx.py UAU --edge cWW_cHW #; echo 1 # good
            #GCU
