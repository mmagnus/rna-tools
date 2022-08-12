from trx import *
seqs = get_seqs()
inst = 0
exists = 0

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='> ')
p = ic

seq = get_seqs()
i = 0
ii = 64
p = ic

for s in seq:
    if tdb_score(s) > 0:
      i += 1
    if tdb_score(s) == -1:
      ii -= 1
p(i, ii)

ic(tdb_score('gcu'), tdb_score('cgc'), tdb_score('ggg'))
# growth! here
# based on rmsd?

sys.exit(1)

for s in seqs:
    #if tdb_score(s):
    if tdb_score(s, type='exists'):
        exists += 1
        if tdb_score(s):
            inst += 1
    print(s, tdb_score(s), tdb_score(s, type='exists'))
print(inst, exists)
