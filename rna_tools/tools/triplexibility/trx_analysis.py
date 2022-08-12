from trx import *
seqs = get_seqs()
inst = 0
exists = 0
for s in seqs:
    #if tdb_score(s):
    if tdb_score(s, type='exists'):
        exists += 1
        if tdb_score(s):
            inst += 1
    print(s, tdb_score(s), tdb_score(s, type='exists'))
print(inst, exists)
