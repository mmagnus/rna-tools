import random
import pandas as pd

values = [x / 10.0 for x in range(-10, 11)]
#print(values)
all = []
def err(g, c):
    return abs(g - c) # -2
    
def log(t, c15, e):
        tx = ''
        if t[0] < 0:
            tx += str(t[3]) + ' inhibits splicing by ' + str(t[0]) + ' and '#.rjust(10)
        if t[0] > 0:
            tx += str(t[3]) + ' promote splicing by ' + str(t[0]) + ' and '#.rjust(12)
        if t[0] == 0:
            tx += str(t[3]) + ' does nothing to '#.rjust(10)

        if c15 < 0:
            tx += 'cwc15 inhibits splicing by ' + str(c15)
        if c15 > 0:
            tx += 'cwc15 promote splicing by ' + str(c15)
        if c15 == 0:
            tx += 'cwc15 does nothing to '
        tx += ' with error ' + str(round(e, 2))
        print(tx)


for i in range(0, 100000000):
    done = False
    ggg = (random.choice(values), 0, -2, 'ggg')
    c15 = random.choice(values)
    triples = [ggg]
    errg = 0
    for t in triples:
        c = t[0] + c15 # -1 + 1 = 0
        g = t[1]
        e1 = err(g, c)
        
        c2 = t[0] - c15 # -1 - +1 = -2
        g2 = t[2]
        e2 = err(g2, c2)
        # print(t[3], t[0], c15, e1, e2)
        errg += e1 + e2

        log(t, c15, e1 + e2)
        
    if errg == 0:
        print('OK!')
        log(t, c15, errg)
        break
