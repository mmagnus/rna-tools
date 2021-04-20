import random
import pandas as pd

values = [x / 10.0 for x in range(-10, 11)]
print(values)
all = []
def err(g, c):
    return abs(g - c) # -2
    
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

    if errg == 0:
        print("OK!", t[0], c15)
        break
