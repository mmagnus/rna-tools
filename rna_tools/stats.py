G_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split()
A_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split()
U_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split()
C_ATOMS = "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split()

d = {'P' : 0,
     'O' : 0,
     'C' : 0,
     'N' : 0}
for i in G_ATOMS + A_ATOMS + U_ATOMS + C_ATOMS:
    d[i[0]] += 1
print(d)
    
