
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

def spl_color():
 for name in cmd.get_names("all"):
    # cmd.do('color grey50') # off gray
    #if name in ['5zwo', '5gm6', '5lj3', '5mps', '6exn', '5ylz', '5y88', '3jb9', '6icz', '6ff7', '5yzg', '5xjc', '5gan', '6qw6', '3jcr', '6qx9', '6ah0']: # this should be auto
    print(" \ Extracting mode for %s" % name)
    
    for n in cmd.get_names("all"):
        if '5zwo' in n:
            object_name = n
            print(object_name)
            
    if '5zwo' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_B5zwo")
        #BRR2
        cmd.do("color grey60, chain D and " + object_name)
        cmd.do("color grey60, BRR2_B5zwo")
        #LEA1
        cmd.do("color palegreen, chain o and " + object_name)
        cmd.do("color palegreen, LEA1_B5zwo")
        #Msl1
        cmd.do("color palegreen, chain p and " + object_name)
        cmd.do("color palegreen, Msl1_B5zwo")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_B5zwo")
        #U2
        cmd.do("color forest, chain H and " + object_name)
        cmd.do("color forest, U2_B5zwo")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_B5zwo")
        #U6
        cmd.do("color firebrick, chain F and " + object_name)
        cmd.do("color firebrick, U6_B5zwo")
        #U4
        cmd.do("color brown, chain I and " + object_name)
        cmd.do("color brown, U4_B5zwo")
        #Intron
        cmd.do("color black, chain G and " + object_name)
        cmd.do("color black, Intron_B5zwo")
        #PRP4
        cmd.do("color grey50, chain K and " + object_name)
        cmd.do("color grey50, PRP4_B5zwo")
        #PRP31
        cmd.do("color grey50, chain L and " + object_name)
        cmd.do("color grey50, PRP31_B5zwo")
        #PRP6
        cmd.do("color grey50, chain N and " + object_name)
        cmd.do("color grey50, PRP6_B5zwo")
        #PRP3
        cmd.do("color grey50, chain J and " + object_name)
        cmd.do("color grey50, PRP3_B5zwo")
        #DIB1
        cmd.do("color grey50, chain E and " + object_name)
        cmd.do("color grey50, DIB1_B5zwo")
        #SNU13
        cmd.do("color grey50, chain M and " + object_name)
        cmd.do("color grey50, SNU13_B5zwo")
        #LSM8
        cmd.do("color grey50, chain z and " + object_name)
        cmd.do("color grey50, LSM8_B5zwo")
        #LSM2
        cmd.do("color grey50, chain q and " + object_name)
        cmd.do("color grey50, LSM2_B5zwo")
        #LSM3
        cmd.do("color grey50, chain r and " + object_name)
        cmd.do("color grey50, LSM3_B5zwo")
        #LSM6
        cmd.do("color grey50, chain x and " + object_name)
        cmd.do("color grey50, LSM6_B5zwo")
        #LSM5
        cmd.do("color grey50, chain t and " + object_name)
        cmd.do("color grey50, LSM5_B5zwo")
        #LSM7
        cmd.do("color grey50, chain y and " + object_name)
        cmd.do("color grey50, LSM7_B5zwo")
        #LSM4
        cmd.do("color grey50, chain s and " + object_name)
        cmd.do("color grey50, LSM4_B5zwo")
        #SNU66
        cmd.do("color grey50, chain O and " + object_name)
        cmd.do("color grey50, SNU66_B5zwo")
        #BUD13
        cmd.do("color grey60, chain Y and " + object_name)
        cmd.do("color grey60, BUD13_B5zwo")
        #Cus1
        cmd.do("color palegreen, chain 2 and " + object_name)
        cmd.do("color palegreen, Cus1_B5zwo")
        #HSH155
        cmd.do("color smudge, chain 1 and " + object_name)
        cmd.do("color smudge, HSH155_B5zwo")
        #HSH49
        cmd.do("color sand, chain 4 and " + object_name)
        cmd.do("color sand, HSH49_B5zwo")
        #PML1
        cmd.do("color grey60, chain Z and " + object_name)
        cmd.do("color grey60, PML1_B5zwo")
        #PRP11
        cmd.do("color palegreen, chain v and " + object_name)
        cmd.do("color palegreen, PRP11_B5zwo")
        #RDS3
        cmd.do("color palegreen, chain 5 and " + object_name)
        cmd.do("color palegreen, RDS3_B5zwo")
        #RSE1
        cmd.do("color smudge, chain 3 and " + object_name)
        cmd.do("color smudge, RSE1_B5zwo")
        #SNU17
        cmd.do("color grey60, chain X and " + object_name)
        cmd.do("color grey60, SNU17_B5zwo")
        #Ysf3
        cmd.do("color palegreen, chain 6 and " + object_name)
        cmd.do("color palegreen, Ysf3_B5zwo")
        cmd.do("color grey50, chain a and 5zwo")
        cmd.do("color grey50, chain P and 5zwo")
        cmd.do("color grey50, chain h and 5zwo")
        cmd.do("color grey50, chain e and 5zwo")
        cmd.do("color grey50, chain T and 5zwo")
        cmd.do("color grey50, chain i and 5zwo")
        cmd.do("color grey50, chain f and 5zwo")
        cmd.do("color grey50, chain U and 5zwo")
        cmd.do("color grey50, chain j and 5zwo")
        cmd.do("color grey50, chain g and 5zwo")
        cmd.do("color grey50, chain V and 5zwo")
        cmd.do("color grey50, chain k and 5zwo")
        cmd.do("color grey50, chain d and 5zwo")
        cmd.do("color grey50, chain S and 5zwo")
        cmd.do("color grey50, chain l and 5zwo")
        cmd.do("color grey50, chain b and 5zwo")
        cmd.do("color grey50, chain Q and 5zwo")
        cmd.do("color grey50, chain m and 5zwo")
        cmd.do("color grey50, chain c and 5zwo")
        cmd.do("color grey50, chain R and 5zwo")
        cmd.do("color grey50, chain n and 5zwo")
        #PRP9
        cmd.do("color grey50, chain u and " + object_name)
        cmd.do("color grey50, PRP9_B5zwo")
        #PRP21
        cmd.do("color grey50, chain w and " + object_name)
        cmd.do("color grey50, PRP21_B5zwo")
        #SNU23
        cmd.do("color grey50, chain W and " + object_name)
        cmd.do("color grey50, SNU23_B5zwo")
        #PRP38
        cmd.do("color grey50, chain 0 and " + object_name)
        cmd.do("color grey50, PRP38_B5zwo")
        #SPP381
        cmd.do("color grey50, chain 9 and " + object_name)
        cmd.do("color grey50, SPP381_B5zwo")
    
    for n in cmd.get_names("all"):
        if '5gm6' in n:
            object_name = n
            print(object_name)
            
    if '5gm6' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_Ba5gm6")
        #BRR2
        cmd.do("color grey60, chain B and " + object_name)
        cmd.do("color grey60, BRR2_Ba5gm6")
        #BUD31
        cmd.do("color dirtyviolet, chain T and " + object_name)
        cmd.do("color dirtyviolet, BUD31_Ba5gm6")
        #CEF1
        cmd.do("color raspberry, chain c and " + object_name)
        cmd.do("color raspberry, CEF1_Ba5gm6")
        #CWC15
        cmd.do("color orange, chain S and " + object_name)
        cmd.do("color orange, CWC15_Ba5gm6")
        #CWC2_hRBM22
        cmd.do("color ruby, chain R and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_Ba5gm6")
        #CWC21
        cmd.do("color violetpurple, chain X and " + object_name)
        cmd.do("color violetpurple, CWC21_Ba5gm6")
        #CWC22
        cmd.do("color bluewhite, chain Z and " + object_name)
        cmd.do("color bluewhite, CWC22_Ba5gm6")
        #PRP45
        cmd.do("color lightpink, chain P and " + object_name)
        cmd.do("color lightpink, PRP45_Ba5gm6")
        #CDC40
        cmd.do("color dirtyviolet, chain n and " + object_name)
        cmd.do("color dirtyviolet, CDC40_Ba5gm6")
        #PRP19
        cmd.do("color grey70, chain f and " + object_name)
        cmd.do("color grey70, PRP19_Ba5gm6")
        #PRP46
        cmd.do("color lightblue, chain O and " + object_name)
        cmd.do("color lightblue, PRP46_Ba5gm6")
        #SLT11/ECM2
        cmd.do("color chocolate, chain Q and " + object_name)
        cmd.do("color chocolate, SLT11/ECM2_Ba5gm6")
        #SNT309
        cmd.do("color grey70, chain t and " + object_name)
        cmd.do("color grey70, SNT309_Ba5gm6")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_Ba5gm6")
        #SYF2
        cmd.do("color brightorange, chain f and " + object_name)
        cmd.do("color brightorange, SYF2_Ba5gm6")
        #SYF1
        cmd.do("color brightorange, chain v and " + object_name)
        cmd.do("color brightorange, SYF1_Ba5gm6")
        #U2
        cmd.do("color forest, chain 2 and " + object_name)
        cmd.do("color forest, U2_Ba5gm6")
        #U5
        cmd.do("color density, chain 5 and " + object_name)
        cmd.do("color density, U5_Ba5gm6")
        #U6
        cmd.do("color firebrick, chain 6 and " + object_name)
        cmd.do("color firebrick, U6_Ba5gm6")
        #Intron
        cmd.do("color black, chain M and " + object_name)
        cmd.do("color black, Intron_Ba5gm6")
        #Exon
        cmd.do("color yellow, chain N and " + object_name)
        cmd.do("color yellow, Exon_Ba5gm6")
        #BUD13
        cmd.do("color grey60, chain W and " + object_name)
        cmd.do("color grey60, BUD13_Ba5gm6")
        #CLF2
        cmd.do("color rasberry, chain d and " + object_name)
        cmd.do("color rasberry, CLF2_Ba5gm6")
        #Cus1
        cmd.do("color palegreen, chain H and " + object_name)
        cmd.do("color palegreen, Cus1_Ba5gm6")
        #CWC24
        cmd.do("color grey60, chain a and " + object_name)
        cmd.do("color grey60, CWC24_Ba5gm6")
        #CWC27
        cmd.do("color grey60, chain b and " + object_name)
        cmd.do("color grey60, CWC27_Ba5gm6")
        #HSH155
        cmd.do("color smudge, chain G and " + object_name)
        cmd.do("color smudge, HSH155_Ba5gm6")
        #HSH49
        cmd.do("color sand, chain e and " + object_name)
        cmd.do("color sand, HSH49_Ba5gm6")
        #PML1
        cmd.do("color grey60, chain U and " + object_name)
        cmd.do("color grey60, PML1_Ba5gm6")
        #PRP11
        cmd.do("color palegreen, chain I and " + object_name)
        cmd.do("color palegreen, PRP11_Ba5gm6")
        #PRP2
        cmd.do("color palegreen, chain Y and " + object_name)
        cmd.do("color palegreen, PRP2_Ba5gm6")
        #RDS3
        cmd.do("color palegreen, chain J and " + object_name)
        cmd.do("color palegreen, RDS3_Ba5gm6")
        #RSE1
        cmd.do("color smudge, chain F and " + object_name)
        cmd.do("color smudge, RSE1_Ba5gm6")
        #SNU17
        cmd.do("color grey60, chain V and " + object_name)
        cmd.do("color grey60, SNU17_Ba5gm6")
        #Ysf3
        cmd.do("color palegreen, chain K and " + object_name)
        cmd.do("color palegreen, Ysf3_Ba5gm6")
    
    for n in cmd.get_names("all"):
        if '5lj3' in n:
            object_name = n
            print(object_name)
            
    if '5lj3' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_C5lj3")
        #BUD31
        cmd.do("color dirtyviolet, chain L and " + object_name)
        cmd.do("color dirtyviolet, BUD31_C5lj3")
        #CEF1
        cmd.do("color raspberry, chain O and " + object_name)
        cmd.do("color raspberry, CEF1_C5lj3")
        #CLF1
        cmd.do("color raspberry, chain S and " + object_name)
        cmd.do("color raspberry, CLF1_C5lj3")
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_C5lj3")
        #CWC16/YJU2
        cmd.do("color lightteal, chain D and " + object_name)
        cmd.do("color lightteal, CWC16/YJU2_C5lj3")
        #CWC2_hRBM22
        cmd.do("color ruby, chain M and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_C5lj3")
        #CWC21
        cmd.do("color violetpurple, chain R and " + object_name)
        cmd.do("color violetpurple, CWC21_C5lj3")
        #CWC22
        cmd.do("color bluewhite, chain H and " + object_name)
        cmd.do("color bluewhite, CWC22_C5lj3")
        #CWC25
        cmd.do("color deepteal, chain F and " + object_name)
        cmd.do("color deepteal, CWC25_C5lj3")
        #ISY1
        cmd.do("color dirtyviolet, chain G and " + object_name)
        cmd.do("color dirtyviolet, ISY1_C5lj3")
        #LEA1
        cmd.do("color palegreen, chain W and " + object_name)
        cmd.do("color palegreen, LEA1_C5lj3")
        #Msl1
        cmd.do("color palegreen, chain Y and " + object_name)
        cmd.do("color palegreen, Msl1_C5lj3")
        #PRP45
        cmd.do("color lightpink, chain K and " + object_name)
        cmd.do("color lightpink, PRP45_C5lj3")
        #PRP46
        cmd.do("color lightblue, chain J and " + object_name)
        cmd.do("color lightblue, PRP46_C5lj3")
        #SLT11/ECM2
        cmd.do("color chocolate, chain N and " + object_name)
        cmd.do("color chocolate, SLT11/ECM2_C5lj3")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_C5lj3")
        #SYF1
        cmd.do("color brightorange, chain T and " + object_name)
        cmd.do("color brightorange, SYF1_C5lj3")
        #U2
        cmd.do("color forest, chain Z and " + object_name)
        cmd.do("color forest, U2_C5lj3")
        #U5
        cmd.do("color density, chain U and " + object_name)
        cmd.do("color density, U5_C5lj3")
        #U6
        cmd.do("color firebrick, chain V and " + object_name)
        cmd.do("color firebrick, U6_C5lj3")
        #Intron
        cmd.do("color black, chain I and " + object_name)
        cmd.do("color black, Intron_C5lj3")
        #Exon
        cmd.do("color yellow, chain E and " + object_name)
        cmd.do("color yellow, Exon_C5lj3")
        cmd.do("color grey50, chain b and 5lj3")
        cmd.do("color grey50, chain k and 5lj3")
        cmd.do("color grey50, chain e and 5lj3")
        cmd.do("color grey50, chain p and 5lj3")
        cmd.do("color grey50, chain f and 5lj3")
        cmd.do("color grey50, chain q and 5lj3")
        cmd.do("color grey50, chain g and 5lj3")
        cmd.do("color grey50, chain r and 5lj3")
        cmd.do("color grey50, chain d and 5lj3")
        cmd.do("color grey50, chain n and 5lj3")
        cmd.do("color grey50, chain h and 5lj3")
        cmd.do("color grey50, chain l and 5lj3")
        cmd.do("color grey50, chain j and 5lj3")
        cmd.do("color grey50, chain m and 5lj3")
    
    for n in cmd.get_names("all"):
        if '5mps' in n:
            object_name = n
            print(object_name)
            
    if '5mps' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_Cs5mps")
        #BUD31
        cmd.do("color dirtyviolet, chain L and " + object_name)
        cmd.do("color dirtyviolet, BUD31_Cs5mps")
        #CEF1
        cmd.do("color raspberry, chain O and " + object_name)
        cmd.do("color raspberry, CEF1_Cs5mps")
        #CLF1
        cmd.do("color raspberry, chain S and " + object_name)
        cmd.do("color raspberry, CLF1_Cs5mps")
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_Cs5mps")
        #CWC2_hRBM22
        cmd.do("color ruby, chain M and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_Cs5mps")
        #CWC21
        cmd.do("color violetpurple, chain R and " + object_name)
        cmd.do("color violetpurple, CWC21_Cs5mps")
        #CWC22
        cmd.do("color bluewhite, chain H and " + object_name)
        cmd.do("color bluewhite, CWC22_Cs5mps")
        #PRP45
        cmd.do("color lightpink, chain K and " + object_name)
        cmd.do("color lightpink, PRP45_Cs5mps")
        #CDC40
        cmd.do("color dirtyviolet, chain o and " + object_name)
        cmd.do("color dirtyviolet, CDC40_Cs5mps")
        #PRP46
        cmd.do("color lightblue, chain J and " + object_name)
        cmd.do("color lightblue, PRP46_Cs5mps")
        #SLT11/ECM2
        cmd.do("color chocolate, chain N and " + object_name)
        cmd.do("color chocolate, SLT11/ECM2_Cs5mps")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_Cs5mps")
        #SYF2
        cmd.do("color brightorange, chain y and " + object_name)
        cmd.do("color brightorange, SYF2_Cs5mps")
        #SYF1
        cmd.do("color brightorange, chain T and " + object_name)
        cmd.do("color brightorange, SYF1_Cs5mps")
        #U2
        cmd.do("color forest, chain 2 and " + object_name)
        cmd.do("color forest, U2_Cs5mps")
        #U5
        cmd.do("color density, chain 5 and " + object_name)
        cmd.do("color density, U5_Cs5mps")
        #U6
        cmd.do("color firebrick, chain 6 and " + object_name)
        cmd.do("color firebrick, U6_Cs5mps")
        #5EXON
        cmd.do("color grey50, chain E and " + object_name)
        cmd.do("color grey50, 5EXON_Cs5mps")
        #Intron
        cmd.do("color black, chain I and " + object_name)
        cmd.do("color black, Intron_Cs5mps")
        #Exon
        cmd.do("color yellow, chain E and " + object_name)
        cmd.do("color yellow, Exon_Cs5mps")
        #SMB1
        cmd.do("color grey50, chain b and " + object_name)
        cmd.do("color grey50, SMB1_Cs5mps")
        #SME1
        cmd.do("color grey50, chain e and " + object_name)
        cmd.do("color grey50, SME1_Cs5mps")
        #SMX3
        cmd.do("color grey50, chain f and " + object_name)
        cmd.do("color grey50, SMX3_Cs5mps")
        #SMX2
        cmd.do("color grey50, chain g and " + object_name)
        cmd.do("color grey50, SMX2_Cs5mps")
        #SMD3
        cmd.do("color grey50, chain d and " + object_name)
        cmd.do("color grey50, SMD3_Cs5mps")
        #SMD1
        cmd.do("color grey50, chain h and " + object_name)
        cmd.do("color grey50, SMD1_Cs5mps")
        #SMD2
        cmd.do("color grey50, chain j and " + object_name)
        cmd.do("color grey50, SMD2_Cs5mps")
        #PRP18
        cmd.do("color grey50, chain a and " + object_name)
        cmd.do("color grey50, PRP18_Cs5mps")
        #SLU7
        cmd.do("color grey50, chain c and " + object_name)
        cmd.do("color grey50, SLU7_Cs5mps")
    
    for n in cmd.get_names("all"):
        if '6exn' in n:
            object_name = n
            print(object_name)
            
    if '6exn' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_P6exn")
        #BUD31
        cmd.do("color dirtyviolet, chain L and " + object_name)
        cmd.do("color dirtyviolet, BUD31_P6exn")
        #CEF1
        cmd.do("color raspberry, chain O and " + object_name)
        cmd.do("color raspberry, CEF1_P6exn")
        #CLF1
        cmd.do("color raspberry, chain S and " + object_name)
        cmd.do("color raspberry, CLF1_P6exn")
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_P6exn")
        #CWC16/YJU2
        cmd.do("color lightteal, chain D and " + object_name)
        cmd.do("color lightteal, CWC16/YJU2_P6exn")
        #CWC2_hRBM22
        cmd.do("color ruby, chain M and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_P6exn")
        #CWC21
        cmd.do("color violetpurple, chain R and " + object_name)
        cmd.do("color violetpurple, CWC21_P6exn")
        #CWC22
        cmd.do("color bluewhite, chain H and " + object_name)
        cmd.do("color bluewhite, CWC22_P6exn")
        #LEA1
        cmd.do("color palegreen, chain W and " + object_name)
        cmd.do("color palegreen, LEA1_P6exn")
        #Msl1
        cmd.do("color palegreen, chain Y and " + object_name)
        cmd.do("color palegreen, Msl1_P6exn")
        #PRP45
        cmd.do("color lightpink, chain K and " + object_name)
        cmd.do("color lightpink, PRP45_P6exn")
        #CDC40
        cmd.do("color dirtyviolet, chain o and " + object_name)
        cmd.do("color dirtyviolet, CDC40_P6exn")
        cmd.do("color grey70, chain t and 6exn")
        cmd.do("color grey70, chain u and 6exn")
        cmd.do("color grey70, chain v and 6exn")
        cmd.do("color grey70, chain w and 6exn")
        #PRP46
        cmd.do("color lightblue, chain J and " + object_name)
        cmd.do("color lightblue, PRP46_P6exn")
        #SLT11/ECM2
        cmd.do("color chocolate, chain N and " + object_name)
        cmd.do("color chocolate, SLT11/ECM2_P6exn")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_P6exn")
        #SYF1
        cmd.do("color brightorange, chain T and " + object_name)
        cmd.do("color brightorange, SYF1_P6exn")
        #U2
        cmd.do("color forest, chain 2 and " + object_name)
        cmd.do("color forest, U2_P6exn")
        #U5
        cmd.do("color density, chain 5 and " + object_name)
        cmd.do("color density, U5_P6exn")
        #U6
        cmd.do("color firebrick, chain 6 and " + object_name)
        cmd.do("color firebrick, U6_P6exn")
        #Intron
        cmd.do("color black, chain I and " + object_name)
        cmd.do("color black, Intron_P6exn")
        #Exon
        cmd.do("color yellow, chain E and " + object_name)
        cmd.do("color yellow, Exon_P6exn")
        cmd.do("color grey50, chain b and 6exn")
        cmd.do("color grey50, chain k and 6exn")
        cmd.do("color grey50, chain e and 6exn")
        cmd.do("color grey50, chain p and 6exn")
        cmd.do("color grey50, chain f and 6exn")
        cmd.do("color grey50, chain q and 6exn")
        cmd.do("color grey50, chain g and 6exn")
        cmd.do("color grey50, chain r and 6exn")
        cmd.do("color grey50, chain d and 6exn")
        cmd.do("color grey50, chain n and 6exn")
        cmd.do("color grey50, chain h and 6exn")
        cmd.do("color grey50, chain l and 6exn")
        cmd.do("color grey50, chain j and 6exn")
        cmd.do("color grey50, chain m and 6exn")
        #PRP22
        cmd.do("color grey50, chain V and " + object_name)
        cmd.do("color grey50, PRP22_P6exn")
        #PRP18
        cmd.do("color grey50, chain a and " + object_name)
        cmd.do("color grey50, PRP18_P6exn")
        #SLU7
        cmd.do("color grey50, chain c and " + object_name)
        cmd.do("color grey50, SLU7_P6exn")
        #unassigned
        cmd.do("color grey50, chain X and " + object_name)
        cmd.do("color grey50, unassigned_P6exn")
    
    for n in cmd.get_names("all"):
        if '5ylz' in n:
            object_name = n
            print(object_name)
            
    if '5ylz' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_P5ylz")
        #BUD31
        cmd.do("color dirtyviolet, chain L and " + object_name)
        cmd.do("color dirtyviolet, BUD31_P5ylz")
        #CEF1
        cmd.do("color raspberry, chain J and " + object_name)
        cmd.do("color raspberry, CEF1_P5ylz")
        #CLF1
        cmd.do("color raspberry, chain I and " + object_name)
        cmd.do("color raspberry, CLF1_P5ylz")
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_P5ylz")
        #CWC2_hRBM22
        cmd.do("color ruby, chain N and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_P5ylz")
        #CWC21
        cmd.do("color violetpurple, chain R and " + object_name)
        cmd.do("color violetpurple, CWC21_P5ylz")
        #CWC22
        cmd.do("color bluewhite, chain S and " + object_name)
        cmd.do("color bluewhite, CWC22_P5ylz")
        #LEA1
        cmd.do("color palegreen, chain o and " + object_name)
        cmd.do("color palegreen, LEA1_P5ylz")
        #Msl1
        cmd.do("color palegreen, chain p and " + object_name)
        cmd.do("color palegreen, Msl1_P5ylz")
        #PRP45
        cmd.do("color lightpink, chain Q and " + object_name)
        cmd.do("color lightpink, PRP45_P5ylz")
        #CDC40
        cmd.do("color dirtyviolet, chain T and " + object_name)
        cmd.do("color dirtyviolet, CDC40_P5ylz")
        cmd.do("color grey70, chain q and 5ylz")
        cmd.do("color grey70, chain r and 5ylz")
        cmd.do("color grey70, chain s and 5ylz")
        cmd.do("color grey70, chain t and 5ylz")
        #PRP46
        cmd.do("color lightblue, chain O and " + object_name)
        cmd.do("color lightblue, PRP46_P5ylz")
        #SLT11/ECM2
        cmd.do("color chocolate, chain M and " + object_name)
        cmd.do("color chocolate, SLT11/ECM2_P5ylz")
        #SNT309
        cmd.do("color grey70, chain G and " + object_name)
        cmd.do("color grey70, SNT309_P5ylz")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_P5ylz")
        #SYF2
        cmd.do("color brightorange, chain K and " + object_name)
        cmd.do("color brightorange, SYF2_P5ylz")
        #SYF1
        cmd.do("color brightorange, chain H and " + object_name)
        cmd.do("color brightorange, SYF1_P5ylz")
        #U2
        cmd.do("color forest, chain F and " + object_name)
        cmd.do("color forest, U2_P5ylz")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_P5ylz")
        #U6
        cmd.do("color firebrick, chain D and " + object_name)
        cmd.do("color firebrick, U6_P5ylz")
        #Intron
        cmd.do("color black, chain E and " + object_name)
        cmd.do("color black, Intron_P5ylz")
        cmd.do("color grey50, chain a and 5ylz")
        cmd.do("color grey50, chain h and 5ylz")
        cmd.do("color grey50, chain b and 5ylz")
        cmd.do("color grey50, chain i and 5ylz")
        cmd.do("color grey50, chain c and 5ylz")
        cmd.do("color grey50, chain j and 5ylz")
        cmd.do("color grey50, chain d and 5ylz")
        cmd.do("color grey50, chain k and 5ylz")
        cmd.do("color grey50, chain e and 5ylz")
        cmd.do("color grey50, chain l and 5ylz")
        cmd.do("color grey50, chain f and 5ylz")
        cmd.do("color grey50, chain m and 5ylz")
        cmd.do("color grey50, chain g and 5ylz")
        cmd.do("color grey50, chain n and 5ylz")
        #PRP22
        cmd.do("color grey50, chain W and " + object_name)
        cmd.do("color grey50, PRP22_P5ylz")
        #PRP18
        cmd.do("color grey50, chain U and " + object_name)
        cmd.do("color grey50, PRP18_P5ylz")
        #SLU7
        cmd.do("color grey50, chain V and " + object_name)
        cmd.do("color grey50, SLU7_P5ylz")
    
    for n in cmd.get_names("all"):
        if '5y88' in n:
            object_name = n
            print(object_name)
            
    if '5y88' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_I5y88")
        #BUD31
        cmd.do("color dirtyviolet, chain L and " + object_name)
        cmd.do("color dirtyviolet, BUD31_I5y88")
        #CLF1
        cmd.do("color raspberry, chain I and " + object_name)
        cmd.do("color raspberry, CLF1_I5y88")
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_I5y88")
        #CWC16/YJU2
        cmd.do("color lightteal, chain R and " + object_name)
        cmd.do("color lightteal, CWC16/YJU2_I5y88")
        #CWC2_hRBM22
        cmd.do("color ruby, chain N and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_I5y88")
        #CWC25
        cmd.do("color deepteal, chain G and " + object_name)
        cmd.do("color deepteal, CWC25_I5y88")
        #Intron_2
        cmd.do("color black, chain E and " + object_name)
        cmd.do("color black, Intron_2_I5y88")
        #LEA1
        cmd.do("color palegreen, chain o and " + object_name)
        cmd.do("color palegreen, LEA1_I5y88")
        #Msl1
        cmd.do("color palegreen, chain p and " + object_name)
        cmd.do("color palegreen, Msl1_I5y88")
        #PRP45
        cmd.do("color lightpink, chain Q and " + object_name)
        cmd.do("color lightpink, PRP45_I5y88")
        #CDC40
        cmd.do("color dirtyviolet, chain S and " + object_name)
        cmd.do("color dirtyviolet, CDC40_I5y88")
        cmd.do("color grey70, chain q and 5y88")
        cmd.do("color grey70, chain r and 5y88")
        cmd.do("color grey70, chain s and 5y88")
        cmd.do("color grey70, chain t and 5y88")
        #PRP46
        cmd.do("color lightblue, chain O and " + object_name)
        cmd.do("color lightblue, PRP46_I5y88")
        #SLT11/ECM2
        cmd.do("color chocolate, chain M and " + object_name)
        cmd.do("color chocolate, SLT11/ECM2_I5y88")
        #SNT309
        cmd.do("color grey70, chain G and " + object_name)
        cmd.do("color grey70, SNT309_I5y88")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_I5y88")
        #SYF2
        cmd.do("color brightorange, chain K and " + object_name)
        cmd.do("color brightorange, SYF2_I5y88")
        #SYF1
        cmd.do("color brightorange, chain H and " + object_name)
        cmd.do("color brightorange, SYF1_I5y88")
        #U2
        cmd.do("color forest, chain F and " + object_name)
        cmd.do("color forest, U2_I5y88")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_I5y88")
        #U6
        cmd.do("color firebrick, chain D and " + object_name)
        cmd.do("color firebrick, U6_I5y88")
        #Intron
        cmd.do("color black, chain x and " + object_name)
        cmd.do("color black, Intron_I5y88")
        #RNA
        cmd.do("color grey50, chain x and " + object_name)
        cmd.do("color grey50, RNA_I5y88")
        #cwc23
        cmd.do("color grey50, chain T and " + object_name)
        cmd.do("color grey50, cwc23_I5y88")
        #SPP382
        cmd.do("color grey50, chain U and " + object_name)
        cmd.do("color grey50, SPP382_I5y88")
        #NTR2
        cmd.do("color grey50, chain V and " + object_name)
        cmd.do("color grey50, NTR2_I5y88")
        #PRP43
        cmd.do("color grey50, chain W and " + object_name)
        cmd.do("color grey50, PRP43_I5y88")
        cmd.do("color grey50, chain a and 5y88")
        cmd.do("color grey50, chain h and 5y88")
        cmd.do("color grey50, chain b and 5y88")
        cmd.do("color grey50, chain i and 5y88")
        cmd.do("color grey50, chain c and 5y88")
        cmd.do("color grey50, chain j and 5y88")
        cmd.do("color grey50, chain d and 5y88")
        cmd.do("color grey50, chain k and 5y88")
        cmd.do("color grey50, chain e and 5y88")
        cmd.do("color grey50, chain l and 5y88")
        cmd.do("color grey50, chain f and 5y88")
        cmd.do("color grey50, chain m and 5y88")
        cmd.do("color grey50, chain g and 5y88")
        cmd.do("color grey50, chain n and 5y88")
    
    for n in cmd.get_names("all"):
        if '3jb9' in n:
            object_name = n
            print(object_name)
            
    if '3jb9' in name.lower():
        #U2
        cmd.do("color forest, chain P and " + object_name)
        cmd.do("color forest, U2_3jb9")
        #U5
        cmd.do("color density, chain C and " + object_name)
        cmd.do("color density, U5_3jb9")
        #U6
        cmd.do("color firebrick, chain N and " + object_name)
        cmd.do("color firebrick, U6_3jb9")
        cmd.do("color black, chain O and 3jb9")
        cmd.do("color black, chain Q and 3jb9")
        #Spp42_yPrp8
        cmd.do("color grey50, chain A and " + object_name)
        cmd.do("color grey50, Spp42_yPrp8_3jb9")
        #CWF15_yCWC15
        cmd.do("color orange, chain h and " + object_name)
        cmd.do("color orange, CWF15_yCWC15_3jb9")
    
    for n in cmd.get_names("all"):
        if '6icz' in n:
            object_name = n
            print(object_name)
            
    if '6icz' in name.lower():
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_hP_6icz")
        #U2
        cmd.do("color forest, chain H and " + object_name)
        cmd.do("color forest, U2_hP_6icz")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_hP_6icz")
        #U6
        cmd.do("color firebrick, chain F and " + object_name)
        cmd.do("color firebrick, U6_hP_6icz")
        #Intron
        cmd.do("color black, chain G and " + object_name)
        cmd.do("color black, Intron_hP_6icz")
        #cwc23
        cmd.do("color grey50, chain 6ICZ and " + object_name)
        cmd.do("color grey50, cwc23_hP_6icz")
    
    for n in cmd.get_names("all"):
        if '6ff7' in n:
            object_name = n
            print(object_name)
            
    if '6ff7' in name.lower():
        #CWC15
        cmd.do("color orange, chain R and " + object_name)
        cmd.do("color orange, CWC15_hBa_6ff7")
        #CWC2_hRBM22
        cmd.do("color ruby, chain P and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_hBa_6ff7")
        #U2
        cmd.do("color forest, chain 2 and " + object_name)
        cmd.do("color forest, U2_hBa_6ff7")
        #U5
        cmd.do("color density, chain 5 and " + object_name)
        cmd.do("color density, U5_hBa_6ff7")
        #U6
        cmd.do("color firebrick, chain 6 and " + object_name)
        cmd.do("color firebrick, U6_hBa_6ff7")
        #Intron
        cmd.do("color black, chain Z and " + object_name)
        cmd.do("color black, Intron_hBa_6ff7")
    
    for n in cmd.get_names("all"):
        if '5yzg' in n:
            object_name = n
            print(object_name)
            
    if '5yzg' in name.lower():
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_hC_5yzg")
        #CWC2_hRBM22
        cmd.do("color ruby, chain O and " + object_name)
        cmd.do("color ruby, CWC2_hRBM22_hC_5yzg")
        #CWC25
        cmd.do("color deepteal, chain X and " + object_name)
        cmd.do("color deepteal, CWC25_hC_5yzg")
        #PRP16_hDHX38
        cmd.do("color smudge, chain Z and " + object_name)
        cmd.do("color smudge, PRP16_hDHX38_hC_5yzg")
        #U2
        cmd.do("color forest, chain H and " + object_name)
        cmd.do("color forest, U2_hC_5yzg")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_hC_5yzg")
        #U6
        cmd.do("color firebrick, chain F and " + object_name)
        cmd.do("color firebrick, U6_hC_5yzg")
        #Intron
        cmd.do("color black, chain G and " + object_name)
        cmd.do("color black, Intron_hC_5yzg")
    
    for n in cmd.get_names("all"):
        if '5xjc' in n:
            object_name = n
            print(object_name)
            
    if '5xjc' in name.lower():
        #CWC15
        cmd.do("color orange, chain P and " + object_name)
        cmd.do("color orange, CWC15_hX_5xjc")
        #CWC25
        cmd.do("color deepteal, chain X and " + object_name)
        cmd.do("color deepteal, CWC25_hX_5xjc")
        #U2
        cmd.do("color forest, chain H and " + object_name)
        cmd.do("color forest, U2_hX_5xjc")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_hX_5xjc")
        #U6
        cmd.do("color firebrick, chain F and " + object_name)
        cmd.do("color firebrick, U6_hX_5xjc")
        #Intron
        cmd.do("color black, chain G and " + object_name)
        cmd.do("color black, Intron_hX_5xjc")
        #PRKRIP1
        cmd.do("color grey50, chain X and " + object_name)
        cmd.do("color grey50, PRKRIP1_hX_5xjc")
    
    for n in cmd.get_names("all"):
        if '5gan' in n:
            object_name = n
            print(object_name)
            
    if '5gan' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_y3_5gan")
        #BRR2
        cmd.do("color grey60, chain B and " + object_name)
        cmd.do("color grey60, BRR2_y3_5gan")
        #PRP45
        cmd.do("color lightpink, chain H and " + object_name)
        cmd.do("color lightpink, PRP45_y3_5gan")
        #SNU114
        cmd.do("color slate, chain C and " + object_name)
        cmd.do("color slate, SNU114_y3_5gan")
        #U5
        cmd.do("color density, chain U and " + object_name)
        cmd.do("color density, U5_y3_5gan")
        #U6
        cmd.do("color firebrick, chain W and " + object_name)
        cmd.do("color firebrick, U6_y3_5gan")
        #U4
        cmd.do("color brown, chain V and " + object_name)
        cmd.do("color brown, U4_y3_5gan")
        #PRP31
        cmd.do("color grey50, chain F and " + object_name)
        cmd.do("color grey50, PRP31_y3_5gan")
        #PRP6
        cmd.do("color grey50, chain J and " + object_name)
        cmd.do("color grey50, PRP6_y3_5gan")
        #PRP3
        cmd.do("color grey50, chain G and " + object_name)
        cmd.do("color grey50, PRP3_y3_5gan")
        #DIB1
        cmd.do("color grey50, chain D and " + object_name)
        cmd.do("color grey50, DIB1_y3_5gan")
        #SNU13
        cmd.do("color grey50, chain K and " + object_name)
        cmd.do("color grey50, SNU13_y3_5gan")
        #LSM8
        cmd.do("color grey50, chain 8 and " + object_name)
        cmd.do("color grey50, LSM8_y3_5gan")
        #LSM2
        cmd.do("color grey50, chain 2 and " + object_name)
        cmd.do("color grey50, LSM2_y3_5gan")
        #LSM3
        cmd.do("color grey50, chain 3 and " + object_name)
        cmd.do("color grey50, LSM3_y3_5gan")
        #LSM6
        cmd.do("color grey50, chain 6 and " + object_name)
        cmd.do("color grey50, LSM6_y3_5gan")
        #LSM5
        cmd.do("color grey50, chain 5 and " + object_name)
        cmd.do("color grey50, LSM5_y3_5gan")
        #LSM7
        cmd.do("color grey50, chain 7 and " + object_name)
        cmd.do("color grey50, LSM7_y3_5gan")
        #LSM4
        cmd.do("color grey50, chain 4 and " + object_name)
        cmd.do("color grey50, LSM4_y3_5gan")
        #SNU66
        cmd.do("color grey50, chain E and " + object_name)
        cmd.do("color grey50, SNU66_y3_5gan")
        cmd.do("color grey50, chain b and 5gan")
        cmd.do("color grey50, chain k and 5gan")
        cmd.do("color grey50, chain e and 5gan")
        cmd.do("color grey50, chain p and 5gan")
        cmd.do("color grey50, chain f and 5gan")
        cmd.do("color grey50, chain q and 5gan")
        cmd.do("color grey50, chain g and 5gan")
        cmd.do("color grey50, chain r and 5gan")
        cmd.do("color grey50, chain d and 5gan")
        cmd.do("color grey50, chain n and 5gan")
        cmd.do("color grey50, chain h and 5gan")
        cmd.do("color grey50, chain l and 5gan")
        cmd.do("color grey50, chain j and 5gan")
        cmd.do("color grey50, chain m and 5gan")
        #unassigned
        cmd.do("color grey50, chain X and " + object_name)
        cmd.do("color grey50, unassigned_y3_5gan")
    
    for n in cmd.get_names("all"):
        if '6qw6' in n:
            object_name = n
            print(object_name)
            
    if '6qw6' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain 5A and " + object_name)
        cmd.do("color skyblue, PRP8_h3_6qw6")
        #U5
        cmd.do("color density, chain 5 and " + object_name)
        cmd.do("color density, U5_h3_6qw6")
        #U6
        cmd.do("color firebrick, chain 6 and " + object_name)
        cmd.do("color firebrick, U6_h3_6qw6")
        #U4
        cmd.do("color brown, chain 4 and " + object_name)
        cmd.do("color brown, U4_h3_6qw6")
        #Prp28
        cmd.do("color red, chain 5X and " + object_name)
        cmd.do("color red, Prp28_h3_6qw6")
    
    for n in cmd.get_names("all"):
        if '3jcr' in n:
            object_name = n
            print(object_name)
            
    if '3jcr' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_h3_3jcr")
        #BRR2
        cmd.do("color grey60, chain C and " + object_name)
        cmd.do("color grey60, BRR2_h3_3jcr")
        #SNU114
        cmd.do("color slate, chain B and " + object_name)
        cmd.do("color slate, SNU114_h3_3jcr")
        #U5
        cmd.do("color density, chain H and " + object_name)
        cmd.do("color density, U5_h3_3jcr")
        #U6
        cmd.do("color firebrick, chain N and " + object_name)
        cmd.do("color firebrick, U6_h3_3jcr")
        #U4
        cmd.do("color brown, chain M and " + object_name)
        cmd.do("color brown, U4_h3_3jcr")
        #PRP4
        cmd.do("color grey50, chain L and " + object_name)
        cmd.do("color grey50, PRP4_h3_3jcr")
        #PRP31
        cmd.do("color grey50, chain J and " + object_name)
        cmd.do("color grey50, PRP31_h3_3jcr")
        #PRP6
        cmd.do("color grey50, chain G and " + object_name)
        cmd.do("color grey50, PRP6_h3_3jcr")
        #PRP3
        cmd.do("color grey50, chain K and " + object_name)
        cmd.do("color grey50, PRP3_h3_3jcr")
        #SNU13
        cmd.do("color grey50, chain I and " + object_name)
        cmd.do("color grey50, SNU13_h3_3jcr")
        #LSM8
        cmd.do("color grey50, chain 8 and " + object_name)
        cmd.do("color grey50, LSM8_h3_3jcr")
        #LSM2
        cmd.do("color grey50, chain 2 and " + object_name)
        cmd.do("color grey50, LSM2_h3_3jcr")
        #LSM3
        cmd.do("color grey50, chain 3 and " + object_name)
        cmd.do("color grey50, LSM3_h3_3jcr")
        #LSM5
        cmd.do("color grey50, chain 5 and " + object_name)
        cmd.do("color grey50, LSM5_h3_3jcr")
        #LSM7
        cmd.do("color grey50, chain 7 and " + object_name)
        cmd.do("color grey50, LSM7_h3_3jcr")
        #LSM4
        cmd.do("color grey50, chain 4 and " + object_name)
        cmd.do("color grey50, LSM4_h3_3jcr")
        cmd.do("color grey50, chain R and 3jcr")
        cmd.do("color grey50, chain r and 3jcr")
        cmd.do("color grey50, chain P and 3jcr")
        cmd.do("color grey50, chain p and 3jcr")
        cmd.do("color grey50, chain Q and 3jcr")
        cmd.do("color grey50, chain q and 3jcr")
        #U5-40K
        cmd.do("color grey71, chain D and " + object_name)
        cmd.do("color grey71, U5-40K_h3_3jcr")
        #Dim1
        cmd.do("color grey60, chain E and " + object_name)
        cmd.do("color grey60, Dim1_h3_3jcr")
        #Prp28
        cmd.do("color red, chain F and " + object_name)
        cmd.do("color red, Prp28_h3_3jcr")
        cmd.do("color grey60, chain S and 3jcr")
        cmd.do("color grey60, chain s and 3jcr")
        cmd.do("color grey60, chain T and 3jcr")
        cmd.do("color grey60, chain t and 3jcr")
        cmd.do("color grey60, chain U and 3jcr")
        cmd.do("color grey60, chain u and 3jcr")
        #Sad1
        cmd.do("color grey60, chain V and " + object_name)
        cmd.do("color grey60, Sad1_h3_3jcr")
        cmd.do("color grey60, chain O and 3jcr")
        cmd.do("color grey60, chain o and 3jcr")
    
    for n in cmd.get_names("all"):
        if '6qx9' in n:
            object_name = n
            print(object_name)
            
    if '6qx9' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain 5A and " + object_name)
        cmd.do("color skyblue, PRP8_h3_6qx9")
        #U2
        cmd.do("color forest, chain 2 and " + object_name)
        cmd.do("color forest, U2_h3_6qx9")
        #U5
        cmd.do("color density, chain 5 and " + object_name)
        cmd.do("color density, U5_h3_6qx9")
        #U6
        cmd.do("color firebrick, chain 6 and " + object_name)
        cmd.do("color firebrick, U6_h3_6qx9")
        #U4
        cmd.do("color brown, chain 4 and " + object_name)
        cmd.do("color brown, U4_h3_6qx9")
        #Intron
        cmd.do("color black, chain I and " + object_name)
        cmd.do("color black, Intron_h3_6qx9")
        #U1
        cmd.do("color green, chain 1 and " + object_name)
        cmd.do("color green, U1_h3_6qx9")
        #Prp28
        cmd.do("color red, chain 5X and " + object_name)
        cmd.do("color red, Prp28_h3_6qx9")
    
    for n in cmd.get_names("all"):
        if '6ah0' in n:
            object_name = n
            print(object_name)
            
    if '6ah0' in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and " + object_name)
        cmd.do("color skyblue, PRP8_h_Bp_6ah0")
        #U2
        cmd.do("color forest, chain H and " + object_name)
        cmd.do("color forest, U2_h_Bp_6ah0")
        #U5
        cmd.do("color density, chain B and " + object_name)
        cmd.do("color density, U5_h_Bp_6ah0")
        #U6
        cmd.do("color firebrick, chain F and " + object_name)
        cmd.do("color firebrick, U6_h_Bp_6ah0")
        #U4
        cmd.do("color brown, chain I and " + object_name)
        cmd.do("color brown, U4_h_Bp_6ah0")
        #Intron
        cmd.do("color black, chain G and " + object_name)
        cmd.do("color black, Intron_h_Bp_6ah0")
        #Prp28
        cmd.do("color red, chain X and " + object_name)
        cmd.do("color red, Prp28_h_Bp_6ah0")
