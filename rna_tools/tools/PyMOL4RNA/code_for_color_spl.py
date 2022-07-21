
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

def spl_color():
 for name in cmd.get_names("all"):
  # cmd.do('color grey50') # off gray
 
    print(" \ Extracting mode for %s" % name)
    if 'hpCs-II_7W5A'.lower() in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and hpCs-II_7W5A")
        cmd.do("color skyblue, PRP8_hpCs-II_7W5A")
        #BUD31
        cmd.do("color dirtyviolet, chain N and hpCs-II_7W5A")
        cmd.do("color dirtyviolet, BUD31_hpCs-II_7W5A")
        #CWC15
        cmd.do("color orange, chain P and hpCs-II_7W5A")
        cmd.do("color orange, CWC15_hpCs-II_7W5A")
        #CWC2_RBM22
        cmd.do("color ruby, chain O and hpCs-II_7W5A")
        cmd.do("color ruby, CWC2_RBM22_hpCs-II_7W5A")
        #CWC22
        cmd.do("color bluewhite, chain V and hpCs-II_7W5A")
        cmd.do("color bluewhite, CWC22_hpCs-II_7W5A")
        #PRP46
        cmd.do("color lightblue, chain T and hpCs-II_7W5A")
        cmd.do("color lightblue, PRP46_hpCs-II_7W5A")
        #SYF2
        cmd.do("color brightorange, chain M and hpCs-II_7W5A")
        cmd.do("color brightorange, SYF2_hpCs-II_7W5A")
        #SYF1
        cmd.do("color brightorange, chain I and hpCs-II_7W5A")
        cmd.do("color brightorange, SYF1_hpCs-II_7W5A")
        #U2
        cmd.do("color forest, chain H and hpCs-II_7W5A")
        cmd.do("color forest, U2_hpCs-II_7W5A")
        #U5
        cmd.do("color density, chain B and hpCs-II_7W5A")
        cmd.do("color density, U5_hpCs-II_7W5A")
        #U6
        cmd.do("color firebrick, chain F and hpCs-II_7W5A")
        cmd.do("color firebrick, U6_hpCs-II_7W5A")
        #5EXON
        cmd.do("color yellow, chain 4 and hpCs-II_7W5A")
        cmd.do("color yellow, 5EXON_hpCs-II_7W5A")
        #Intron
        cmd.do("color grey40, chain G and hpCs-II_7W5A")
        cmd.do("color grey40, Intron_hpCs-II_7W5A")
        #SLU7
        cmd.do("color grey50, chain 1 and hpCs-II_7W5A")
        cmd.do("color grey50, SLU7_hpCs-II_7W5A")
        #PRKRIP1
        cmd.do("color grey50, chain 2 and hpCs-II_7W5A")
        cmd.do("color grey50, PRKRIP1_hpCs-II_7W5A")
        #AQR
        cmd.do("color raspberry, chain Q and hpCs-II_7W5A")
        cmd.do("color raspberry, AQR_hpCs-II_7W5A")
        #SNRP116
        cmd.do("color gray50, chain C and hpCs-II_7W5A")
        cmd.do("color gray50, SNRP116_hpCs-II_7W5A")
        #PPIL1
        cmd.do("color firebrick, chain S and hpCs-II_7W5A")
        cmd.do("color firebrick, PPIL1_hpCs-II_7W5A")
        #SNRNP200
        cmd.do("color gray51, chain D and hpCs-II_7W5A")
        cmd.do("color gray51, SNRNP200_hpCs-II_7W5A")
        #PPIE
        cmd.do("color firebrick, chain y and hpCs-II_7W5A")
        cmd.do("color firebrick, PPIE_hpCs-II_7W5A")
        #SNRNP40
        cmd.do("color gray52, chain E and hpCs-II_7W5A")
        cmd.do("color gray52, SNRNP40_hpCs-II_7W5A")
        cmd.do("color gray53, chain a and hpCs-II_7W5A")
        cmd.do("color gray53, chain h and hpCs-II_7W5A")
        cmd.do("color gray54, chain b and hpCs-II_7W5A")
        cmd.do("color gray54, chain i and hpCs-II_7W5A")
        #EIF4A3
        cmd.do("color indianred, chain u and hpCs-II_7W5A")
        cmd.do("color indianred, EIF4A3_hpCs-II_7W5A")
        cmd.do("color gray55, chain c and hpCs-II_7W5A")
        cmd.do("color gray55, chain j and hpCs-II_7W5A")
        #MAGOH
        cmd.do("color lavenderblush, chain v and hpCs-II_7W5A")
        cmd.do("color lavenderblush, MAGOH_hpCs-II_7W5A")
        cmd.do("color gray56, chain d and hpCs-II_7W5A")
        cmd.do("color gray56, chain k and hpCs-II_7W5A")
        #RBM8A
        cmd.do("color raspberry, chain w and hpCs-II_7W5A")
        cmd.do("color raspberry, RBM8A_hpCs-II_7W5A")
        cmd.do("color gray57, chain f and hpCs-II_7W5A")
        cmd.do("color gray57, chain m and hpCs-II_7W5A")
        #CASC3
        cmd.do("color white, chain x and hpCs-II_7W5A")
        cmd.do("color white, CASC3_hpCs-II_7W5A")
        cmd.do("color gray58, chain e and hpCs-II_7W5A")
        cmd.do("color gray58, chain l and hpCs-II_7W5A")
        cmd.do("color gray59, chain g and hpCs-II_7W5A")
        cmd.do("color gray59, chain n and hpCs-II_7W5A")
        #SNRPA1
        cmd.do("color gray60, chain o and hpCs-II_7W5A")
        cmd.do("color gray60, SNRPA1_hpCs-II_7W5A")
        #SNRPB2
        cmd.do("color gray61, chain p and hpCs-II_7W5A")
        cmd.do("color gray61, SNRPB2_hpCs-II_7W5A")
        #CRNKL1
        cmd.do("color gray62, chain J and hpCs-II_7W5A")
        cmd.do("color gray62, CRNKL1_hpCs-II_7W5A")
        #CDC5L
        cmd.do("color gray63, chain L and hpCs-II_7W5A")
        cmd.do("color gray63, CDC5L_hpCs-II_7W5A")
        cmd.do("color gray64, chain q and hpCs-II_7W5A")
        cmd.do("color gray64, chain r and hpCs-II_7W5A")
        cmd.do("color gray64, chain s and hpCs-II_7W5A")
        cmd.do("color gray64, chain t and hpCs-II_7W5A")
        #BCAS2
        cmd.do("color gray65, chain K and hpCs-II_7W5A")
        cmd.do("color gray65, BCAS2_hpCs-II_7W5A")
        #SKIP
        cmd.do("color gray66, chain R and hpCs-II_7W5A")
        cmd.do("color gray66, SKIP_hpCs-II_7W5A")
        #SRRM2
        cmd.do("color gray68, chain U and hpCs-II_7W5A")
        cmd.do("color gray68, SRRM2_hpCs-II_7W5A")
        #CDC40
        cmd.do("color gray69, chain W and hpCs-II_7W5A")
        cmd.do("color gray69, CDC40_hpCs-II_7W5A")
        #DHX8
        cmd.do("color gray50, chain Y and hpCs-II_7W5A")
        cmd.do("color gray50, DHX8_hpCs-II_7W5A")
        #FAM32A
        cmd.do("color gray50, chain z and hpCs-II_7W5A")
        cmd.do("color gray50, FAM32A_hpCs-II_7W5A")
        #NKAP
        cmd.do("color gray50, chain 3 and hpCs-II_7W5A")
        cmd.do("color gray50, NKAP_hpCs-II_7W5A")
        #Cactin
        cmd.do("color gray55, chain Z and hpCs-II_7W5A")
        cmd.do("color gray55, Cactin_hpCs-II_7W5A")
    if 'hpCs-I_7W59'.lower() in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and hpCs-I_7W59")
        cmd.do("color skyblue, PRP8_hpCs-I_7W59")
        #BUD31
        cmd.do("color dirtyviolet, chain N and hpCs-I_7W59")
        cmd.do("color dirtyviolet, BUD31_hpCs-I_7W59")
        #CWC15
        cmd.do("color orange, chain P and hpCs-I_7W59")
        cmd.do("color orange, CWC15_hpCs-I_7W59")
        #CWC2_RBM22
        cmd.do("color ruby, chain O and hpCs-I_7W59")
        cmd.do("color ruby, CWC2_RBM22_hpCs-I_7W59")
        #CWC22
        cmd.do("color bluewhite, chain V and hpCs-I_7W59")
        cmd.do("color bluewhite, CWC22_hpCs-I_7W59")
        #PRP46
        cmd.do("color lightblue, chain T and hpCs-I_7W59")
        cmd.do("color lightblue, PRP46_hpCs-I_7W59")
        #SYF2
        cmd.do("color brightorange, chain M and hpCs-I_7W59")
        cmd.do("color brightorange, SYF2_hpCs-I_7W59")
        #SYF1
        cmd.do("color brightorange, chain I and hpCs-I_7W59")
        cmd.do("color brightorange, SYF1_hpCs-I_7W59")
        #U2
        cmd.do("color forest, chain H and hpCs-I_7W59")
        cmd.do("color forest, U2_hpCs-I_7W59")
        #U5
        cmd.do("color density, chain B and hpCs-I_7W59")
        cmd.do("color density, U5_hpCs-I_7W59")
        #U6
        cmd.do("color firebrick, chain F and hpCs-I_7W59")
        cmd.do("color firebrick, U6_hpCs-I_7W59")
        #5EXON
        cmd.do("color yellow, chain 4 and hpCs-I_7W59")
        cmd.do("color yellow, 5EXON_hpCs-I_7W59")
        #Intron
        cmd.do("color grey40, chain G and hpCs-I_7W59")
        cmd.do("color grey40, Intron_hpCs-I_7W59")
        #SLU7
        cmd.do("color grey50, chain 1 and hpCs-I_7W59")
        cmd.do("color grey50, SLU7_hpCs-I_7W59")
        #PRKRIP1
        cmd.do("color grey50, chain 2 and hpCs-I_7W59")
        cmd.do("color grey50, PRKRIP1_hpCs-I_7W59")
        #AQR
        cmd.do("color raspberry, chain Q and hpCs-I_7W59")
        cmd.do("color raspberry, AQR_hpCs-I_7W59")
        #SNRP116
        cmd.do("color gray50, chain C and hpCs-I_7W59")
        cmd.do("color gray50, SNRP116_hpCs-I_7W59")
        #PPIL1
        cmd.do("color firebrick, chain S and hpCs-I_7W59")
        cmd.do("color firebrick, PPIL1_hpCs-I_7W59")
        #SNRNP200
        cmd.do("color gray51, chain D and hpCs-I_7W59")
        cmd.do("color gray51, SNRNP200_hpCs-I_7W59")
        #PPIE
        cmd.do("color firebrick, chain y and hpCs-I_7W59")
        cmd.do("color firebrick, PPIE_hpCs-I_7W59")
        #SNRNP40
        cmd.do("color gray52, chain E and hpCs-I_7W59")
        cmd.do("color gray52, SNRNP40_hpCs-I_7W59")
        cmd.do("color gray53, chain a and hpCs-I_7W59")
        cmd.do("color gray53, chain h and hpCs-I_7W59")
        cmd.do("color gray54, chain b and hpCs-I_7W59")
        cmd.do("color gray54, chain i and hpCs-I_7W59")
        #EIF4A3
        cmd.do("color indianred, chain u and hpCs-I_7W59")
        cmd.do("color indianred, EIF4A3_hpCs-I_7W59")
        cmd.do("color gray55, chain c and hpCs-I_7W59")
        cmd.do("color gray55, chain j and hpCs-I_7W59")
        #MAGOH
        cmd.do("color lavenderblush, chain v and hpCs-I_7W59")
        cmd.do("color lavenderblush, MAGOH_hpCs-I_7W59")
        cmd.do("color gray56, chain d and hpCs-I_7W59")
        cmd.do("color gray56, chain k and hpCs-I_7W59")
        #RBM8A
        cmd.do("color raspberry, chain w and hpCs-I_7W59")
        cmd.do("color raspberry, RBM8A_hpCs-I_7W59")
        cmd.do("color gray57, chain f and hpCs-I_7W59")
        cmd.do("color gray57, chain m and hpCs-I_7W59")
        #CASC3
        cmd.do("color white, chain x and hpCs-I_7W59")
        cmd.do("color white, CASC3_hpCs-I_7W59")
        cmd.do("color gray58, chain e and hpCs-I_7W59")
        cmd.do("color gray58, chain l and hpCs-I_7W59")
        cmd.do("color gray59, chain g and hpCs-I_7W59")
        cmd.do("color gray59, chain n and hpCs-I_7W59")
        #SNRPA1
        cmd.do("color gray60, chain o and hpCs-I_7W59")
        cmd.do("color gray60, SNRPA1_hpCs-I_7W59")
        #SNRPB2
        cmd.do("color gray61, chain p and hpCs-I_7W59")
        cmd.do("color gray61, SNRPB2_hpCs-I_7W59")
        #CRNKL1
        cmd.do("color gray62, chain J and hpCs-I_7W59")
        cmd.do("color gray62, CRNKL1_hpCs-I_7W59")
        #CDC5L
        cmd.do("color gray63, chain L and hpCs-I_7W59")
        cmd.do("color gray63, CDC5L_hpCs-I_7W59")
        cmd.do("color gray64, chain q and hpCs-I_7W59")
        cmd.do("color gray64, chain r and hpCs-I_7W59")
        cmd.do("color gray64, chain s and hpCs-I_7W59")
        cmd.do("color gray64, chain t and hpCs-I_7W59")
        #BCAS2
        cmd.do("color gray65, chain K and hpCs-I_7W59")
        cmd.do("color gray65, BCAS2_hpCs-I_7W59")
        #SKIP
        cmd.do("color gray66, chain R and hpCs-I_7W59")
        cmd.do("color gray66, SKIP_hpCs-I_7W59")
        #SRRM2
        cmd.do("color gray68, chain U and hpCs-I_7W59")
        cmd.do("color gray68, SRRM2_hpCs-I_7W59")
        #CDC40
        cmd.do("color gray69, chain W and hpCs-I_7W59")
        cmd.do("color gray69, CDC40_hpCs-I_7W59")
        #DHX8
        cmd.do("color gray50, chain Y and hpCs-I_7W59")
        cmd.do("color gray50, DHX8_hpCs-I_7W59")
        #FAM32A
        cmd.do("color gray50, chain z and hpCs-I_7W59")
        cmd.do("color gray50, FAM32A_hpCs-I_7W59")
        #NKAP
        cmd.do("color gray50, chain 3 and hpCs-I_7W59")
        cmd.do("color gray50, NKAP_hpCs-I_7W59")
        #Cactin
        cmd.do("color gray55, chain Z and hpCs-I_7W59")
        cmd.do("color gray55, Cactin_hpCs-I_7W59")
    if 'hpCs_7W5B'.lower() in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and hpCs_7W5B")
        cmd.do("color skyblue, PRP8_hpCs_7W5B")
        #BUD31
        cmd.do("color dirtyviolet, chain N and hpCs_7W5B")
        cmd.do("color dirtyviolet, BUD31_hpCs_7W5B")
        #CWC15
        cmd.do("color orange, chain P and hpCs_7W5B")
        cmd.do("color orange, CWC15_hpCs_7W5B")
        #CWC2_RBM22
        cmd.do("color ruby, chain O and hpCs_7W5B")
        cmd.do("color ruby, CWC2_RBM22_hpCs_7W5B")
        #CWC22
        cmd.do("color bluewhite, chain V and hpCs_7W5B")
        cmd.do("color bluewhite, CWC22_hpCs_7W5B")
        #PRP46
        cmd.do("color lightblue, chain T and hpCs_7W5B")
        cmd.do("color lightblue, PRP46_hpCs_7W5B")
        #SYF2
        cmd.do("color brightorange, chain M and hpCs_7W5B")
        cmd.do("color brightorange, SYF2_hpCs_7W5B")
        #SYF1
        cmd.do("color brightorange, chain I and hpCs_7W5B")
        cmd.do("color brightorange, SYF1_hpCs_7W5B")
        #U2
        cmd.do("color forest, chain H and hpCs_7W5B")
        cmd.do("color forest, U2_hpCs_7W5B")
        #U5
        cmd.do("color density, chain B and hpCs_7W5B")
        cmd.do("color density, U5_hpCs_7W5B")
        #U6
        cmd.do("color firebrick, chain F and hpCs_7W5B")
        cmd.do("color firebrick, U6_hpCs_7W5B")
        #5EXON
        cmd.do("color yellow, chain 4 and hpCs_7W5B")
        cmd.do("color yellow, 5EXON_hpCs_7W5B")
        #Intron
        cmd.do("color grey40, chain G and hpCs_7W5B")
        cmd.do("color grey40, Intron_hpCs_7W5B")
        #SLU7
        cmd.do("color grey50, chain 1 and hpCs_7W5B")
        cmd.do("color grey50, SLU7_hpCs_7W5B")
        #PRKRIP1
        cmd.do("color grey50, chain 2 and hpCs_7W5B")
        cmd.do("color grey50, PRKRIP1_hpCs_7W5B")
        #AQR
        cmd.do("color raspberry, chain Q and hpCs_7W5B")
        cmd.do("color raspberry, AQR_hpCs_7W5B")
        #SNRP116
        cmd.do("color gray50, chain C and hpCs_7W5B")
        cmd.do("color gray50, SNRP116_hpCs_7W5B")
        #PPIL1
        cmd.do("color firebrick, chain S and hpCs_7W5B")
        cmd.do("color firebrick, PPIL1_hpCs_7W5B")
        #SNRNP200
        cmd.do("color gray51, chain D and hpCs_7W5B")
        cmd.do("color gray51, SNRNP200_hpCs_7W5B")
        #PPIE
        cmd.do("color firebrick, chain y and hpCs_7W5B")
        cmd.do("color firebrick, PPIE_hpCs_7W5B")
        #SNRNP40
        cmd.do("color gray52, chain E and hpCs_7W5B")
        cmd.do("color gray52, SNRNP40_hpCs_7W5B")
        cmd.do("color gray53, chain a and hpCs_7W5B")
        cmd.do("color gray53, chain h and hpCs_7W5B")
        cmd.do("color gray54, chain b and hpCs_7W5B")
        cmd.do("color gray54, chain i and hpCs_7W5B")
        #EIF4A3
        cmd.do("color indianred, chain u and hpCs_7W5B")
        cmd.do("color indianred, EIF4A3_hpCs_7W5B")
        cmd.do("color gray55, chain c and hpCs_7W5B")
        cmd.do("color gray55, chain j and hpCs_7W5B")
        #MAGOH
        cmd.do("color lavenderblush, chain v and hpCs_7W5B")
        cmd.do("color lavenderblush, MAGOH_hpCs_7W5B")
        cmd.do("color gray56, chain d and hpCs_7W5B")
        cmd.do("color gray56, chain k and hpCs_7W5B")
        #RBM8A
        cmd.do("color raspberry, chain w and hpCs_7W5B")
        cmd.do("color raspberry, RBM8A_hpCs_7W5B")
        cmd.do("color gray57, chain f and hpCs_7W5B")
        cmd.do("color gray57, chain m and hpCs_7W5B")
        #CASC3
        cmd.do("color white, chain x and hpCs_7W5B")
        cmd.do("color white, CASC3_hpCs_7W5B")
        cmd.do("color gray58, chain e and hpCs_7W5B")
        cmd.do("color gray58, chain l and hpCs_7W5B")
        cmd.do("color gray59, chain g and hpCs_7W5B")
        cmd.do("color gray59, chain n and hpCs_7W5B")
        #SNRPA1
        cmd.do("color gray60, chain o and hpCs_7W5B")
        cmd.do("color gray60, SNRPA1_hpCs_7W5B")
        #SNRPB2
        cmd.do("color gray61, chain p and hpCs_7W5B")
        cmd.do("color gray61, SNRPB2_hpCs_7W5B")
        #CRNKL1
        cmd.do("color gray62, chain J and hpCs_7W5B")
        cmd.do("color gray62, CRNKL1_hpCs_7W5B")
        #CDC5L
        cmd.do("color gray63, chain L and hpCs_7W5B")
        cmd.do("color gray63, CDC5L_hpCs_7W5B")
        cmd.do("color gray64, chain q and hpCs_7W5B")
        cmd.do("color gray64, chain r and hpCs_7W5B")
        cmd.do("color gray64, chain s and hpCs_7W5B")
        cmd.do("color gray64, chain t and hpCs_7W5B")
        #BCAS2
        cmd.do("color gray65, chain K and hpCs_7W5B")
        cmd.do("color gray65, BCAS2_hpCs_7W5B")
        #SKIP
        cmd.do("color gray66, chain R and hpCs_7W5B")
        cmd.do("color gray66, SKIP_hpCs_7W5B")
        #SRRM2
        cmd.do("color gray68, chain U and hpCs_7W5B")
        cmd.do("color gray68, SRRM2_hpCs_7W5B")
        #CDC40
        cmd.do("color gray69, chain W and hpCs_7W5B")
        cmd.do("color gray69, CDC40_hpCs_7W5B")
        #DHX8
        cmd.do("color gray50, chain Y and hpCs_7W5B")
        cmd.do("color gray50, DHX8_hpCs_7W5B")
        #FAM32A
        cmd.do("color gray50, chain z and hpCs_7W5B")
        cmd.do("color gray50, FAM32A_hpCs_7W5B")
        #NKAP
        cmd.do("color gray50, chain 3 and hpCs_7W5B")
        cmd.do("color gray50, NKAP_hpCs_7W5B")
        #Cactin
        cmd.do("color gray55, chain Z and hpCs_7W5B")
        cmd.do("color gray55, Cactin_hpCs_7W5B")
