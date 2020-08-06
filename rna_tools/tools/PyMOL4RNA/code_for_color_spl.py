
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

def spl_color():
 for name in cmd.get_names("all"):
  # cmd.do('color grey50') # off gray
 
    print(" \ Extracting mode for %s" % name)
    if 'h3_6QW6'.lower() in name.lower():
        #PRP8
        cmd.do("color skyblue, chain 5A and h3_6QW6")
        cmd.do("color skyblue, PRP8_h3_6QW6")
        #BRR2
        cmd.do("color lightsteelblue, chain 5B and h3_6QW6")
        cmd.do("color lightsteelblue, BRR2_h3_6QW6")
        #SNU114
        cmd.do("color cornflowerblue, chain 5C and h3_6QW6")
        cmd.do("color cornflowerblue, SNU114_h3_6QW6")
        #U5
        cmd.do("color density, chain 5 and h3_6QW6")
        cmd.do("color density, U5_h3_6QW6")
        #U6
        cmd.do("color firebrick, chain 6 and h3_6QW6")
        cmd.do("color firebrick, U6_h3_6QW6")
        #U4
        cmd.do("color brown, chain 4 and h3_6QW6")
        cmd.do("color brown, U4_h3_6QW6")
        #PRP4
        cmd.do("color lightorange, chain 4B and h3_6QW6")
        cmd.do("color lightorange, PRP4_h3_6QW6")
        #PRP31
        cmd.do("color yelloworange, chain 4C and h3_6QW6")
        cmd.do("color yelloworange, PRP31_h3_6QW6")
        #PRP6
        cmd.do("color palecyan, chain 5J and h3_6QW6")
        cmd.do("color palecyan, PRP6_h3_6QW6")
        #PRP3
        cmd.do("color tv_yellow, chain 4A and h3_6QW6")
        cmd.do("color tv_yellow, PRP3_h3_6QW6")
        #SNU13
        cmd.do("color wheat, chain 4D and h3_6QW6")
        cmd.do("color wheat, SNU13_h3_6QW6")
        #SNU66
        cmd.do("color violetpurple, chain S and h3_6QW6")
        cmd.do("color violetpurple, SNU66_h3_6QW6")
        #U5-40K
        cmd.do("color lightsteelblue, chain 5O and h3_6QW6")
        cmd.do("color lightsteelblue, U5-40K_h3_6QW6")
        #DIM1
        cmd.do("color orange, chain 5D and h3_6QW6")
        cmd.do("color orange, DIM1_h3_6QW6")
        #PRP28
        cmd.do("color raspberry, chain 5X and h3_6QW6")
        cmd.do("color raspberry, PRP28_h3_6QW6")
        #SAD1
        cmd.do("color lightorange, chain U and h3_6QW6")
        cmd.do("color lightorange, SAD1_h3_6QW6")
        #RBM42
        cmd.do("color lemonchiffon, chain R and h3_6QW6")
        cmd.do("color lemonchiffon, RBM42_h3_6QW6")
        #SNRNP27_27K
        cmd.do("color violetBlue, chain X and h3_6QW6")
        cmd.do("color violetBlue, SNRNP27_27K_h3_6QW6")
        cmd.do("color bluewhite, chain 51 and h3_6QW6")
        cmd.do("color bluewhite, chain 52 and h3_6QW6")
        cmd.do("color bluewhite, chain 53 and h3_6QW6")
        cmd.do("color bluewhite, chain 5b and h3_6QW6")
        cmd.do("color bluewhite, chain 5e and h3_6QW6")
        cmd.do("color bluewhite, chain 5f and h3_6QW6")
        cmd.do("color bluewhite, chain 5g and h3_6QW6")
        cmd.do("color wheat, chain 41 and h3_6QW6")
        cmd.do("color wheat, chain 42 and h3_6QW6")
        cmd.do("color wheat, chain 43 and h3_6QW6")
        cmd.do("color wheat, chain 4b and h3_6QW6")
        cmd.do("color wheat, chain 4e and h3_6QW6")
        cmd.do("color wheat, chain 4f and h3_6QW6")
        cmd.do("color wheat, chain 4g and h3_6QW6")
        cmd.do("color salmon, chain 62 and h3_6QW6")
        cmd.do("color salmon, chain 63 and h3_6QW6")
        cmd.do("color salmon, chain 64 and h3_6QW6")
        cmd.do("color salmon, chain 65 and h3_6QW6")
        cmd.do("color salmon, chain 66 and h3_6QW6")
        cmd.do("color salmon, chain 67 and h3_6QW6")
        cmd.do("color salmon, chain 68 and h3_6QW6")
    if 'hBpre_6QX9'.lower() in name.lower():
        #PRP8
        cmd.do("color skyblue, chain 5A and hBpre_6QX9")
        cmd.do("color skyblue, PRP8_hBpre_6QX9")
        #SNU114
        cmd.do("color cornflowerblue, chain 5C and hBpre_6QX9")
        cmd.do("color cornflowerblue, SNU114_hBpre_6QX9")
        #U2
        cmd.do("color forest, chain 2 and hBpre_6QX9")
        cmd.do("color forest, U2_hBpre_6QX9")
        #U5
        cmd.do("color density, chain 5 and hBpre_6QX9")
        cmd.do("color density, U5_hBpre_6QX9")
        #U6
        cmd.do("color firebrick, chain 6 and hBpre_6QX9")
        cmd.do("color firebrick, U6_hBpre_6QX9")
        #U4
        cmd.do("color brown, chain 4 and hBpre_6QX9")
        cmd.do("color brown, U4_hBpre_6QX9")
        #Intron
        cmd.do("color grey40, chain I and hBpre_6QX9")
        cmd.do("color grey40, Intron_hBpre_6QX9")
        #U1
        cmd.do("color green, chain 1 and hBpre_6QX9")
        cmd.do("color green, U1_hBpre_6QX9")
        #PRP4
        cmd.do("color lightorange, chain 4B and hBpre_6QX9")
        cmd.do("color lightorange, PRP4_hBpre_6QX9")
        #PRP31
        cmd.do("color yelloworange, chain 4C and hBpre_6QX9")
        cmd.do("color yelloworange, PRP31_hBpre_6QX9")
        #PRP6
        cmd.do("color palecyan, chain 5J and hBpre_6QX9")
        cmd.do("color palecyan, PRP6_hBpre_6QX9")
        #PRP3
        cmd.do("color tv_yellow, chain 4A and hBpre_6QX9")
        cmd.do("color tv_yellow, PRP3_hBpre_6QX9")
        #SNU66
        cmd.do("color violetpurple, chain S and hBpre_6QX9")
        cmd.do("color violetpurple, SNU66_hBpre_6QX9")
        #U5-40K
        cmd.do("color lightsteelblue, chain 5O and hBpre_6QX9")
        cmd.do("color lightsteelblue, U5-40K_hBpre_6QX9")
        #DIM1
        cmd.do("color orange, chain 5D and hBpre_6QX9")
        cmd.do("color orange, DIM1_hBpre_6QX9")
        #PRP28
        cmd.do("color raspberry, chain 5X and hBpre_6QX9")
        cmd.do("color raspberry, PRP28_hBpre_6QX9")
        #SAD1
        cmd.do("color lightorange, chain U and hBpre_6QX9")
        cmd.do("color lightorange, SAD1_hBpre_6QX9")
        #RBM42
        cmd.do("color lemonchiffon, chain R and hBpre_6QX9")
        cmd.do("color lemonchiffon, RBM42_hBpre_6QX9")
        #SNRNP27_27K
        cmd.do("color violetBlue, chain X and hBpre_6QX9")
        cmd.do("color violetBlue, SNRNP27_27K_hBpre_6QX9")
        cmd.do("color bluewhite, chain 51 and hBpre_6QX9")
        cmd.do("color bluewhite, chain 52 and hBpre_6QX9")
        cmd.do("color bluewhite, chain 53 and hBpre_6QX9")
        cmd.do("color bluewhite, chain 5b and hBpre_6QX9")
        cmd.do("color bluewhite, chain 5e and hBpre_6QX9")
        cmd.do("color bluewhite, chain 5f and hBpre_6QX9")
        cmd.do("color bluewhite, chain 5g and hBpre_6QX9")
        cmd.do("color wheat, chain 41 and hBpre_6QX9")
        cmd.do("color wheat, chain 42 and hBpre_6QX9")
        cmd.do("color wheat, chain 43 and hBpre_6QX9")
        cmd.do("color wheat, chain 4b and hBpre_6QX9")
        cmd.do("color wheat, chain 4e and hBpre_6QX9")
        cmd.do("color wheat, chain 4f and hBpre_6QX9")
        cmd.do("color wheat, chain 4g and hBpre_6QX9")
        cmd.do("color salmon, chain 62 and hBpre_6QX9")
        cmd.do("color salmon, chain 63 and hBpre_6QX9")
        cmd.do("color salmon, chain 64 and hBpre_6QX9")
        cmd.do("color salmon, chain 65 and hBpre_6QX9")
        cmd.do("color salmon, chain 66 and hBpre_6QX9")
        cmd.do("color salmon, chain 67 and hBpre_6QX9")
        cmd.do("color salmon, chain 68 and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 11 and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 12 and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 13 and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 1b and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 1e and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 1f and hBpre_6QX9")
        cmd.do("color lavenderblush, chain 1g and hBpre_6QX9")
        #U1A
        cmd.do("color lavenderblush, chain 1A and hBpre_6QX9")
        cmd.do("color lavenderblush, U1A_hBpre_6QX9")
        #U170k
        cmd.do("color thistle, chain 1K and hBpre_6QX9")
        cmd.do("color thistle, U170k_hBpre_6QX9")
        #U1C
        cmd.do("color warmpink, chain 1C and hBpre_6QX9")
        cmd.do("color warmpink, U1C_hBpre_6QX9")
        #SF3B1
        cmd.do("color lightgreen, chain B1 and hBpre_6QX9")
        cmd.do("color lightgreen, SF3B1_hBpre_6QX9")
        #SF3B3
        cmd.do("color forest, chain B3 and hBpre_6QX9")
        cmd.do("color forest, SF3B3_hBpre_6QX9")
        #SF3B2
        cmd.do("color darkgreen, chain B2 and hBpre_6QX9")
        cmd.do("color darkgreen, SF3B2_hBpre_6QX9")
        #SF3B4
        cmd.do("color limon, chain B4 and hBpre_6QX9")
        cmd.do("color limon, SF3B4_hBpre_6QX9")
        #PHF5A
        cmd.do("color limegreen, chain BP and hBpre_6QX9")
        cmd.do("color limegreen, PHF5A_hBpre_6QX9")
        #SF3A3
        cmd.do("color lightteal, chain A3 and hBpre_6QX9")
        cmd.do("color lightteal, SF3A3_hBpre_6QX9")
        #SF3A2
        cmd.do("color bluewhite, chain A2 and hBpre_6QX9")
        cmd.do("color bluewhite, SF3A2_hBpre_6QX9")
        #SF3A1
        cmd.do("color greencyan, chain A1 and hBpre_6QX9")
        cmd.do("color greencyan, SF3A1_hBpre_6QX9")
        #U2A
        cmd.do("color lightgreen, chain 2A and hBpre_6QX9")
        cmd.do("color lightgreen, U2A_hBpre_6QX9")
        #U2B
        cmd.do("color lightgreen, chain 2B and hBpre_6QX9")
        cmd.do("color lightgreen, U2B_hBpre_6QX9")
        #SF3B5
        cmd.do("color tv_green, chain B5 and hBpre_6QX9")
        cmd.do("color tv_green, SF3B5_hBpre_6QX9")
        cmd.do("color palegreen, chain 21 and hBpre_6QX9")
        cmd.do("color palegreen, chain 22 and hBpre_6QX9")
        cmd.do("color palegreen, chain 23 and hBpre_6QX9")
        cmd.do("color palegreen, chain 2b and hBpre_6QX9")
        cmd.do("color palegreen, chain 2e and hBpre_6QX9")
        cmd.do("color palegreen, chain 2f and hBpre_6QX9")
        cmd.do("color palegreen, chain 2g and hBpre_6QX9")
        #Prp4Kinase
        cmd.do("color mediumpurple, chain K and hBpre_6QX9")
        cmd.do("color mediumpurple, Prp4Kinase_hBpre_6QX9")
    if 'hB_6AHD'.lower() in name.lower():
        #PRP8
        cmd.do("color skyblue, chain A and hB_6AHD")
        cmd.do("color skyblue, PRP8_hB_6AHD")
        #BRR2
        cmd.do("color lightsteelblue, chain D and hB_6AHD")
        cmd.do("color lightsteelblue, BRR2_hB_6AHD")
        #SNU114
        cmd.do("color cornflowerblue, chain C and hB_6AHD")
        cmd.do("color cornflowerblue, SNU114_hB_6AHD")
        #U2
        cmd.do("color forest, chain H and hB_6AHD")
        cmd.do("color forest, U2_hB_6AHD")
        #U5
        cmd.do("color density, chain B and hB_6AHD")
        cmd.do("color density, U5_hB_6AHD")
        #U6
        cmd.do("color firebrick, chain F and hB_6AHD")
        cmd.do("color firebrick, U6_hB_6AHD")
        #U4
        cmd.do("color brown, chain I and hB_6AHD")
        cmd.do("color brown, U4_hB_6AHD")
        #Intron
        cmd.do("color grey40, chain G and hB_6AHD")
        cmd.do("color grey40, Intron_hB_6AHD")
        #Exon
        cmd.do("color yellow, chain G and resi \-30-\-1 and hB_6AHD")
        cmd.do("color yellow, Exon_hB_6AHD")
        #BP
        cmd.do("color purple, chain G and resi 144 and hB_6AHD")
        cmd.do("color purple, BP_hB_6AHD")
        #PRP4
        cmd.do("color lightorange, chain K and hB_6AHD")
        cmd.do("color lightorange, PRP4_hB_6AHD")
        #PRP31
        cmd.do("color yelloworange, chain L and hB_6AHD")
        cmd.do("color yelloworange, PRP31_hB_6AHD")
        #PRP6
        cmd.do("color palecyan, chain N and hB_6AHD")
        cmd.do("color palecyan, PRP6_hB_6AHD")
        #PRP3
        cmd.do("color tv_yellow, chain J and hB_6AHD")
        cmd.do("color tv_yellow, PRP3_hB_6AHD")
        #DIB1
        cmd.do("color brightorange, chain O and hB_6AHD")
        cmd.do("color brightorange, DIB1_hB_6AHD")
        #SNU13
        cmd.do("color wheat, chain M and hB_6AHD")
        cmd.do("color wheat, SNU13_hB_6AHD")
        #SNU66
        cmd.do("color violetpurple, chain 9 and hB_6AHD")
        cmd.do("color violetpurple, SNU66_hB_6AHD")
        #U5-40K
        cmd.do("color lightsteelblue, chain E and hB_6AHD")
        cmd.do("color lightsteelblue, U5-40K_hB_6AHD")
        cmd.do("color bluewhite, chain a and hB_6AHD")
        cmd.do("color bluewhite, chain b and hB_6AHD")
        cmd.do("color bluewhite, chain c and hB_6AHD")
        cmd.do("color bluewhite, chain d and hB_6AHD")
        cmd.do("color bluewhite, chain e and hB_6AHD")
        cmd.do("color bluewhite, chain f and hB_6AHD")
        cmd.do("color bluewhite, chain g and hB_6AHD")
        cmd.do("color wheat, chain U and hB_6AHD")
        cmd.do("color wheat, chain V and hB_6AHD")
        cmd.do("color wheat, chain P and hB_6AHD")
        cmd.do("color wheat, chain Q and hB_6AHD")
        cmd.do("color wheat, chain R and hB_6AHD")
        cmd.do("color wheat, chain S and hB_6AHD")
        cmd.do("color wheat, chain T and hB_6AHD")
        cmd.do("color salmon, chain q and hB_6AHD")
        cmd.do("color salmon, chain r and hB_6AHD")
        cmd.do("color salmon, chain s and hB_6AHD")
        cmd.do("color salmon, chain t and hB_6AHD")
        cmd.do("color salmon, chain x and hB_6AHD")
        cmd.do("color salmon, chain y and hB_6AHD")
        cmd.do("color salmon, chain z and hB_6AHD")
        #SF3B1
        cmd.do("color lightgreen, chain 1 and hB_6AHD")
        cmd.do("color lightgreen, SF3B1_hB_6AHD")
        #SF3B3
        cmd.do("color forest, chain 3 and hB_6AHD")
        cmd.do("color forest, SF3B3_hB_6AHD")
        #SF3B2
        cmd.do("color darkgreen, chain 2 and hB_6AHD")
        cmd.do("color darkgreen, SF3B2_hB_6AHD")
        #SF3B4
        cmd.do("color limon, chain 4 and hB_6AHD")
        cmd.do("color limon, SF3B4_hB_6AHD")
        #PHF5A
        cmd.do("color limegreen, chain 6 and hB_6AHD")
        cmd.do("color limegreen, PHF5A_hB_6AHD")
        #SF3A3
        cmd.do("color lightteal, chain w and hB_6AHD")
        cmd.do("color lightteal, SF3A3_hB_6AHD")
        #SF3A2
        cmd.do("color bluewhite, chain v and hB_6AHD")
        cmd.do("color bluewhite, SF3A2_hB_6AHD")
        #SF3A1
        cmd.do("color greencyan, chain u and hB_6AHD")
        cmd.do("color greencyan, SF3A1_hB_6AHD")
        #U2A
        cmd.do("color lightgreen, chain o and hB_6AHD")
        cmd.do("color lightgreen, U2A_hB_6AHD")
        #U2B
        cmd.do("color lightgreen, chain p and hB_6AHD")
        cmd.do("color lightgreen, U2B_hB_6AHD")
        #SF3B5
        cmd.do("color tv_green, chain 7 and hB_6AHD")
        cmd.do("color tv_green, SF3B5_hB_6AHD")
        cmd.do("color palegreen, chain i and hB_6AHD")
        cmd.do("color palegreen, chain j and hB_6AHD")
        cmd.do("color palegreen, chain k and hB_6AHD")
        cmd.do("color palegreen, chain l and hB_6AHD")
        cmd.do("color palegreen, chain m and hB_6AHD")
        cmd.do("color palegreen, chain n and hB_6AHD")
        cmd.do("color palegreen, chain h and hB_6AHD")
        #FBP21
        cmd.do("color raspberry, chain X and hB_6AHD")
        cmd.do("color raspberry, FBP21_hB_6AHD")
        #PPIH
        cmd.do("color mistyrose, chain W and hB_6AHD")
        cmd.do("color mistyrose, PPIH_hB_6AHD")
        #UBL5
        cmd.do("color orchid, chain A0 and hB_6AHD")
        cmd.do("color orchid, UBL5_hB_6AHD")
        #MFAP1
        cmd.do("color lavenderblush, chain 0 and hB_6AHD")
        cmd.do("color lavenderblush, MFAP1_hB_6AHD")
        #PRPF38A
        cmd.do("color thistle, chain Z and hB_6AHD")
        cmd.do("color thistle, PRPF38A_hB_6AHD")
        #ZMAT2
        cmd.do("color violetpurple, chain 8 and hB_6AHD")
        cmd.do("color violetpurple, ZMAT2_hB_6AHD")
        #SMU1
        cmd.do("color lemonchiffon, chain Y and hB_6AHD")
        cmd.do("color lemonchiffon, SMU1_hB_6AHD")
        #SF3B6
        cmd.do("color tv_green, chain 5 and hB_6AHD")
        cmd.do("color tv_green, SF3B6_hB_6AHD")
