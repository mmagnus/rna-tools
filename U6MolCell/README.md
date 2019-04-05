Rearrangements within the U6 snRNA core during the transition between the two catalytic steps of splicing
Eysmont K., Matylla-Kulińska K., Jaskulska A., Magnus M., Konarska M.

Alignments
-------------------------------------------------------------------------------

- U6_RF00026.stockholm.sto
- Intron_gpII_RF00029.stockholm.sto
- U6atac_minor_RF00619.stockholm.sto

R2R
-------------------------------------------------------------------------------

    ../../src/r2r --GSC-weighted-consensus u6-only.stk u6-only-cons.sto 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1

    R2R-1.0.5
    [dhcp177-lan203] u6$ ../../src/r2r u6-only-cons.sto u6_cons.pdf
    PROCESSING: u6-only-cons
    NOTE: no #=GC R2R_LABEL line.  But that's okay.
    <<<.<<<<.....>>>..>>>>
    ssContextList:
    [0,3;19,22) T,F,T Pair
    [3,4;19,19) F,F,T InternalLoop
    [4,5;18,19) F,F,T Pair
    [5,5;16,18) F,F,T InternalLoop
    [5,8;13,16) F,F,T Pair
    [8,13;13,13) F,T,T TerminalLoop
    ssContext
            [41,43;62,64] {raw [0,3;19,22) }  Pair  openHairpin
        link
            posFrom: 2 , [41,43;62,64] {raw [0,3;19,22) }
            posTo  : 3 , [44,44;62,=] {raw [3,4;19,19) }
            default rule
            involvesCircularLayout==true
    ssContext
            [44,44;62,=] {raw [3,4;19,19) }  InternalLoop
        link
            posFrom: 3 , [44,44;62,=] {raw [3,4;19,19) }
            posTo  : 2 , [41,43;62,64] {raw [0,3;19,22) }
            default rule
            involvesCircularLayout==true
        link
            posFrom: 3 , [44,44;62,=] {raw [3,4;19,19) }
            posTo  : 4 , [46,46;61,61] {raw [4,5;18,19) }
            default rule
            involvesCircularLayout==true
    ssContext
            [46,46;61,61] {raw [4,5;18,19) }  Pair
        link
            posFrom: 4 , [46,46;61,61] {raw [4,5;18,19) }
            posTo  : 3 , [44,44;62,=] {raw [3,4;19,19) }
            default rule
            involvesCircularLayout==true
        link
            posFrom: 4 , [46,46;61,61] {raw [4,5;18,19) }
            posTo  : 17 , [47,=;59,60] {raw [5,5;16,18) }
            default rule
            involvesCircularLayout==true
    ssContext
            [47,=;59,60] {raw [5,5;16,18) }  InternalLoop
        link
            posFrom: 17 , [47,=;59,60] {raw [5,5;16,18) }
            posTo  : 4 , [46,46;61,61] {raw [4,5;18,19) }
            default rule
            involvesCircularLayout==true
        link
            posFrom: 17 , [47,=;59,60] {raw [5,5;16,18) }
            posTo  : 5 , [47,49;55,58] {raw [5,8;13,16) }
            default rule
            involvesCircularLayout==true
    ssContext
            [47,49;55,58] {raw [5,8;13,16) }  Pair
        link
            posFrom: 5 , [47,49;55,58] {raw [5,8;13,16) }
            posTo  : 17 , [47,=;59,60] {raw [5,5;16,18) }
            default rule
            involvesCircularLayout==true
        link
            posFrom: 7 , [47,49;55,58] {raw [5,8;13,16) }
            posTo  : 8 , [50,54;55,=] {raw [8,13;13,13) }
            default rule
            involvesCircularLayout==true
    ssContext
            [50,54;55,=] {raw [8,13;13,13) }  TerminalLoop
        link
            posFrom: 8 , [50,54;55,=] {raw [8,13;13,13) }
            posTo  : 7 , [47,49;55,58] {raw [5,8;13,16) }
            default rule
            involvesCircularLayout==true
    going through place_explicit's in queue
    got place_explicit from queue: calling PositionBackboneElement
        link
            posFrom: dummy ssContext
            posTo  : 0 , [41,43;62,64] {raw [0,3;19,22) }
            priorityClass=FirstRule
            default rule
            dummy rule to position first element
                at: (0,0) -90  no flip
    got place_explicit from queue: calling PositionBackboneElement
        link
            posFrom: 2 , [41,43;62,64] {raw [0,3;19,22) }
            posTo  : 3 , [44,44;62,=] {raw [3,4;19,19) }
            default rule
            involvesCircularLayout==true
    got place_explicit from queue: calling PositionBackboneElement
        link
            posFrom: 3 , [44,44;62,=] {raw [3,4;19,19) }
            posTo  : 4 , [46,46;61,61] {raw [4,5;18,19) }
            default rule
            involvesCircularLayout==true
    got place_explicit from queue: calling PositionBackboneElement
        link
            posFrom: 4 , [46,46;61,61] {raw [4,5;18,19) }
            posTo  : 17 , [47,=;59,60] {raw [5,5;16,18) }
            default rule
            involvesCircularLayout==true
    got place_explicit from queue: calling PositionBackboneElement
        link
            posFrom: 17 , [47,=;59,60] {raw [5,5;16,18) }
            posTo  : 5 , [47,49;55,58] {raw [5,8;13,16) }
            default rule
            involvesCircularLayout==true
    got place_explicit from queue: calling PositionBackboneElement
        link
            posFrom: 7 , [47,49;55,58] {raw [5,8;13,16) }
            posTo  : 8 , [50,54;55,=] {raw [8,13;13,13) }
            default rule
            involvesCircularLayout==true
        done: going through place_explicit's in queue

    doing place_defer (bulges)
    DONE PARSING: u6-only-cons

