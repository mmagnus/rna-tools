"""
pucker.py is a PyMol script that returns the sugar pucker information (phase, amplitude, pucker) for a given selection.

This script uses its own dihedral calculation scheme rather than the get_dihedral command. Thus, it is lightning fast!

If a selection does not contain any ribose sugars then an error message is returned.

Author: Sean Law
Institute: University of Michigan
E-mail: seanlaw@umich.edu

Adapted for Python 3 by magnus <mmagnus@fas.harvard.edu>

Source: <https://pymolwiki.org/index.php/Pucker>
"""
from pymol.cgo import *        # get constants
from math import *
from pymol import cmd

def pucker(selection,state=1):

    #Comparison to output from 3DNA is identical
    #Note that the 3DNA output has the second chain
    #printed out in reverse order
    state=int(state)
    if state == 0:
        for state in range(1,cmd.count_states()+1):
            sel=cmd.get_model(selection,state)
            if state == 1:
                print(" ") #Blank line to separate chain output
                print("%5s  %6s      %6s  %8s  Residue" % ("State","Phase","Amp", "Pucker"))
            get_pucker(sel,iterate=state)
    else:
        sel=cmd.get_model(selection,state)
        get_pucker(sel,iterate=0)
    return

def get_pucker (sel,iterate=0):
    iterate=int(iterate)
    first=1
    old=" "
    oldchain=" "
    residue={}
    theta={}
    count=0
    for atom in sel.atom:
        new=atom.chain+" "+str(atom.resi)
        newchain=atom.chain+" "+atom.segi
        if oldchain != newchain and first:
            if iterate == 0:
                print(" ") #Blank line to separate chain output
                print("%6s  %6s      %8s  Residue" % ("Phase", "Amp", "Pucker"))
        if new != old and not first:
            #Check that all 5 atoms exist
            if len(residue) == 15:
                #Construct planes
                get_dihedrals(residue,theta)
                #Calculate pucker
                info = pseudo(residue,theta)
                print(info+"  "+old)
            else:
                print("There is no sugar in this residue")
            if oldchain != newchain:
                print(" ") #Blank line to separate chain output
                print("%6s  %6s      %8s  Residue" % ("Phase", "Amp", "Pucker"))
            #Clear values
            residue={}
            dihedrals={}
            theta={}
            #Store new value
            store_atom(atom,residue)
        else:
            store_atom(atom,residue)
        first=0
        old=new
        oldchain=newchain
    #Final Residue
    #Calculate dihedrals for final residue
    if len(residue) == 15:
        #Construct planes
        get_dihedrals(residue,theta)
        #Calculate pucker for final residue
        info = pseudo(residue,theta)
        if iterate == 0:
            print(info+"  "+old)
        else:
            print("%5d  %s" % (iterate,info+"  "+old))
    else:
        print("There is no sugar in this residue")
    return

def sele_exists(sele):
    return sele in cmd.get_names("selections");

def pseudo(residue, theta):
    print(theta)

    other=2*(sin(math.radians(36.0))+sin(math.radians(72.0)))
    
    #phase=atan2((theta4+theta1)-(theta3+theta0),2*theta2*(sin(math.radians(36.0))+sin(math.radians(72.0))))
    phase=atan2((theta['4']+theta['1'])-(theta['3']+theta['0']),theta['2']*other)
    amplitude=theta['2']/cos(phase)
    phase=math.degrees(phase)
    if phase < 0:
        phase+=360 # 0 <= Phase < 360
    #Determine pucker
    if phase < 36:
        pucker = "C3'-endo"
    elif phase < 72:
        pucker = "C4'-exo"
    elif phase < 108:
        pucker = "O4'-endo"
    elif phase < 144:
        pucker = "C1'-exo"
    elif phase < 180:
        pucker = "C2'-endo"
    elif phase < 216:
        pucker = "C3'-exo"
    elif phase < 252:
        pucker = "C4'-endo"
    elif phase < 288:
        pucker = "O4'-exo"
    elif phase < 324:
        pucker = "C1'-endo"
    elif phase < 360:
        pucker = "C2'-exo"
    else:
        pucker = "Phase is out of range"
    info = "%6.2f  %6.2f  %8s" % (phase, amplitude, pucker)
    return info
    

def store_atom(atom,residue):
    if atom.name == "O4'" or atom.name == "O4*":
        residue['O4*X'] = atom.coord[0]
        residue['O4*Y'] = atom.coord[1]
        residue['O4*Z'] = atom.coord[2]
    elif atom.name == "C1'" or atom.name == "C1*":
        residue['C1*X'] = atom.coord[0]
        residue['C1*Y'] = atom.coord[1]
        residue['C1*Z'] = atom.coord[2]
    elif atom.name == "C2'" or atom.name == "C2*":
        residue['C2*X'] = atom.coord[0]
        residue['C2*Y'] = atom.coord[1]
        residue['C2*Z'] = atom.coord[2]
    elif atom.name == "C3'" or atom.name == "C3*":
        residue['C3*X'] = atom.coord[0]
        residue['C3*Y'] = atom.coord[1]
        residue['C3*Z'] = atom.coord[2]
    elif atom.name == "C4'" or atom.name == "C4*":
        residue['C4*X'] = atom.coord[0]
        residue['C4*Y'] = atom.coord[1]
        residue['C4*Z'] = atom.coord[2]
    return

def get_dihedrals(residue,theta):

    C = []
    ribose = list(residue.keys())
    ribose.sort()
    
    shift_up(ribose,6)
    for i in range (0,12):
        C.append(residue[ribose[i]])
    theta['0']=calcdihedral(C)
    
    C = []
    shift_down(ribose,3)
    for i in range (0,12):
        C.append(residue[ribose[i]])
    theta['1']=calcdihedral(C)


    C = []
    shift_down(ribose,3)
    for i in range (0,12):
        C.append(residue[ribose[i]])
    theta['2']=calcdihedral(C)


    C = []
    shift_down(ribose,3)
    for i in range (0,12):
        C.append(residue[ribose[i]])
    theta['3']=calcdihedral(C)

    C = []
    shift_down(ribose,3)
    for i in range (0,12):
        C.append(residue[ribose[i]])
    theta['4']=calcdihedral(C)
    
    return

def shift_up(list,value):
    for i in range (0,value):
        list.insert(0,list.pop())
    return

def shift_down(list,value):
    for i in range (0,value):
        list.append(list.pop(0))
    return

def calcdihedral(C):
    
    DX12=C[0]-C[3]
    DY12=C[1]-C[4]
    DZ12=C[2]-C[5]
    
    DX23=C[3]-C[6]
    DY23=C[4]-C[7]
    DZ23=C[5]-C[8]
    
    DX43=C[9]-C[6];
    DY43=C[10]-C[7];
    DZ43=C[11]-C[8];

    #Cross product to get normal
    
    PX1=DY12*DZ23-DY23*DZ12;
    PY1=DZ12*DX23-DZ23*DX12;
    PZ1=DX12*DY23-DX23*DY12;

    NP1=sqrt(PX1*PX1+PY1*PY1+PZ1*PZ1);

    PX1=PX1/NP1
    PY1=PY1/NP1
    PZ1=PZ1/NP1

    PX2=DY43*DZ23-DY23*DZ43;
    PY2=DZ43*DX23-DZ23*DX43;
    PZ2=DX43*DY23-DX23*DY43;

    NP2=sqrt(PX2*PX2+PY2*PY2+PZ2*PZ2);
    
    PX2=PX2/NP2
    PY2=PY2/NP2
    PZ2=PZ2/NP2

    DP12=PX1*PX2+PY1*PY2+PZ1*PZ2

    TS=1.0-DP12*DP12

    if TS < 0:
        TS=0
    else:
        TS=sqrt(TS)
    
    ANGLE=math.pi/2.0-atan2(DP12,TS)

    PX3=PY1*PZ2-PY2*PZ1
    PY3=PZ1*PX2-PZ2*PX1
    PZ3=PX1*PY2-PX2*PY1

    DP233=PX3*DX23+PY3*DY23+PZ3*DZ23

    if DP233 > 0:
        ANGLE=-ANGLE

    ANGLE=math.degrees(ANGLE)

    if ANGLE > 180:
        ANGLE-=360
    if ANGLE < -180:
        ANGLE+=360

    return ANGLE

cmd.extend("pucker",pucker)
