#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
::

    rna_md.py -s 100 -n 500 --box-size 8 post2ndXXXXXXXXX_min.pdb # mac 4:13

hz::

    rna_md.py -s 100 -n 200 --box-size 8 1st-m-U57a_rpr_min.pdb
    1st-m-U57a_rpr_min.pdb ...
    #"Progress (%)","Step","Potential Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Time Remaining"
    50.0%,100,-939919.8341049199,-916122.9981389735,57.542401419204666,--
    100.0%,200,-921820.8285070366,-879767.8151791963,101.68710568331676,0:00
    saved  1st-m-U57a_rpr_min_b8.0_MD.pdb
    Cmd: rna_md.py -s 100 -n 200 --box-size 8 1st-m-U57a_rpr_min.pdb

notes::

     20000 / 100 =  200
     10000 / 100 =  100
    100000 / 100 = 1000

    # 500 000 * 0.002 ps = 1 ns      # * 1000 * 500 -> 1 ns
    #  50 000 * 0.002 ps = 0.1 ns
    #  30000

"""
from __future__ import print_function
from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout
import argparse
from rna_tools.tools.mq.lib.timex import timex


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', "--steps", type=int, help="", default=100)
    parser.add_argument('-n', "--nsim", type=int, help="", default=500000)
    # 500000 to get 1 ns
    # 50000 
    parser.add_argument('-r', "--run", help="", default='_')
    parser.add_argument('-p', "--prefix", help="", default='')
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--dcd",
                        action="store_true", help="save also as dcd")
    parser.add_argument("file", help="", default="") # nargs='+')

    parser.add_argument("--min",
                        action="store_true", help="be verbose")
    parser.add_argument("--pymol",
                        action="store_true", help="be verbose")
    parser.add_argument("--box-size", help="", default=1, type=float) # nargs='+')
    parser.add_argument("--solv-padding",
                        action="store_true", help="be verbose")

    return parser


if __name__ == '__main__':
    # print(0.002*picoseconds * 500000, '1000 ps is 1 ns :P')

    parser = get_parser()
    args = parser.parse_args()

    # / 1000 ps is 1 ns :P
    print('your sim is', 0.002*picoseconds * args.nsim, round((0.002*picoseconds * args.nsim) / nanoseconds, 2) * 100, '% of 1 ns')
    import time
    time.sleep(3)

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        print(f, '...')
        t = timex.Timex()
        t.start()

        print(args)

        pdb = PDBFile(f)
        modeller = Modeller(pdb.topology, pdb.positions)
        ff = 'ff14SB.xml' #amber14sb.xml' # 'amber14-all.xml'
        ff = 'amberfb15.xml'#amber14-all.xml'
        ff = 'amber14-all.xml'
        forcefield = ForceField(ff, 'amber14/tip3pfb.xml') # 'tip3p.xml') #

        modeller.addHydrogens(forcefield)
        #modeller.addSolvent(forcefield, ionicStrength=0.1*molar)
        # modeller.addSolvent(forcefield, model='tip5p')

        if args.solv_padding:
            modeller.addSolvent(forcefield, padding=0.5*nanometers)
            out = f.replace('.pdb','') + '_padd0.5' + '_MD.pdb'
            log = f.replace('.pdb','') + '_padd0.5' + '_log.csv'

        elif args.box_size:
            bs = args.box_size
            modeller.addSolvent(forcefield, boxSize=Vec3(bs, bs, bs)*nanometers) # nanometers)
            b = '_b' + str(bs)
            out = f.replace('.pdb','') + b + '_MD.pdb'
            log = f.replace('.pdb','') + b + '_log.csv'
        else:
            #print('boxSize=Vec3(5.0, 3.5, 3.5)*nanometers')
            modeller.addSolvent(forcefield, boxSize=Vec3(5.0, 3.5, 3.5)*nanometers)
            out = f.replace('.pdb','') + '_MD.pdb'
            log = f.replace('.pdb','') + '_log.csv'

        out = out.replace('.pdb', '_s%sn%s' % (args.steps, args.nsim) + '.pdb')
        if args.prefix:
            out = out.replace('.pdb', '_' + args.prefix + '.pdb')

        log = out.replace('.pdb', '.csv')

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, #nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        print(0.002*picoseconds * args.nsim)
        
        if args.min:
            simulation.minimizeEnergy()

            position = simulation.context.getState(getPositions=True).getPositions()
            energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
            pdbout = f.replace('.pdb', '_min.pdb')
            app.PDBFile.writeFile(simulation.topology, position,
                              open(pdbout, 'w'))
            print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
            print('saved ', pdbout)

        nstep = args.steps

        if args.dcd:
            simulation.reporters.append(DCDReporter(out.replace('.pdb', '.dcd'), nstep))
        simulation.reporters.append(PDBReporter(out, nstep))#  keepIds=True))
        
        simulation.reporters.append(StateDataReporter(stdout, nstep, step=True,
                potentialEnergy=True, temperature=True, progress=True, totalSteps=args.nsim, remainingTime=True, totalEnergy=True))
        simulation.reporters.append(StateDataReporter(log, nstep, step=True,
                potentialEnergy=True, temperature=True, progress=True, totalSteps=args.nsim, remainingTime=True,
                                                      totalEnergy=True))
        simulation.step(args.nsim)
        print('saved ', out)
        if args.pymol:
            os.system('open %s' % out)

        print(t.end())
