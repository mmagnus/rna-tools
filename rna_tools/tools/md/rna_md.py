#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

    # 500000 # 0.002 ps * 1000 * 500 -> 1 ns
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
    parser.add_argument('-n', "--nsim", type=int, help="", default=50000) # 500000)
    parser.add_argument('-r', "--run", help="", default='_')
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
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
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        print(f, '...')
        t = timex.Timex()
        t.start()

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
            out = f.replace('.pdb','') + 'padd0.5' + '_MD.pdb'
            log = f.replace('.pdb','') + 'padd0.5' + '_log.csv'

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

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, #nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

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
