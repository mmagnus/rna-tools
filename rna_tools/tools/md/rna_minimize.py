#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run for parallel:

    parallel rna_minimize.py {} ::: *mdr.pdb

"""
from __future__ import print_function
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except:
    print('pip intsall openmm')

from sys import stdout
import argparse
from rna_tools.tools.mq.lib.timex import timex


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    parser.add_argument("-b", "--box-size", help="", default=1, type=float) # nargs='+')
    parser.add_argument("-sp", "--solv-padding", action="store_true")
    parser.add_argument("--score",
                        action="store_true", help="be verbose")
    parser.add_argument("--pymol",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    score = args.score

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        if not score:
            print(f, '...')
            t = timex.Timex()
            t.start()

        pdbout = f.replace('.pdb','') + '_min.pdb'
        pdb = PDBFile(f)
        log = f.replace('.pdb','') + '.log'
        modeller = Modeller(pdb.topology, pdb.positions)

        #ff = 'charmm36.xml' #ff14SB.xml' #amber14sb.xml' # 'amber14-all.xml'
        ff = 'amber14-all.xml'
        #ff = 'amberfb15.xml'
        #ff = 'amber14/RNA.OL3.xml'
        #ff = 'amber99sb.xml'

        forcefield = ForceField(ff, 'amber14/tip3pfb.xml')
        modeller.addHydrogens(forcefield)
        #modeller.addSolvent(forcefield, ionicStrength=0.1*molar)
        # modeller.addSolvent(forcefield, model='tip5p')

        bs = args.box_size
        args.solv_padding = True
        if args.solv_padding:
            # print(1*nanometers)
            modeller.addSolvent(forcefield, padding=1*nanometers)
        else:
            modeller.addSolvent(forcefield, boxSize=Vec3(bs, bs, bs)*nanometers)
        # 5.0, 3.5, 3.5
        # modeller.addSolvent(forcefield, boxSize=Vec3(2, 2, 2)*nanometers)

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, #nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
        #simulation.reporters.append(PDBReporter("output.pdb", 1))
        simulation.minimizeEnergy(maxIterations=1000)#, verbose=True)#1.0, verbose=True) #)
        # from http://zarbi.chem.yale.edu/ligpargen/openMM_tutorial.html
        position = simulation.context.getState(getPositions=True).getPositions()
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        app.PDBFile.writeFile(simulation.topology, position,
                          open(pdbout, 'w'))
        if score:
            #sprint(energy)
            print(f + ',' + str(energy._value * KcalPerKJ))
        else:
            print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
            print('saved ', pdbout)
            if args.pymol:
                os.system('open %s' % out)
            print(t.end())
