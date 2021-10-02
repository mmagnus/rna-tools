#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    parser.add_argument("--pymol",
                        action="store_true", help="be verbose")
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        print(f, '...')
        pdbout = f.replace('.pdb','') + '_min.pdb'
        pdb = PDBFile(f)
        log = f.replace('.pdb','') + '.log'
        modeller = Modeller(pdb.topology, pdb.positions)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller.addHydrogens(forcefield)
        #modeller.addSolvent(forcefield, ionicStrength=0.1*molar)
        # modeller.addSolvent(forcefield, model='tip5p')
        #modeller.addSolvent(forcefield, padding=0.5*nanometers)
        modeller.addSolvent(forcefield, boxSize=Vec3(5.0, 3.5, 3.5)*nanometers)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff, #nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy()
        # from http://zarbi.chem.yale.edu/ligpargen/openMM_tutorial.html
        position = simulation.context.getState(getPositions=True).getPositions()
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        app.PDBFile.writeFile(simulation.topology, position,
                          open(pdbout, 'w'))
        print('Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ))
        print('saved ', pdbout)
        if args.pymol:
            os.system('open %s' % out)
