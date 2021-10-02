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
        pdb = PDBFile(f)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, #nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print(f, energy)
