#!/usr/bin/env python
import sys
import os
import re
import gzip
from optparse import OptionParser
# biopython
from Bio import PDB

# TODO: normalize atom names, ie.
# O1P => OP1, O2P => OP2, itd.
# O3* => O3', itd.

def parse_args():
    """setup program options parsing"""
    parser = OptionParser(description="normalize PDB file (change residue numbers, chains, reverse)")
    parser.add_option("-o", "--output", dest="output",
                  help="save output to FILE", metavar="FILE")
    parser.add_option("-i", "--input", dest="input",
                  help="read input from FILE", metavar="FILE")
    parser.add_option("--reverse", dest="reverse", action="store_true",
                  help="reverse the order of the residues",
                  default=False)
    parser.add_option("--reverse-id", dest="reverse_id", action="store_true",
                  help="reverse the id of the residues",
                  default=False)
    parser.add_option("--residues-limit", dest="residues_limit",
                  help="save only N residues", metavar="N")
    (options, args)  = parser.parse_args()
    return (parser, options, args)

def reverse_structure(s):
    new_s = PDB.Structure.Structure("reversed")
    for model in s:
        new_model = PDB.Model.Model(model.id)
        new_chains = []
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            residues = []
            for r in chain:
                residues.append(r)
            residues.reverse()
            for r in residues:
                new_chain.add(r)
            new_chains.append(new_chain)
        new_chains.reverse()

        for new_chain in new_chains:
            new_model.add(new_chain)
        new_s.add(new_model)

    return new_s


def normalize_pdb(inp_fn, out_fn, reverse=False, reverse_id=False, limit=None):
    """generate out_fn and returns array of tuples
       (old_chain,old_res_id), (new_chain, new_res_id)"""
    numbers = []

    parser = PDB.PDBParser()
    if re.match("^.*.gz",inp_fn):
        inp = gzip.open(inp_fn)
    else:
        inp = open(inp_fn)
    s = parser.get_structure("c",inp)
    inp.close()
    total = 0
    for model in list(s):
        n = len([r for r in model.get_residues()])
        for chain in list(model):
            for i, r in enumerate(list(chain)):
                total += 1
                if limit is not None and total>limit:
                    chain.detach_child(r.get_id())
                    continue
                old_id = r.get_id()
                old_num = (str(old_id[1])+str(old_id[2])).strip()
                old_chain = r.get_parent().get_id()
                if reverse_id:
                    new_num = n-i
                else:
                    new_num = i
                new_id = (old_id[0], new_num, ' ') # we ignore old_id[2]
                numbers.append(((old_chain, old_num), (old_chain, str(new_num))))
                # print old_id, new_id
                r.id = new_id
            if len(chain)==0:
                model.detach_child(chain.id)
        if len(model)==0:
            s.detach_child(model.id)
    if reverse:
        s = reverse_structure(s)
    if limit is not None:
        c = len(list(s.get_residues()))
        print("number of residues=%s" % c)
        assert c<=limit

    # from http://biopython.org/wiki/Remove_PDB_disordered_atoms
    class NotDisordered(PDB.Select):
        def accept_atom(self, atom):
            return not atom.is_disordered() or atom.get_altloc()=='A'

    io = PDB.PDBIO()
    io.set_structure(s)
    io.save(out_fn, select=NotDisordered())

    return numbers

def main():
    (parser, options, args) = parse_args()
    if options.input and options.output:
        limit = None
        if options.residues_limit is not None:
            limit = int(options.residues_limit)
        normalize_pdb(options.input, options.output, options.reverse, options.reverse_id, limit=limit)
    else:
        print("specify input and output")
        parser.print_help()
        exit(1)

if __name__ == '__main__':
    main()

