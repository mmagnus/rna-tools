#!/usr/bin/env python
from get_pairs import ClarnaGraph

clarna_graph = ClarnaGraph("GAAA_tetraloop.pdb")
print "\n".join(clarna_graph.get_ss_from_graph())

