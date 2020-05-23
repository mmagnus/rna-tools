#!/usr/bin/env python
"""
A simple wrapper that returns ClaRNA results as a list of paired residues
and/or secondary structure in dot-bracket format.
"""
import os, sys

CLARNA_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, CLARNA_PATH)
gtss = __import__("graph-to-secondary-structure")

import clarna
import utils
import networkx as nx
import tempfile
import subprocess

class DummyOptions(object):
    pdb = ""
    input = ""

dummy = DummyOptions


class ClarnaException(Exception):
    pass


class ClarnaGraph(object):
    def __init__(self, pdb_fname, graph_fname=None):
        if graph_fname is None:
            temp_file = tempfile.NamedTemporaryFile(prefix="clarna")
            graph_fname = temp_file.name
        self.graph_fname = self.make_graph(pdb_fname, graph_fname)
        self.graph = utils.load_graph(self.graph_fname)

    def make_graph(self, pdb_fname, graph_fname):
        """
        Executes ClaRNA and generates a graph file from an input .pdb file.
        Returns the name of the generated file.
        """
        command = (os.path.join(CLARNA_PATH, "clarna.py"),
                   "--save-graph=" + graph_fname,
                   pdb_fname)
        exitcode = subprocess.call(command, stdout=open(os.devnull, "w"))
        if exitcode != 0:
            raise ClarnaException
        return graph_fname

    def get_pairs_from_graph(self):
        """
        Converts the graph to a list of paired residue IDs.
        """
        nodes = gtss.get_nodes_ordering(self.graph, dummy)
        nodes2idx = dict([(n, i) for i, n in enumerate(nodes)])
        gg = nx.DiGraph()
        canonical_pairs = ('CG', 'GC', 'AU', 'UA')
        for u, v, data in self.graph.edges(data=True):
            if data.get('desc') == 'WW_cis' and data.get('n_type')\
            in canonical_pairs:
                gg.add_edge(u, v)
        brackets = ['.' for n in nodes]
        return nx.algorithms.maximal_matching(gg)

    def get_ss_from_graph(self):
        """
        Returns the sequence and dot-bracket secondary structure.
        """
        ss = gtss.gen_secondary_structure(self.graph, dummy)
        if len(ss) > 3:
            return ss.split("\n")[1:-1]
        return None


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Please specify input file"
        exit()

    gr_file = sys.argv[2] if len(sys.argv) > 2 else None

    clarna_graph = ClarnaGraph(sys.argv[1], gr_file)
    print "pairs:"
    print clarna_graph.get_pairs_from_graph()
    print "secondary structure:"
    print clarna_graph.get_ss_from_graph()

