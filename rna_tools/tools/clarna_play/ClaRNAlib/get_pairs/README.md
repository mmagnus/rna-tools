Before you use this wrapper, you need to change CLARNA_PATH variable in get_pairs.py to reflect the actual location of ClaRNA on your machine.

Usage:

	from get_pairs import ClarnaGraph
	clarna_graph = ClarnaGraph("input_file.pdb", ["graph_file"])
	graph = clarna_graph.graph  # generated graph (a networkx's MultiDiGraph object)
	pairs = clarna_graph.get_pairs_from_graph()  # list of paired residues
	ss = clarna_graph.get_ss_from_graph()  # secondary structure as a two-element list

