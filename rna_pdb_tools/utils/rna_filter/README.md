    python rna_filter.py -r restraints.txt -s test_data/CG.pdb -v
    (d:A1-A2 <  10.0  1)|(d:A2-A1 <= 10 1)
     restraints [('A1', 'A2', '<', '10.0', '1'), ('A2', 'A1', '<=', '10', '1')]

    test_data/CG.pdb
     mb for  A1 [ 54.729   28.9375  41.421 ]
     mb for  A2 [ 55.3425  35.3605  42.7455]
      d:A1-A2 6.58677550096
      d:A2-A1 6.58677550096

    python rna_filter.py -r restraints.txt -t test_data/CG.trafl -v
    (d:A1-A2 <  10.0  1)|(d:A2-A1 <= 10 1)
     restraints [('A1', 'A2', '<', '10.0', '1'), ('A2', 'A1', '<=', '10', '1')]

    Frame #1 e:1252.26
      mb for A1 [ 54.729   28.9375  41.421 ]
      mb for A2 [ 55.3425  35.3605  42.7455]
       d:A1-A2 6.58677550096
      mb for A2 [ 55.3425  35.3605  42.7455]
      mb for A1 [ 54.729   28.9375  41.421 ]
       d:A2-A1 6.58677550096
