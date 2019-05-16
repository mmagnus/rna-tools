#!/usr/bin/env python
from rna_tools.tools.simrna_trajectory.simrna_trajectory import SimRNATrajectory
s = SimRNATrajectory()
s.load_from_file('rp14_aa22-6d8fb934_ALL.trafl', top_level=True)
s.plot_energy('plot1.png')
s.sort()
s.plot_energy('plot2.png')
