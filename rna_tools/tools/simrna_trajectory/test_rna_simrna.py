from simrna_trajectory import SimRNATrajectory


def test():
    s = SimRNATrajectory()
    s.load_from_file('test_data/mini.trafl', top_level=True)
    s.plot_energy('plot.png')
    print('OK')

if __name__=="__main__":
    test()
