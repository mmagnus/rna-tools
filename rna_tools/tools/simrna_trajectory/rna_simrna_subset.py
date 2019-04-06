from .simrna_trajectory import SimRNATrajectory

if __name__ == '__main__':
    fn_trafl = 'test_data/mini.trafl'
    fn = (line for line in open(fn_trafl).xreadlines())
    c = 0
    while 1:
        try:
            header = fn.next().strip()
        except StopIteration:  # not nice
            break
        c += 1
        coords = fn.next().strip()

        traj = SimRNATrajectory()
        traj.load_from_string(c, header + '\n' + coords)
        f = traj.frames[0]
        print(f.header)
        idx = range(13, 16 + 1) + range(34, 36 + 1) + range(45, 49 + 1)
        for r in [f.residues[i] for i in idx]:
            # print 'atoms of ', r, r.get_atoms_txt()
            print(r.get_atoms_txt(),)
        print
