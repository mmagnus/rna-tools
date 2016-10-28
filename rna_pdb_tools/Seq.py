from cogent.app.vienna_package import RNAfold

class Seq:
    def __init__(self, seq):
        self.seq = seq
        self.ss = ''
        self.ss_log = ''
        
    def __repr__(self):
        return self.seq
    def predict_ss(self):
        r = RNAfold(WorkingDir="/tmp")
        res = r([self.seq])
        self.ss_log = res['StdOut'].read()
        return self.ss_log.strip().split('\n')[-1].split()[0]

if __name__ == '__main__':
    seq = Seq("CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG")
    print seq.get_ss()
    
