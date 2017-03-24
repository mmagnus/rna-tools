"""Seq and secondary structure prediction"""

import subprocess

class Seq:
    def __init__(self, seq):
        self.seq = seq
        self.ss = ''
        self.ss_log = ''
        
    def __repr__(self):
        return self.seq

    def predict_ss(self, method="RNAfold"):
        """it creats /tmp/ss.fa and runs various methods for ss prediction"""
        # make tmp file
        f = open('/tmp/ss.fa','w')
        f.write('>test\n')
        f.write(self.seq)
        f.close()
        # run prediction
        if method == "RNAsubopt":
            from cogent.app.vienna_package import RNAfold, RNAsubopt
            r = RNAsubopt(WorkingDir="/tmp")
            res = r([self.seq])
            return str(res['StdOut'].read()).strip()

        if method == "ipknot":
            self.ss_log = subprocess.getoutput('ipknot /tmp/ss.fa')
            return '\n'.join(self.ss_log.split('\n')[2:])

        if method == "contextfold":
            cmd = "cd /home/magnus/work/opt/ContextFold_1_00 && java -cp bin contextFold.app.Predict in:" + self.seq
            self.ss_log = subprocess.getoutput(cmd)
            return '\n'.join(self.ss_log.split('\n')[1:])
        
        if method == "centroid_fold":
            self.ss_log = subprocess.getoutput('centroid_fold /tmp/ss.fa')
            return '\n'.join(self.ss_log.split('\n')[2:])

        if method == 'RNAfold':
            from cogent.app.vienna_package import RNAfold, RNAsubopt
            r = RNAfold(WorkingDir="/tmp")
            res = r([self.seq])
            self.ss_log = res['StdOut'].read()
            return self.ss_log.strip().split('\n')[-1].split()[0]

    def get_ss():
        if self.ss:
            return self.ss
        else:
            return self.predict_ss()
        
#main
if __name__ == '__main__':
    seq = Seq("CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG")
    print(seq.predict_ss(method="ipknot"))
    
