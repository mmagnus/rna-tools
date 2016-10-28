import subprocess
import os

class RfamSearch():
    def __init__(self):
        pass

    def cmscan(self, seq):
        print seq
        # make tmp file
        f = open('/tmp/ss.fa','w')
        f.write('>test\n')
        f.write(seq.seq)
        f.close()

        old_pwd = os.getcwd()
        os.chdir('/home/magnus/work/rfamdb')
        cmd = 'cmscan -E 1 Rfam.cm /tmp/ss.fa > /tmp/cmscan.txt'
        subprocess.Popen(cmd, shell=True)
        self.output = open('/tmp/cmscan.txt').read()
        os.chdir(old_pwd)
        return self.output
#main
if __name__ == '__main__':
    import Seq
    seq = Seq.Seq("GGCGCGGCACCGUCCGCGGAACAAACGG")
    rs = RfamSearch()
    rs.cmscan(seq)
