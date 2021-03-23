#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This module contains functions for computing RNA3DCNN potential

Output::

    Trainable params: 4,282,801
    Non-trainable params: 0
    _________________________________________________________________
    Scores for each nucleotide in /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmptg6jy2ud/query.pdb:
    [[ 0.02462262]
     [ 0.03271335]
     [ 0.06199259]
     [ 0.02006263]
     [ 0.05937254]
     [ 0.12025979]
     [ 0.20201728]
     [ 0.24463326]
     [ 0.43518737]
     [ 0.7260638 ]
     [ 0.6140108 ]
     [ 0.6588027 ]
     [ 0.7668936 ]
     [ 0.4776191 ]
     [ 0.39859247]
     [ 0.572009  ]
     [ 0.64892375]
     [ 0.11587611]
     [ 0.0560993 ]
     [ 0.05285829]
     [ 0.0167731 ]
     [ 0.01759553]
     [ 0.02143204]
     [-0.01818037]]
    Total score for /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmptg6jy2ud/query.pdb is  6.3262305

If missing atoms::

    Total params: 4,282,801
    Trainable params: 4,282,801
    Non-trainable params: 0
    _________________________________________________________________
    There is no atom O5' in residue 620A in chain  A in PDB /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpx87uus6x/query.pdb.
    There is no atom O5' in residue 635A in chain  B in PDB /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpx87uus6x/query.pdb.
    There is no atom O5' in residue 1750G in chain  C in PDB /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpx87uus6x/query.pdb.
    Scores for each nucleotide in /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpx87uus6x/query.pdb:
    []
    Total score for /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpx87uus6x/query.pdb is  0.0
   
"""
import os
import re
from shutil import copyfile
from rna_tools.tools.mq.lib.wrappers.SubprocessUtils import run_command
from rna_tools.tools.pdb_formatix.PDBFile import PDBFile#resname_check_and_3to1, set_residues_bfactor
from rna_tools.tools.mq.lib.wrappers.base_wrappers import ProgramWrapper
from rna_tools.rna_tools_config import RNA3DCNN_PATH, PYTHON3_PATH


class RNA3DCNN(ProgramWrapper):
    """
    Wrapper class for RNA3DCNN.
    """
    max_seq_len = 100000  # I don't know about any restriction

    def __init__(self):
        super(RNA3DCNN, self).__init__()

    def run(self, path_to_pdb, verbose=False):
        copyfile(path_to_pdb, self.sandbox_dir + os.sep + 'query.pdb')
        old_pwd = os.getcwd()

        os.chdir(self.sandbox_dir)

        self.log('start for %s' % self.sandbox_dir + '/query.pdb', level="debug")

        # 4. To print scores of each nucleotide and total scores, use flag "-local 1"
        # 5. To print only total scores, use flag "-local 0"
        # For example:<br />
        # python Main.py -pl pdblist -model RNA3DCNN_MD.hdf5 -local 0<br />
        cmd = PYTHON3_PATH + ' ' + RNA3DCNN_PATH + '/Main.py ' + \
        ' -pn ' + self.sandbox_dir + '/query.pdb ' + \
        ' -model ' + RNA3DCNN_PATH + '/RNA3DCNN_MD.hdf5 ' + \
        ' -local 1 2>>  '  + self.sandbox_dir + '/log.txt >>' + self.sandbox_dir + '/log.txt'
        if verbose:
            print(cmd)
        os.system(cmd)

        """
        Total params: 4,282,801
        Trainable params: 4,282,801
        Non-trainable params: 0
        _________________________________________________________________
        Total score for /var/folders/yc/ssr9692s5fzf7k165grnhpk80000gp/T/tmpO_jRVR/query.pdb is  6.3262305
        """
        self.log('Run finished')
        output = open(self.sandbox_dir + '/log.txt').read()
        score = output.split()[-1]
        lscore = re.search('\[\[.*\]\]', output, re.M | re.DOTALL)
        try:
            lscore = lscore.group().replace('[','').replace(']\n', ',').replace(']]','')
        except:
            print(output)
            return error
        # 0.02462262,  0.03271335,  0.06199259,  0.02006263,  0.05937254,  0.12025979, 
        lscore = [float(x.strip()) for x in lscore.split(',')] # [0.02462262, 0.03271335, ...
        if verbose:
            print(lscore)
            print(score)
        os.chdir(old_pwd)
        try:
            return float(score)
        except:
            error = 'Error: problem with the file'
            with open(self.sandbox_dir + '/log.txt') as f:
                error += f.read()
            return error

def main():
    wrapper = RNA3DCNN()
    try:
        result = wrapper.run('../test' + os.sep + '1a9n.pdb', verbose=True)
        if result:
            print(result)
    except Exception as e:
        print(e)
    finally:
        #wrapper.cleanup()
        pass

if '__main__' == __name__:
    main()
