"""Install by adding the Python file to you ~/.pymolrc, e.g.

    run /Users/magnus/work/src/rna-tools/rna_tools/tools/mq/RASP/rasp_for_pymol.py
    
Usage: 

    rasp all # all is selection here
    
"""

from pymol import cmd
import os
import subprocess
from tempfile import gettempdir

from rna_tools.rna_tools_config import RASP_PATH, RASP_PATH_DIR

def rasmol2pymol(profile_file):
    f = open(profile_file)
    old_lines = f.readlines()
    new_lines = []
    line = ''
    res_num = 0
    for l in old_lines:
        if l.startswith('color'):
            color = eval(l.split()[1])
            color = [float(i)/256 for i in color]
            colname = 'c' + str(color).replace(' ', '').replace('[', '').replace(']', '').replace(',', '').replace('.', '')
            line = 'set_color %s, [%f , %f , %f]' % (colname, color[0], color[1], color[2])
            new_lines.append(line)
            line = 'color ' + colname + ', resi ' + res_num
            new_lines.append(line)
        elif l.startswith('select'):
            res_num = l.split()[1]
        line = ''
    return new_lines

def rasp(selection):
    tempdir = gettempdir()
    cmd.save(os.path.join(tempdir, 'rasp_input.pdb'), selection, format='pdb')
    cmd1 = (
            os.path.join(RASP_PATH, '', 'bin', 'rasp_fd') + 
            ' -e all -p ' + os.path.join(tempdir, 'rasp_input.pdb')
            ).split()
    print(' '.join(cmd1))
    cmd2 = (
            os.path.join(RASP_PATH, '', 'bin', 'rasp_profile_fd') +
            ' -e all -r -p ' + os.path.join(tempdir, 'rasp_input.pdb')
            ).split()

    rasp_dir = os.path.join(RASP_PATH, '')
    curr_dir = os.getcwd()
    os.chdir(tempdir)
    print(tempdir)
    subprocess.Popen(cmd1, env={'RASP': rasp_dir})
    subprocess.Popen(cmd2, env={'RASP': rasp_dir})
    print(' '.join(cmd1))
    print(' '.join(cmd2))
    os.chdir(curr_dir)
    cmds = rasmol2pymol(tempdir + '/profile.scr')
    for c in cmds:
        print(c)
        cmd.do(c)



cmd.extend("rasp", rasp)
