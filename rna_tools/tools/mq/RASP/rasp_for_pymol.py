"""Install by adding the Python file to you ~/.pymolrc, e.g.

    run /Users/magnus/work/src/rna-tools/rna_tools/tools/mq/RASP/rasp_for_pymol.py
    
"""

from pymol import cmd
import os
import subprocess
from tempfile import gettempdir
RASP_PATH = "/Users/magnus/work/src/rna-tools/opt/" # rasp-fd-1.0/bin"
PROFILE_FILE = '/Users/magnus/Desktop/profile.scr'

def rasmol2pymol():
    f = open(PROFILE_FILE)
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
            os.path.join(RASP_PATH, 'rasp-fd-1.0', 'bin', 'rasp_fd') + 
            ' -e all -p ' + os.path.join(tempdir, 'rasp_input.pdb')
            ).split()
    print(' '.join(cmd1))
    cmd2 = (
            os.path.join(RASP_PATH, 'rasp-fd-1.0', 'bin', 'rasp_profile_fd') +
            ' -e all -r -p ' + os.path.join(tempdir, 'rasp_input.pdb')
            ).split()

    rasp_dir = os.path.join(RASP_PATH, 'rasp-fd-1.0')
    subprocess.Popen(cmd1, env={'RASP': rasp_dir})
    subprocess.Popen(cmd2, env={'RASP': rasp_dir})
    print(' '.join(cmd1))
    print(' '.join(cmd2))
    cmds = rasmol2pymol()
    for c in cmds:
        print(c)
        cmd.do(c)



cmd.extend("rasp", rasp)
