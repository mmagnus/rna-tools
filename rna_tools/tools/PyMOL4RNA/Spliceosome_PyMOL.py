import matplotlib.pyplot as plt  # to start
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
import numpy as np
import sys

"""version when the same protein gets SMD2_2 and SMD2_1"""

f = "/Users/magnus/work-src/rna-tools/rna_tools/tools/PyMOL4RNA/Spliceosome_PyMOL.xlsx"
df = pd.read_excel(f)

strucs = [
          ['S1_B_5zwo', '5zwo', '_B5zwo'],
          ['S2_Bact_5gm6', '5gm6', '_Ba5gm6'],
          ['S3_C_5lj3', '5lj3', '_C5lj3'],
          ['S4_Cstar_5mps', '5mps', '_Cs5mps'],
          ['S5_P_6exn', '6exn', '_P6exn'],
          ['S5_P_5ylz', '5ylz', '_P5ylz'],
          ['S6_ILS_5y88', '5y88', '_I5y88'],
          ['SX_pX_3jb9', '3jb9', '_3jb9'],
          ['S6_hP_6icz', '6icz', '_hP_6icz'],
          ['hBa_6ff7', '6ff7', '_hBa_6ff7'],
          ['hC_5yzg', '5yzg', '_hC_5yzg'],
          ['hX_5xjc', '5xjc', '_hX_5xjc'],
          ['yCs_5mq0', '5mq0', '_yCs_5mq0'],
          ['yC_5lj5', '5lj5', '_yC_5lj5'],
          ['yPre_6g90', '6g90', '_yPre_6g90'],
          ['yE_6n7r', '6n7r', '_yE_6n7r'],
          ['yE_6n7p', '6n7p', '_yE_6n7p'],
         ]


strs = []
for s in strucs:
    strs.append(s[1])
strs = str(strs)    #strs = "['5gm6', '5zwo', '5lj3', '5gm6', '5mps', '6exn', '5ylz', '5y88', '3jb9', '6icz']"

#strucs = [['S3_C_5lj3', '5lj3']] # for debugging
txt = """
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

def spl_extract():
 for name in cmd.get_names("all"):
    #if name in """ + strs + """: # this should be auto
    print(" \ Extracting mode for %s" % name)
"""

for struc_col, pdb_id, complex_code in strucs:
    txt += "    if '" + pdb_id + "' in name.lower():\n"
    for index, i in df.iterrows():
        #print(i)
        ref = str(i['chain']).strip()
        struc = str(i[struc_col]).strip()
        # print(ref, struc)
        # todo fix  cmd.extract("CDC40 (PRP17, SLU4, XRS2)", "chain n")
        i['protein'] = i['protein'].split()[0]
        if struc != 'nan' and struc != '-': # if ref != 'nan' # i dont' use reference chain any more
            if ',' in struc:
                subchains = struc.split(',') # a,P,h
                for index , sc in enumerate(subchains):
                    sc = sc.strip()
                    # cmd.extract(name, selection, source_state, target_state)
                    txt += '        cmd.extract("' + i['protein'] + '_' + str(index + 1)                     + complex_code + '", "chain ' + sc + ' and ' + pdb_id + '")\n'
                # one mode to see these two proteins
                #print('        cmd.extract("' + i['protein'] + '", "chain ' + '+'.join(subchains) + '")')
            else:
                #if ref != struc:
                txt += '        cmd.extract("' + i['protein'] + complex_code + '", "chain ' + struc + ' and ' + pdb_id + '")\n'
    # the rest make unkown at the very end
    txt += '        cmd.set_name("' + pdb_id + '", "unknown_other' + complex_code + '")\n'
    txt += '        cmd.group("' + complex_code.replace('_', '') + '", "*' + complex_code + '")\n'
    txt += '        cmd.do("order *, yes")\n'
print(txt)
with open('/home/magnus/work-src/rna-tools/rna_tools/tools/PyMOL4RNA/code_for_spl.py', 'w') as f:
    f.write(txt)


# In[113]:


txt = """
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

def spl_color():
 for name in cmd.get_names("all"):
  # cmd.do('color grey50') # off gray
  if name in """ + strs + """: # this should be auto
    print(" \ Extracting mode for %s" % name)
"""

for struc_col, pdb_id, complex_code in strucs:
    txt += "    if '" + pdb_id + "' in name.lower():\n"
    for index, i in df.iterrows():
        ref = str(i['chain']).strip()
        struc = str(i[struc_col]).strip()
        # print(ref, struc)
        # todo fix  cmd.extract("CDC40 (PRP17, SLU4, XRS2)", "chain n")
        i['protein'] = i['protein'].split()[0]
        if struc != 'nan' or struc != '-': # if ref != 'nan' # i dont' use reference chain any more
            if ',' in struc:
                subchains = struc.split(',') # a,P,h
                for index , sc in enumerate(subchains):
                    sc = sc.strip()
                    # cmd.extract(name, selection, source_state, target_state)
                    txt += '        cmd.do("color ' + i['color'] +                     ', chain ' + sc + ' and ' + pdb_id + '")\n'
                # one mode to see these two proteins
                #print('        cmd.extract("' + i['protein'] + '", "chain ' + '+'.join(subchains) + '")')
            else:
                #if ref != struc:
                txt += '        #' + i['protein'] +'\n'
                try:
                    txt += '        cmd.do("color ' + i['color'] + ', chain ' + struc + ' and ' +                 pdb_id + '")' + '\n'
                except TypeError:
                    print('Color entry is missing for ' + i['protein'])
                    sys.exit(1)

                # color <object>
                txt += '        cmd.do("color ' + i['color'] + ', ' + i['protein'] +                 complex_code + '")\n'
print(txt)

with open('/home/magnus/work-src/rna-tools/rna_tools/tools/PyMOL4RNA/code_for_color_spl.py', 'w') as f:
    f.write(txt)

#strucs = [['S3_C_5lj3', '5lj3']] # for debugging
txt = """
try:
    from pymol import cmd
except ImportError:
    print("PyMOL Python lib is missing")
    # sys.exit(0)

def spl_color():
  #cmd.do('color grey50') # PRP8
"""

for struc_col, pdb_id, complex_code in strucs:
    txt += '  if True: # fake if, just a quick hack\n'
    for index, i in df.iterrows():
        ref = str(i['chain']).strip()
        struc = str(i[struc_col]).strip()
        # print(ref, struc)
        # todo fix  cmd.extract("CDC40 (PRP17, SLU4, XRS2)", "chain n")
        i['protein'] = i['protein'].split()[0]
        if struc != 'nan' and struc != '-': # if ref != 'nan' # i dont' use reference chain any more
            if ',' in struc:
                subchains = struc.split(',') # a,P,h
                for index , sc in enumerate(subchains):
                    sc = sc.strip()
                    # cmd.extract(name, selection, source_state, target_state)
                    txt += '        cmd.do("color ' + i['color'] + ', ' + i['protein'] + '_' + str(index + 1) +                     complex_code + '")\n'
                # one mode to see these two proteins
                #print('        cmd.extract("' + i['protein'] + '", "chain ' + '+'.join(subchains) + '")')
            else:
                #if ref != struc:
                txt += '        cmd.do("color ' + i['color'] + ', ' + i['protein'] + complex_code + '")\n'
print(txt)
with open('/home/magnus/work-src/rna-tools/rna_tools/tools/PyMOL4RNA/code_for_color_spl_objects.py', 'w') as f:
    f.write(txt)
