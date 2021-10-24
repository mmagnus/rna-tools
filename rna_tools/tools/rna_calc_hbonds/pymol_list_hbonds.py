# Copyright (c) 2010 Robert L. Campbell
# @mmagnus this is not my code
# taken from http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/list_hbonds.py
from pymol import cmd

def list_hbonds(selection,selection2=None,cutoff=3.2,angle=45,mode=1,hb_list_name='hbonds'):
  """
  USAGE

  list_hb selection, [selection2 (default=None)], [cutoff (default=3.2)],
                     [angle (default=45)], [mode (default=1)],
                     [hb_list_name (default='hbonds')]

  The script automatically adds a requirement that atoms in the
  selection (and selection2 if used) must be either of the elements N or
  O.

  If mode is set to 0 instead of the default value 1, then no angle
  cutoff is used, otherwise the angle cutoff is used and defaults to 45
  degrees.

  e.g.
  To get a list of all H-bonds within chain A of an object
    list_hb 1abc & c. a &! r. hoh, cutoff=3.2, hb_list_name=abc-hbonds

  To get a list of H-bonds between chain B and everything else:
    list_hb 1tl9 & c. b, 1tl9 &! c. b

  """
  cutoff=float(cutoff)
  angle=float(angle)
  mode=float(mode)
# ensure only N and O atoms are in the selection
  selection = selection + " & e. N+O"
  if not selection2:
    hb = cmd.find_pairs(selection,selection,mode=mode,cutoff=cutoff,angle=angle)
  else:
    selection2 = selection2 + " & e. N+O"
    hb = cmd.find_pairs(selection,selection2,mode=mode,cutoff=cutoff,angle=angle)

# sort the list for easier reading
  hb.sort(key=lambda x:x[0][1])

  for pairs in hb:
    cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'print("%1s/%3s`%s/%-4s " % (chain,resn,resi,name),end=" ")')
    cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'print("%1s/%3s`%s/%-4s " % (chain,resn,resi,name),end=" ")')
    print("%.2f" % cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1])))

cmd.extend("hbl", list_hbonds)
