#!/usr/bin/python
import Bio.PDB 
import sys
import re


# Defines the usage of command line arguments associated with calls to
# functions in class Entropy
class Usage_PDBSeq:
  def __init__(self):
     # define a usage statement
     PROGRAM = 'get_pdb_seq'
     USAGE =  PROGRAM + ' -i <infile>.pdb [-split/-linear]\n'
     USAGE +='\n'
     USAGE +='flag             status        purpose                                                      \n'
     USAGE +='------         ---------       ------------------------------------------------------------ \n'
     USAGE +=' -i <f>        required        specify input PDB file name <f> \'pdb\' (required!)          \n' 
     USAGE +=' -split            -           produce a sequence [[\'M\',\'E\',\'T\',\'A\',\'L\'...], ...] \n' 
     USAGE +=' -linear           -           default [\'METAL....\', ...]                                 \n' 
     USAGE +=' -short            -           write as \'M\' instead of \'MET\', required for linear format\n'
     USAGE +=' -fasta            -           output the sequence part as a fasta file\n'
     self.u = USAGE
  #
  
  def get_Usage_PDBSeq(self):
     return 'Usage: %s\n' % self.u
  #
#

#
#  Reads the command line associated with programs using class Entropy
#
class PDBSeq_Opt:
  
  # NOTE: __everything__ is read or assumed to be a STRING here
  # whether it is a number or a string must be decided elsewhere.
  
  def __init__(self, command_line):
    self.command_line = command_line
    self.infile   = ''
    self.split = False
    self.short_form = False
    self.debug_PDBSeq_Opt = False
    self.fasta = False
  #
    
  # parses the command line arguments for Entropy calls
  def parse_PDBSeq_CL(self):
    #
    n = len(self.command_line)
    if (n == 0):
      print 'ERROR: command line undefined. '
      sys.exit(1)
    #
    flag_flnm = False
    split_cnt = 0
    i = 1
    while i < n:
      if self.debug_PDBSeq_Opt:
        print "index %d: \'%s\'" % (i, self.command_line[i])
      
      
      # is it the input file?
      if self.command_line[i] == '-i':
        #  input file
        i = self.check_end(i, n)
        self.infile = self.command_line[i]
        
        #
        if (self.infile == ''):
          emsg = "ERROR: input file name undefined."
          self.error(emsg)
          #
        #
        if self.debug_PDBSeq_Opt:
          print "infile name: %s" % self.infile 
        use_old_approach = False
        if use_old_approach:
          p = re.compile('.pdb$')
          a = p.findall(self.infile)
          # print len(a)
          if (len(a) == 0):
            emsg = 'ERROR: output file must have the extension \'.pdb\'.'
            self.error(emsg)
          #
        else:
          ext = self.infile[-3:]
          if not (ext == "pdb" or ext == "ent"):
            emsg = 'ERROR: output file must have the extension \'.pdb\' or \'.ent\'.'
            self.error(emsg)
            #
        #
        flag_flnm = True
           
      # should the output be in the format [['M', 'E', 'T', 'A', 'L', ....], .....] 
      elif self.command_line[i] == '-split':
        self.split = True
        split_cnt += 1
        #
           
        if self.debug_PDBSeq_Opt:
          print "Requested separated format [[\'M\',\'E\',\'T\',\'A\',\'L\'...], ...]"
           
        #
        #
        #
      # should the output be in the format ['METAL....', .....] 
      elif self.command_line[i] == '-linear':
        self.split = False
        self.short_form = True
        split_cnt += 1
        #
           
        if self.debug_PDBSeq_Opt:
          print "Requested linear format [\'METAL..\', ...]"
          print "short form of sequence required"
          
        #
        #
      #
      elif self.command_line[i] == '-short':
        self.short_form = True
        #
          
        if self.debug_PDBSeq_Opt:
          print "Requested short form for protein sequences"
        
        #
        #
      #
      elif self.command_line[i] == '-fasta':
        self.fasta = True
        #
          
        if self.debug_PDBSeq_Opt:
          print "Requested fasta formated sequence output"
        
        #
        #
      #
      # is it a plea for help? (in one of its many different versions ...)
      elif (self.command_line[i] == '-h') or \
           (self.command_line[i] == '--help') or \
           (self.command_line[i] == '-help'):
        # default
        a = Usage_PDBSeq()
        print '%s\n' % a.get_Usage_PDBSeq()
        sys.exit(0)
        #
      #
      else:
        emsg = 'ERROR: command line is corrupted:\n'
        emsg += self.get_CL()
        emsg += '\noffending command \'%s\'' % self.command_line[i]
        self.error(emsg)
        #
      #
      i += 1
    #
    # error checking
    # print "flag_flnm = ", flag_flnm
    if flag_flnm == False:
      emsg = "ERROR, no file defined."
      self.error(emsg)
    if split_cnt > 1:
      emsg = "ERROR, command line corrupted with more than one format request."
      self.error(emsg)
        
     
    if self.debug_PDBSeq_Opt:
      print 'infile:    %s' % self.infile
    return 0
  #

  def check_end(self, i, n):
     if (i+1 < n):
        i += 1
     else:
        emsg = "ERROR, insufficient number of command line arguments."
        self.error(emsg)
        # terminates
     return i
      

  # delivers error/exit messages for this class
  def error(self, e_msg):
     print '\n'
     print e_msg
     ss = 'input command: ' + self.get_CL() + '\n'
     print ss
     a = Usage_PDBSeq()
     print '%s\n' % a.get_Usage_PDBSeq()
     sys.exit(1)
     #
  #

  # separates the command line arguments into specific strings
  def get_CL(self):
     n = len(self.command_line)
     ss = ''
     for i in range(n):
        if (i < n-1):
           ss += '%s ' % self.command_line[i]
        else:
           ss += '%s' % self.command_line[i]
     return ss
  #
  #
  #
#



class PDBSeq:
   def __init__(self):
      self.molecule = []
      self.short_form = False
      self.split = False
      self.lens = []
      self.fasta = False
   #
   
   def set_short_form(self, b):
      self.short_form = b
      return 0
   
   def set_split(self, b):
      self.split = b
      return 0
   
   def set_fasta(self, b):
      self.fasta = b
      return 0
   
   #
   def get_seq(self, pdbflnm):
      # print "get_seq: short_form = %s, split = %s" % (self.short_form, self.split)
      # print pdbflnm
      p = Bio.PDB.PDBParser(QUIET=True)
      try:
         structure = p.get_structure('X', pdbflnm)
      except:
         print 'PDBseq: problems opening file %s' % pdbflnm
         sys.exit(1)
      
      # write the pdb file
      self.molecule = []
      conv = SeqFormat()
      for model in structure:
         # print 'model'
         for chain in model:
            ndx = 1
            cid = chain.get_id()
            # print 'chain'
            if self.split:
               seq_a = []
               for residue in chain:
                  res = conv.long2short(residue.get_resname().strip(), cid, ndx)
                  if not (res == ''):
                     seq_a += [res]
                  ndx += 1
            else:
               seq_a = ''
               for residue in chain:
                  res = conv.long2short(residue.get_resname().strip(), cid, ndx)
                  seq_a += res
                  ndx += 1
            #
            self.molecule += [seq_a]
         #
      #
      for i in range(0,len(self.molecule)):
        self.lens += [len(self.molecule[i])]
      output = ''
      if self.fasta:
        chains = ''
        for ch_i in range(0,len(self.molecule)):
          s = "> %s ch %d\n" % (pdbflnm, ch_i + 1)
          n = len(self.molecule[ch_i])
          sq_k = 0
          while sq_k < n:
            dd = min(60, n - sq_k)
            for kk in range(0,dd):
              s += self.molecule[ch_i][sq_k + kk]
            s += '\n'
            sq_k += 60
          # print s
          chains += s
          chains += '\n'
        #
        output = chains
      else:
        output = self.molecule
      return output
   #
#


# label formats, especially for proteins 
class SeqFormat:  
   def __init__(self):
      self.debug = False
   #
   
   def long2short(self, res, cid, ndx):
      r = ''
      if   res == 'ALA' or res == 'Ala':
         r = 'A'
      elif res == 'CYS' or res == 'Cys':
         r = 'C'
      elif res == 'ASP' or res == 'Asp':
         r = 'D'
      elif res == 'GLU' or res == 'Glu':
         r = 'E'
      elif res == 'PHE' or res == 'Phe':
         r = 'F'
      elif res == 'GLY' or res == 'Gly':
         r = 'G'
      elif res == 'HIS' or res == 'His':
         r = 'H'
      elif res == 'HID' or res == 'Hid':
         r = 'H'
      elif res == 'HIE' or res == 'Hie':
         r = 'H'
      elif res == 'HIP' or res == 'Hip':
         r = 'H'
      elif res == 'ILE' or res == 'Ile':
         r = 'I'
      elif res == 'LYS' or res == 'Lys':
         r = 'K'
      elif res == 'LEU' or res == 'Leu':
         r = 'L'
      elif res == 'MET' or res == 'Met':
         r = 'M'
      elif res == 'ASN' or res == 'Asn':
         r = 'N'
      elif res == 'PRO' or res == 'Pro':
         r = 'P'
      elif res == 'GLN' or res == 'Gln':
         r = 'Q'
      elif res == 'ARG' or res == 'Arg':
         r = 'R'
      elif res == 'SER' or res == 'Ser':
         r = 'S'
      elif res == 'THR' or res == 'Thr':
         r = 'T'
      elif res == 'VAL' or res == 'Val':
         r = 'V'
      elif res == 'TRP' or res == 'Trp':
         r = 'W'
      elif res == 'TYR' or res == 'Tyr':
         r = 'Y'
      else:
         # no change for Nucleic acids
         if not (res == 'A' or res == 'C' or res == 'G' or res == 'T' or res == 'U'):
            print 'WARNING: chain %s: residue %5d (\'%s\') not recognized' % (cid, ndx, res)
         else:
            r = res
      return r
   #
   def short2long(self, res, cid, ndx):
      r = ''
      if   res == 'A':
         r = 'ALA'
      elif res == 'C':
         r = 'CYS'
      elif res == 'D':
         r = 'ASN'
      elif res == 'E':
         r = 'GLU'
      elif res == 'F':
         r = 'PHE'
      elif res == 'G':
         r = 'GLY'
      elif res == 'H':
         r = 'HIS'
      elif res == 'I':
         r = 'ILE'
      elif res == 'K':
         r = 'LYS'
      elif res == 'L':
         r = 'LEU'
      elif res == 'M':
         r = 'MET'
      elif res == 'N':
         r = 'ASN'
      elif res == 'P':
         r = 'PRO'
      elif res == 'Q':
         r = 'GLN'
      elif res == 'R':
         r = 'ARG'
      elif res == 'S':
         r = 'SER'
      elif res == 'T':
         r = 'THR'
      elif res == 'V':
         r = 'VAL'
      elif res == 'W':
         r = 'TRP'
      elif res == 'Y':
         r = 'TYR'
      else:
         # no change for Nucleic acids
         if not (res == 'A' or res == 'C' or res == 'G' or res == 'T' or res == 'U'):
            print 'WARNING: chain %s: residue %5d (\'%s\') not recognized' % (cid, ndx, res)
         else:
            r = res
      return r
   #
#

class ResInfo:  
   def __init__(self):
      self.debug = False
   #
   def get_ResLen(self, res):
      len = 0
      if   res == 'ALA':
         len = 2
      elif res == 'CYS':
         len = 2
      elif res == 'ASP':
         len = 2
      elif res == 'GLU':
         len = 3
      elif res == 'PHE':
         len = 3
      elif res == 'GLY':
         len = 1
      elif res == 'HIS':
         len = 4
      elif res == 'ILE':
         len = 3
      elif res == 'LYS':
         len = 4
      elif res == 'LEU':
         len = 3
      elif res == 'MET':
         len = 3
      elif res == 'ASN':
         len = 2
      elif res == 'PRO':
         len = 2
      elif res == 'GLN':
         len = 3
      elif res == 'ARG':
         len = 4
      elif res == 'SER':
         len = 2
      elif res == 'THR':
         len = 3
      elif res == 'VAL':
         len = 3
      elif res == 'TRP':
         len = 6
      elif res == 'TYR':
         len = 5
      else:
         # no change for Nucleic acids
         if not (res == 'A' or res == 'C' or res == 'G' or res == 'T' or res == 'U'):
            print 'WARNING: residue \'%s\' not recognized' % res
            len = -1
         else:
            len = 5
         #
      #
      return len
   #
   
   def get_names(self, res):
      aname = []
      if   res == 'ALA':
         aname += ['CA']
         aname += ['CB']
      elif res == 'CYS':
         aname += ['CA']
         aname += ['CB']
      elif res == 'ASP':
         aname += ['CA']
         aname += ['CB']
      elif res == 'GLU':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'PHE':
         aname += ['CA']
         aname += ['CB']
         aname += ['CE']
      elif res == 'GLY':
         aname += ['CA']
      elif res == 'HIS':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
         aname += ['CE1']
         aname += ['CE2']
      elif res == 'ILE':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'LYS':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
         aname += ['CD']
         aname += ['CE']
      elif res == 'LEU':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'MET':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
         aname += ['CE']
      elif res == 'ASN':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'PRO':
         aname += ['CA']
         aname += ['CB']
      elif res == 'GLN':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'ARG':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
         aname += ['CE']
         aname += ['CZ']
      elif res == 'SER':
         aname += ['CA']
         aname += ['CB']
      elif res == 'THR':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'VAL':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
      elif res == 'TRP':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
         aname += ['CE']
         aname += ['CZ']
         aname += ['CT']
      elif res == 'TYR':
         aname += ['CA']
         aname += ['CB']
         aname += ['CG']
         aname += ['CE']
         aname += ['CZ']
      else:
         # no change for Nucleic acids
         if not (res == 'A' or res == 'C' or res == 'G' or res == 'T' or res == 'U'):
            print 'WARNING: residue \'%s\' not recognized' % res
         else:
            aname += ['P']
            aname += ['C4']
            aname += ['N']
            aname += ['B1']
            aname += ['B2']
         #
      #
      return aname
   #

#
