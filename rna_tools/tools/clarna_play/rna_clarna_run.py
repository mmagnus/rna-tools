#!/usr/bin/env python
# this is a basic way to access the clarna program with a minimum of
# pain.

PROGRAM = "run_clarna.py"

import sys, os, re
import urllib.request, urllib.error, urllib.parse, urllib.request, urllib.parse, urllib.error
import string
from Bio import SeqIO
from optparse import OptionParser

####### ClaRNA related imports  ########
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
########################################
from Bio import PDB

from rna_tools.tools.clarna_play.ClaRNAlib.structure_ciach import StructureCiachCiach
from rna_tools.tools.clarna_play.lib.rna_pdb_tools.pdb_parser_lib import StrucFile
from rna_tools import rna_tools_lib

# The Original settings of this import operation before adding
# MC-Annotate and rnaview

from rna_tools.tools.clarna_play.ClaRNAlib.utils import simplify_residue, compute_close_doublets, bench_start, bench_stop, load_json
# from utils import *
# ....there are quite a few functions involved in this module.

from rna_tools.tools.clarna_play.ClaRNAlib.distances import Residue
from rna_tools.tools.clarna_play.ClaRNAlib.clarna import MAX_RES_DIST, _find_contacts_using_library, ContactProgram, run_rnaview, parse_rnaview, run_mc_annotate, parse_mc_annotate, run_fr3d, parse_fr3d 
# NOTE: the programs run_rnaview, parse_rnaview, run_mc_annotate and
# parse_mc_annotate are actually from clarna_utils. The reason we can
# call it from here is because it was imported in clarna by the
# command "from clarna_utils import *", so this is why we don't have
# to specify clarna_utils here.


########################################
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
####### ClaRNA related imports  ########



class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__

#################

# Provides a definition of the command line arguments used in this program.
class Usage_Clarna:

  def __init__(self):
     # define a usage statement
     USAGE =  PROGRAM + ' -ipdb <infile>.pdb [-Clarna] [ -thresh f ] [ -bps ] [ -PS ] [ -stack ] [ -other ]\n'
     USAGE += PROGRAM + ' -ipdb <infile>.pdb -rnaview\n'
     USAGE += PROGRAM + ' -ipdb <infile>.pdb -mc_annotate\n'
     USAGE += PROGRAM + ' -ipdb <infile>.pdb -fr3d\n'
     USAGE +='\n'
     USAGE +='flag             status          function                             typical values\n'
     USAGE +='------         ---------       --------------------------------    -----------------\n'
     USAGE +=' -ipdb  f       required       input pdb file name f                 f   \'17re.pdb\' (required!)\n' 
     USAGE +=' -Clarna        (DEFAULT)      option to use Clarna itself               N/A\n' 
     USAGE +=' -thresh p     (optional)      minimum score allowed                 p   \'0.5\'\n' 
     USAGE +=' -bps          (optional)      analyze bps                               N/A\n'
     USAGE +=' -PS           (optional)      analyze base phosphate/sugar              N/A\n'
     USAGE +=' -stack        (optional)      analyze stacking                          N/A\n'
     USAGE +=' -other        (optional)      analyze other types of interactions       N/A\n'
     USAGE +=' -bp+stack     (optional)      analyzes both bps and stacking together   N/A\n'
     USAGE +='------         ---------       --------------------------------    -----------------\n'
     USAGE +=' -rnaview      required!       use rnaview, with clarna notation         N/A\n' 
     USAGE +='------         ---------       --------------------------------    -----------------\n'
     USAGE +=' -mc_annotate  required!       use MC-Annotate with clarna notation      N/A\n' 
     USAGE +='------         ---------       --------------------------------    -----------------\n'
     USAGE +=' -fr3d         required!       use FR3D with clarna notation             N/A\n' 
     USAGE +='------         ---------       --------------------------------    -----------------\n'
     USAGE +=' -h                -           prints this message                          \n' 
     USAGE +=' -help             -           ditto                                        \n' 
     USAGE +=' --help            -           ditto                                        \n' 
     
     self.u = USAGE
  #
  
  def get_Usage_Clarna(self):
     return 'Usage: %s\n' % self.u
  #
#

#################

#################

class CommandLine:
        
    def __init__(self):
        self.command_line = ''
        self.option_list  = ['-ipdb', '-iref', '-ichk', '-thresh', '-bps', '-stack', '-bp+stack', '-PS', '-other', '-Clarna', '-rnaview', '-mc_annotate']
        self.opt_cnt      = 0  # the option requested or the default (0)
        self.flnm_pdb     = '' # originally, I had "17ra_00_A.pdb" as a default
        self.chkpdb       = ''
        self.refpdb       = ''
        self.min_score    = 0.6
        self.clarna_opts    = "bps" 
        self.analPDB_opts = "Clarna"
        self.flag_for_SimRNA = False
    
    
    def parse_Clarna_CL(self, command_line):
        """parses the command line arguments for Entropy calls"""
        self.command_line = command_line
        debug_Clarna_CL_Opt = False
        n = len(self.command_line)
        if (n == 0):
            print('ERROR: command line undefined. ')
            sys.exit(1)
        flag_refPDB = False
        flag_chkPDB = False
        
        i = 1
        while i < n:
            if debug_Clarna_CL_Opt:
                print("index %d: \'%s\'" % (i, self.command_line[i]))
            # is it the input file?
            if self.command_line[i] == '-ipdb':
                #
                #  input file
                i = self.check_end(i, n) # increments i
                self.flnm_pdb = self.command_line[i]
           
                #
                if (self.flnm_pdb == ''):
                    emsg = "ERROR: input file name undefined."
                    self.error(emsg)
                #
                #
                if debug_Clarna_CL_Opt:
                    print("flnm_pdb name: %s" % self.flnm_pdb) 
                p = re.compile('.pdb$')
                a = p.findall(self.flnm_pdb)
                # print len(a)
                if (len(a) == 0):
                    emsg = "ERROR: input file (%s) must have the extension '.pdb'." % self.flnm_pdb
                    self.error(emsg)
                    #
                #
            #  how many
            elif self.command_line[i] == '-thresh':
                #
                i = self.check_end(i, n)
                
                try:
                    ms = self.command_line[i] 
                    self.min_score = float(ms)
                except:
                    print('ERROR: Must enter an float in requests using option \'-thresh\'.')
                    print('       You entered \'', ms, '\'.')
                    sys.exit(1)
                
                if self.min_score <= 0.0 or self.min_score >= 1.0: 
                    print('ERROR: \'-thresh\' only accepts 0.0 < thresh < 1.0.') 
                    print('       You entered \'%f\'.' % self.min_score)
                    sys.exit(1)
                
                if debug_Clarna_CL_Opt:
                    print("threshold to be displayed ", self.min_score)
                #
                #
            
            
            # option for analyzing base pair interactions
            # self.bps    = 'lib/classifier.bp.json.gz'
            elif self.command_line[i] == '-bps':
                self.clarna_opts = "bps"
                self.opt_cnt += 1
                if debug_Clarna_CL_Opt:
                    print("search base pairs requested.")
                #
            # option for stacking interactions
            # self.stack  = 'lib/classifier.stacking.json.gz'
            elif self.command_line[i] == '-stack':
                self.clarna_opts = "stack"
                self.opt_cnt += 1
                
                if debug_Clarna_CL_Opt:
                    print("search for stacking interactions requested.")
                #
            elif self.command_line[i] == '-bp+stack':
                self.clarna_opts = "bp+stack"
                self.opt_cnt += 1
                
                if debug_Clarna_CL_Opt:
                    print("search for stacking interactions requested.")
                #
            # option for analyzing phosphate/sugar base interactions (currently not *greatly* supported)
            # self.basePS = 'lib/classifier.base-phosphate.json.gz,lib/classifier.base-ribose.json.gz'
            elif self.command_line[i] == '-PS':
                self.clarna_opts = "PS"
                self.opt_cnt += 1
                
                if debug_Clarna_CL_Opt:
                    print("search base phosphate/sugar requested.")
                #
            # option for "other" types of interactions (currently not *greatly* supported)
            # self.other  = 'lib/classifier.other.json.gz,lib/classifier.other2.json.gz,lib/classifier.other3.json.gz'
            elif self.command_line[i] == '-other':
                self.clarna_opts = "other"
                
                self.opt_cnt += 1
                if debug_Clarna_CL_Opt:
                    print("search for \'all other\' interactions requested.")
                #
            # use Clarna (default)
            elif self.command_line[i] == '-Clarna':
                self.analPDB_opts = "Clarna"
                if debug_Clarna_CL_Opt:
                    print("requested to analyze using \'Clarna\'.")
                #
            # use RNAVIEW
            elif self.command_line[i] == '-rnaview':
                self.analPDB_opts = "rnaview"
                if debug_Clarna_CL_Opt:
                    print("requested to analyze using \'RNAVIEW\'.")
                #
            # use MC-Annotate
            elif self.command_line[i] == '-mc_annotate':
                self.analPDB_opts = "mc_annotate"
                if debug_Clarna_CL_Opt:
                    print("requested to analyze using \'MC-Annotate\'.")
                #
            # use FR3D
            elif self.command_line[i] == '-fr3d':
                self.analPDB_opts = "fr3d"
                if debug_Clarna_CL_Opt:
                    print("requested to analyze using \'FR3D\'.")
                #
            # is it a plea for help? (in one of its many different versions ...)
            elif (self.command_line[i] == '-h') or \
                 (self.command_line[i] == '--help') or \
                 (self.command_line[i] == '-help'):
                # default
                a = Usage_Clarna()
                print('%s\n' % a.get_Usage_Clarna())
                sys.exit(0)
                #
            #
            else:
                emsg = 'ERROR: command line is corrupted:\n'
                emsg += self.show_CL()
                emsg += '\noffending command \'%s\'' % self.command_line[i]
                self.error(emsg)
            #
            #
            i += 1
        #
        # error checking
        
        if self.opt_cnt > 1:
            emsg = 'ERROR: presently only one interaction option allowed\n' 
            for ol in self.option_list:
                emsg += '%s ' % ol
            emsg += '\n'
            emsg += self.show_CL()
            self.error(emsg)
        #
        if debug_Clarna_CL_Opt:
            print("INFORMATION: classification will be done using %s" % self.analPDB_opts)
            print("             all Clarna options will be ignored")
        
        print("Classifier: %s" % self.analPDB_opts)
        #
        
        # sys.exit(0)
        return 0
    #
    
    
    def check_end(self, i, n):
        if (i+1 < n):
            i += 1
        else:
            emsg = "ERROR, insufficient number of command line arguments.\n"
            self.error(emsg)
            # terminates
        return i
    
    
    def error(self, e_msg):
        """delivers error/exit messages for this class"""
        print('\n')
        print(e_msg)
        ss = 'input command: ' + self.show_CL() + '\n'
        print(ss)
        a = Usage_Clarna()
        print('%s' % a.get_Usage_Clarna())
        sys.exit(1)

        
    def show_CL(self):
        """separates the command line arguments into specific strings"""
        n = len(self.command_line)
        ss = ''
        for i in range(n):
            if (i < n-1):
                ss += '%s ' % self.command_line[i]
            else:
                ss += '%s' % self.command_line[i]
        return ss

class Clarna_utils:

    def __init__(self):
        # types of libraries available
        # main vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        # bp = base pairs   ! definitely needed
        # stacking          ! perhaps needed 
        # base-phosphate    ! not sure if needed or not
        # base-ribose       ! not sure if needed or not
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # other  = if the above is not enough
        # other2 = if the above is not enough
        # other3 = if even the above is not enough
        clarna_libs = ['lib/classifier.bp.json.gz', \
                       'lib/classifier.stacking.json.gz', \
                       'lib/classifier.base-phosphate.json.gz', \
                       'lib/classifier.base-ribose.json.gz', \
                       'lib/classifier.other.json.gz', \
                       'lib/classifier.other2.json.gz', \
                       'lib/classifier.other3.json.gz']
        
        self.init_opts = {'input_doublet_id': None, \
                          'pdb_id': None, \
                          'ignore_bad_doublets': None, \
                          'compare_with': None, \
                          'descriptions_dict': '../ClaRNA/descriptions-dict.json', \
                          'disable_align_to': True, \
                          'input': None, \
                          'show_scores_for': None, \
                          'use_old_params': False, \
                          'normalize_graph': False, \
                          'save_graph': None, \
                          'input_doublet_ids': None, \
                          'threads': 1, \
                          'lib': None, \
                          'disable_compare_distances': True, \
                          'data_dir': 'gc-data/pdb_files', \
                          'save_scores': None, \
                          'show_pair': None}
        
        
        self.bps    = 'lib/classifier.bp.json.gz'
        self.stack  = 'lib/classifier.stacking.json.gz'
        self.basePS = 'lib/classifier.base-phosphate.json.gz,lib/classifier.base-ribose.json.gz'
        self.other  = 'lib/classifier.other.json.gz,lib/classifier.other2.json.gz,lib/classifier.other3.json.gz'
        self.other1 = 'lib/classifier.other.json.gz'
        self.other2 = 'lib/classifier.other2.json.gz'
        self.other3 = 'lib/classifier.other3.json.gz'

        self.options = AttributeDict()
        
        # sets up the default library to lookup
        self.set_library("bps")

        # assign options for varies keys (not so helpful I know)
        for key,val in list(self.init_opts.items()):
            self.options[key] = val
        
        #bench_start("loading classifier libraries")
        self.libs = []
        for lib_fn in self.options['lib'].split(","):
            if not os.path.isfile(lib_fn):
                print("library %s does not exists" % lib_fn)
                exit(1)
            lib = ContactProgram()
            lib.load_from_json(lib_fn)
            self.libs.append(lib)
        
        #bench_stop()
        self.min_score    = 0.5
        self.clarna_opts    = "bps" 
    
    def set_min_score(self, min_score):
        self.min_score = min_score
    
    def set_clarna_opts(self, ch):
        self.clarna_opts = ch
    
    def set_library(self, libnm):
        if libnm == "bps":
            self.init_opts['lib'] = self.bps
            self.clarna_opts        = "bps" 
        elif libnm == "stack":
            self.init_opts['lib'] = self.stack
            self.clarna_opts        = "stack" 
        elif libnm == "bp+stack":
            self.init_opts['lib'] = self.bps + "," + self.stack
            # print self.init_opts['lib']
            # sys.exit(0)
            self.clarna_opts        = "bp+stack" 
        elif libnm == "PS":
            self.init_opts['lib'] = self.basePS
            self.clarna_opts        = "PS" 
        elif libnm == "other":
            self.init_opts['lib'] = self.other
            self.clarna_opts        = "other" 
        else: 
            print("don't recognize library (%s)" % libnm)
            a = Usage_Clarna()
            print('%s' % a.get_Usage_Clarna())
            sys.exit(1)
        #
        # reduce classifier database (standard baspairs and stacking only)
        # self.init_opts['lib'] = 'lib/classifier.bp.json.gz,lib/classifier.stacking.json.gz'  # !!!!!
        # reduce classifier database (standard baspairs only)
        # self.init_opts['lib'] = 'lib/classifier.bp.json.gz'  
        
        # update ClaRNA paths

        CLARNA_PATH = rna_tools_lib.get_rna_tools_path() + '/tools/clarna_play/ClaRNAlib/'
        self.init_opts['lib'] = ",".join([CLARNA_PATH+item for item in self.init_opts['lib'].split(',')])
        
        for key,val in list(self.init_opts.items()):
            self.options[key] = val
        
        #bench_start("loading classifier libraries")
        self.libs = []
        for lib_fn in self.options['lib'].split(","):
            if not os.path.isfile(lib_fn):
                print("library %s does not exists" % lib_fn)
                exit(1)
            lib = ContactProgram()
            lib.load_from_json(lib_fn)
            self.libs.append(lib)
        #
    #
    
    def load_pdb_fobject(self, fobject):
        parser = PDB.PDBParser(QUIET = True)
        res = parser.get_structure("c",fobject)
        for a in res.get_atoms():
            if re.match(r'^[A-Z]{1,2}[0-9]?\*$',a.id):
                a.id = a.id.replace("*","'")
        return res

    # -------------------------------------------------------------------------

    def extract_close_doublets_from_fobject(self, fpdbObj, ids=None, use_old_params=False):
        # TODO: this code should be rewritten, code duplications with StructureCiachCiach
        res_num = {}

        structure = self.load_pdb_fobject(fpdbObj)
        ciach = StructureCiachCiach(structure,dont_normalize=True)
        close_doublets = []
        s = []
        for model in structure:
            s_begin = len(s)
            residues = []
            for chain in model:
                for r in chain:
                    # XXX
                    r.resname = r.resname.upper()
                    rr = simplify_residue(r)
                    if rr is not None:
                        next_o3p = ciach._locate_backbone(r,"P")
                        if next_o3p is not None:
                            for a in next_o3p:
                                rr['NEXT:%s'%a.id] = a.get_coord().tolist() 
                        _id = chain.id+str(r.get_id()[1])
                        _num = r.get_id()[1]
                        _resname = r.get_resname().strip()
                        if ids is not None:
                            res_num[_id] = len(s) 
                        if not _resname in ['A','C','G','U']:
                            #print "ignoring %s (bad resname: %s)" % (_id,_resname)
                            residues.append(None)
                        else:
                            residues.append(r)
                        if not use_old_params:
                            s.append([chain.id, _resname, _num, Residue(_id,_resname,rr), _id])
                        else:
                            s.append([chain.id, _resname, _num, rr, _id])
                    else:
                        print("ignoring %s (missing atoms)" % _id)
            s_end = len(s)

            if ids is None:
                for (i,j) in compute_close_doublets(residues, max_distance=MAX_RES_DIST):
                    close_doublets.append((i,j,0))
                    close_doublets.append((j,i,0))

        if ids is not None:
            d_ids = [re.sub("^[^:]+:","",tmp) for tmp in ids]
            for id in d_ids:
                (id1,id2) = id.split(":")
                i = res_num.get(id1)
                j = res_num.get(id2)
                if i is None or j is None:
                    print("UNKNOWN doublet %s" % id)
                    continue
                close_doublets.append((i,j,0))
        return ("dummy_id", s, close_doublets)

    # runs clarna classification scripts
    def start_clarna(self, flnm):
        debug_start_clarna = False

        # open the file
        fpdbObj = open(flnm,"r")
        
        # fileObj is just the contents of the PDB file
        #bench_start("running classifier")
        (pdb_id, s, close_doublets) = self.extract_close_doublets_from_fobject(fpdbObj)
        
        res = _find_contacts_using_library(pdb_id, s, close_doublets, self.libs, self.options)
        #bench_stop()
        return res
    #
    
    # runs rnaview classification scripts
    def start_rnaview(self, flnm):
        debug_start_rnaview = False
        flag_show_warning = False

        # first read the PDB file and analyze it using rnaview
        rv_data = run_rnaview(flnm) 
        # rv_data corresponds to the standard output of rnaview with
        # default settings: you run it without any other options.
        
        # 
        rv_res = parse_rnaview(rv_data, stackings=True,use_chain=True)
        if debug_start_rnaview:
            print('rv_res:')
            print(rv_res)
        # rv_res is a list
            
        # [[('A1', 'A41'), ('G', 'C'), 'PP_cis_XIX'], 
        # [('A2', 'A40'), ('G', 'C'), 'PP_cis_XIX'], 
        #                    ....,
        # [('A35', 'A36'), ('A', 'A'), 'stacking'], 
        # [('A36', 'A37'), ('A', 'U'), 'stacking']]
        #
        
        # translate into something that can be compared with Clarna
        cb = ClarnaBabbel()
        
        if "RV" not in cb.descriptions: 
            print("ERROR: no definitions available for RNAview")
            print("       .... don't know what to do.")
            sys.exit(1)
        #
        
        rv_dict = {}
        for rv_res_k in rv_res:
            bb = rv_res_k[2]
            if bb in cb.descriptions["RV"]:
                # print cb.descriptions["RV"][bb]
                rv2c = cb.rnaview2clarna(rv_res_k)
                if debug_start_rnaview:
                    print(rv2c)
                #
                rv_dict.update(rv2c)
            else:
                if flag_show_warning:
                    print("INFORMATION: ignoring MC-Annotate key {%s:}..." % bb)
                #
                # print "trouble!!!, no key like \"%s\" exists!" % bb
                # sys.exit(1)
            #
        #
        rv_dict = self.rnaview_get_stacks(rv_res, rv_dict)
        return rv_dict
    #

    def rnaview_get_stacks(self, rv_res, bbclarna):
        debug_rnaview_get_stacks = False
        flag_show_warnings = False
        if debug_rnaview_get_stacks:
            print("stacking")
            print(rv_res)
        #
        
        # There is the VERY REAL PROBLEM that there is overlap here.
        # For example, consider 1a60_04_A.pdb, 1a60_13_A.pdb, or
        # 1anr_08_A.pdb
        
        # from 1a60_13_A.pdb:
        
        # Warning: found a key already in the dictionary
        # bb_5p:  {('A', 1, 'A', 2): {'bp': ('G', 'G'), '>>': 1.0}}
        # Warning: found a key already in the dictionary
        # bb_3p:  {('A', 16, 'A', 17): {'bp': ('C', 'C'), '>>': 1.0}}
        # Warning: found a key already in the dictionary
        # bb_5p:  {('A', 10, 'A', 11): {'bp': ('A', 'C'), '>>': 1.0}}
        # Warning: found a key already in the dictionary
        # bb_3p:  {('A', 30, 'A', 31): {'bp': ('G', 'G'), '>>': 1.0}}
        # Warning: found a key already in the dictionary
        # bb_3p:  {('A', 36, 'A', 37): {'bp': ('U', 'C'), '>>': 1.0}}
        
        if flag_show_warnings:
            # list all the keys already associated with stacking
            for bb in bbclarna:
                if ">>" in bbclarna[bb]:
                    print(bb, bbclarna[bb])
        
        # data format of rv_res: 
        
        # [[('A1', 'A41'), ('G', 'C'), 'PP_cis_XIX'], 
        # [('A2', 'A40'), ('G', 'C'), 'PP_cis_XIX'], 
        #                    ....,
        # [('A35', 'A36'), ('A', 'A'), 'stacking'], 
        # [('A36', 'A37'), ('A', 'U'), 'stacking']]
        
        for k in range(0, len(rv_res)-1):
            bbk1 = rv_res[k]
            bbk2 = rv_res[k+1]
            if debug_rnaview_get_stacks:
                print("bbk12: ", bbk1, bbk2)
            
            # RNAview does, fortunately, list the base-base
            # interactions sequentially in the keys.  Therefore, it is
            # not a big job to simply scan through the list of keys
            # and verify that there is actually a nearest neighbor
            # interaction there. If there is, then there is also
            # stacking, since the base-base interaction and a nearest
            # neighbor base-base clearly constitute the definition of
            # a stack.
            
            # base-base 1
            dt15     = bbk1[0][0]
            ch1_5p = dt15[0]
            rs1_5p = int(dt15[1:])
            nt1_5p = bbk1[1][0]
            dt13     = bbk1[0][1]
            ch1_3p = dt13[0]
            rs1_3p = int(dt13[1:])
            nt1_3p = bbk1[1][1]
            
            # base-base 2
            dt25     = bbk2[0][0]
            ch2_5p = dt25[0]
            rs2_5p = int(dt25[1:])
            nt2_5p = bbk2[1][0]
            dt23     = bbk2[0][1]
            ch2_3p = dt23[0]
            rs2_3p = int(dt23[1:])
            nt2_3p = bbk2[1][1]
            
            if debug_rnaview_get_stacks:
                print("dt15/13/25/23: ", dt15, dt13, dt25, dt23)
            
            # add stacking on the 5' end of the chain if the nearest
            # neighbor is one residue away.
            if debug_rnaview_get_stacks:
                print("rs1_5p, rs2_5p: ", rs1_5p, rs2_5p)
            if rs1_5p == rs2_5p - 1:
                key = (ch1_5p, rs1_5p, ch2_5p, rs2_5p)
                if key not in bbclarna:
                    bb_5p = { key : { "bp": (nt1_5p, nt2_5p), ">>": 1.0 } }
                    bbclarna.update(bb_5p)
                    if debug_rnaview_get_stacks:
                        print("bb_5p: ",  bb_5p)
                else:
                    if flag_show_warnings:
                        # Option to show this warning if the stacking
                        # keys already exists
                        print("Warning: found a key already in the dictionary")
                        bb_5p = { key : { "bp": (nt1_5p, nt2_5p), ">>": 1.0 } }
                        print("bb_5p: ",  bb_5p)
            #
            # add stacking on the 3' end of the chain if the nearest
            # neighbor is one residue away.
            if debug_rnaview_get_stacks:
                print("rs1_3p, rs2_3p: ", rs1_3p, rs2_3p)
            if rs1_3p - 1 == rs2_3p:
                key = (ch2_3p, rs2_3p, ch1_3p, rs1_3p)
                if key not in bbclarna:
                    bb_3p = { key : { "bp": (nt2_3p, nt1_3p), ">>": 1.0 } }
                    bbclarna.update(bb_3p)
                    if debug_rnaview_get_stacks:
                        print("bb_3p: ",  bb_3p)
                else:
                    if flag_show_warnings:
                        # Option to show this warning if the stacking
                        # keys already exists
                        print("Warning: found a key already in the dictionary")
                        bb_3p = { key : { "bp": (nt2_3p, nt1_3p), ">>": 1.0 } }
                        print("bb_3p: ",  bb_3p)
            #
        #
        return bbclarna
    #
    
    
    # runs MC-Annotate classification scripts
    def start_mc_annotate(self, flnm):
        debug_start_mc_annotate = False
        flag_show_warning = False

        # first read the PDB file and analyze it using mc_annotate
        mc_data = run_mc_annotate(flnm) 
        # mc_data corresponds to the standard output of mc_annotate which
        # you run it without any other settings.
        
        # print 'mc_data:'
        # print mc_data
        
        mc_res = parse_mc_annotate(mc_data, stackings=True,use_chain=True)
        # mc_res is a list of the above data: 
        
        # (
        #  [[('1', '41'), ('G', 'C'), 'WwWw_pairing_antiparallel_cis_XIX'], 
        #   [('2', '40'), ('G', 'C'), 'WwWw_pairing_antiparallel_cis_XIX'], ..... ],
        #  { {'24': {'resname': 'G', 'conf': 'anti'}, 
        #     '25': {'resname': 'A', 'conf': 'anti'}, ..... }  )
        # 
        
        if debug_start_mc_annotate:
            print('mc_res:', len(mc_res))
            for mc_i in mc_res:
                print("mcx: ", mc_i)
            #
        #
        
        # for MC-Annotate, there is additional information about the
        # particular arrangement of the base that seems to be not so
        # plainly mentioned in Clarna.  Therefore there are two data
        # sets generated when MC-Annotate is called. What can be
        # actually assessed in a given round with Clarna is the first
        # set of information.

        mc_bbres  = mc_res[0]
        mc_ntinfo = mc_res[1]

        # translate into something that can be compared with Clarna
        cb = ClarnaBabbel()
        
        if "MC" not in cb.descriptions: 
            print("ERROR: no definitions available for MC-Annotate")
            print("       .... don't know what to do.")
            sys.exit(1)
        #
        
        mc_dict = {}
        # first go through the base base interactions
        for mc_bbres_k in mc_bbres:
            bb = mc_bbres_k[2]
            if bb in cb.descriptions["MC"]:
                # print cb.descriptions["MC"][bb]
                mc2c = cb.mc_annotate2clarna(mc_bbres_k, mc_ntinfo)
                if debug_start_mc_annotate:
                    print(mc2c)
                #
                mc_dict.update(mc2c)
            #
            else:
                if flag_show_warning:
                    print("INFORMATION: ignoring MC-Annotate key {%s:}..." % bb)
                #
            #
        #
        return mc_dict
    #
    # runs MC-Annotate classification scripts
    def start_fr3d(self, flnm):
        debug_start_fr3d = False
        flag_show_warning = False

        # first read the PDB file and analyze it using fr3d
        fr_data = run_fr3d(flnm) 
        # fr_data corresponds to the standard output of fr3d which
        # you run it without any other settings.
        
        # print 'fr_data:'
        # print fr_data
        
        fr_res = parse_fr3d(fr_data, stackings=True,use_chain=True)
        # fr_res is a list of the above data: 
        
        # [[('A258', 'A257'), ('A', 'G'), 'stacking_s53'], 
        # [('A299', 'A257'), ('C', 'G'), 'cWW'], 
        # [('A257', 'A258'), ('G', 'A'), 'stacking_s35'], ....
        #              .....
        # [('A296', 'A295'), ('C', 'A'), 'n0BR'], 
        # [('A297', 'A296'), ('C', 'C'), 'n0BR'], 
        # [('A298', 'A297'), ('U', 'C'), 'n0BR'], 
        # [('A299', 'A298'), ('C', 'U'), 'n0BR']]
        
        if debug_start_fr3d:
            print('fr_res:', len(fr_res))
            for fr_i in fr_res:
                print("frx: ", fr_i)
            #
        #
        
        
        # translate into something that can be compared with Clarna
        cb = ClarnaBabbel()
        
        if "RV" not in cb.descriptions: 
            print("ERROR: no definitions available for FR3D")
            print("       .... don't know what to do.")
            sys.exit(1)
        #
        
        fr_dict = {}
        for fr_res_k in fr_res:
            bb = fr_res_k[2]
            if bb in cb.descriptions["FR"]:
                if re.match('^[HSW][HSW]_', cb.descriptions["FR"][bb]) \
                   or re.match('^[<>][<>]$', cb.descriptions["FR"][bb]):
                    # print cb.descriptions["FR"][bb]
                    fr2c = cb.fr3d2clarna(fr_res_k)
                    if debug_start_fr3d:
                        print(fr2c)
                    #
                    fr_dict.update(fr2c)
                    #
                else:
                    if flag_show_warning:
                        print("INFORMATION: skipping unrelated FR3D key {%s:}..." % bb)
                    #
                #
            else:
                if flag_show_warning:
                    print("INFORMATION: ignoring untranslatable FR3D key {%s:}..." % bb)
                #
                # print "trouble!!!, no key like \"%s\" exists!" % bb
                # sys.exit(1)
        #
        return fr_dict
    #

class ClarnaBabbel:
    """ClaRNABabbel"""
    def __init__(self):
        # the file "descriptions-dict.json" contains all the critical
        # definitions that allow translation between various other
        # classifier methods into clarna notation: basically rnaview,
        # MC-Annotate, and FR3D
        self.descriptions = load_json(os.path.join(CLARNA_PATH,"descriptions-dict.json"))
        #
    #
    def rnaview2clarna(self, bbk):
        # bbk = base-base information for item k 
        
        # data format after running parse_rnaview
        
        # [[('A1', 'A41'), ('G', 'C'), 'PP_cis_XIX'], 
        # [('A2', 'A40'), ('G', 'C'), 'PP_cis_XIX'], 
        #                    ....,
        # [('A35', 'A36'), ('A', 'A'), 'stacking'], 
        # [('A36', 'A37'), ('A', 'U'), 'stacking']]
        
        # In the current information format, rnaview does not provide
        # any weight for the reliability score.  Therefore, since
        # clarna requires it, we have to make it up.
        
        nt  = bbk[0][0]
        ch1 = nt[0]
        rs1 = nt[1:]
        nt =bbk[0][1]
        ch2 = nt[0]
        rs2 = nt[1:]
        
        bbk2c = self.descriptions["RV"][bbk[2]]
        # add the required stacking information 
        if bbk[2] == "stacking":
            bbk2c = ">>"
        #
        bbclarna = {(ch1, int(rs1), ch2, int(rs2)) : \
                    { "bp": (bbk[1][0], bbk[1][1]), bbk2c: 1.0 }  }
        return bbclarna
    #
    
    
    
    def mc_annotate2clarna(self, bbk, ntinfo):
        # bbk = base-base information for item k 
        
        # data format after running parse_mc_annotate
        
        # (
        #  [[('A1', 'A41'), ('G', 'C'), 'WwWw_pairing_antiparallel_cis_XIX'], 
        #   [('A2', 'A40'), ('G', 'C'), 'WwWw_pairing_antiparallel_cis_XIX'], ..... ],
        #  { {'A24': {'resname': 'G', 'conf': 'anti'}, 
        #     'A25': {'resname': 'A', 'conf': 'anti'}, ..... }  )
        # 
        
        
        nm1 = nm2 = ''
        nt1  = bbk[0][0]
        ch1 = nt1[0]
        rs1 = nt1[1:]
        nt2 =bbk[0][1]
        ch2 = nt2[0]
        rs2 = nt2[1:]
        
        if bbk[1][0] == '-' or bbk[1][1] == '-':
            nm1 = ntinfo[nt1]['resname']
            nm2 = ntinfo[nt2]['resname']
        else:
            nm1 = bbk[1][0]            
            nm2 = bbk[1][1]            
        #
        bbk2c = self.descriptions["MC"][bbk[2]]
        bbclarna = {(ch1, int(rs1), ch2, int(rs2)) : \
                    { "bp": (nm1, nm2), bbk2c: 1.0 }  }
        
        return bbclarna
    #
    def fr3d2clarna(self, bbk):
        # bbk = base-base information for item k 
        
        # data format after running parse_rnaview
        
        # [[('A258', 'A257'), ('A', 'G'), 'stacking_s53'], 
        # [('A299', 'A257'), ('C', 'G'), 'cWW'], 
        # [('A257', 'A258'), ('G', 'A'), 'stacking_s35'], ....
        #              .....
        # [('A296', 'A295'), ('C', 'A'), 'n0BR'], 
        # [('A297', 'A296'), ('C', 'C'), 'n0BR'], 
        # [('A298', 'A297'), ('U', 'C'), 'n0BR'], 
        # [('A299', 'A298'), ('C', 'U'), 'n0BR']]
        
        # In the current information format, FR3D does not provide
        # any weight for the reliability score.  Therefore, since
        # clarna requires it, we have to make it up.
        
        # print "bbk: ", bbk
        nt  = bbk[0][0]
        ch1 = nt[0]
        rs1 = nt[1:]
        nt =bbk[0][1]
        ch2 = nt[0]
        rs2 = nt[1:]
        
        bbk2c = self.descriptions["FR"][bbk[2]]
        
        bbclarna = {(ch1, int(rs1), ch2, int(rs2)) : \
                    { "bp": (bbk[1][0], bbk[1][1]), bbk2c: 1.0 }  }
        return bbclarna

class SeeClarna:
    """SeeClarna"""
    def __init__(self):
        self.pdbclasslist={} # I am not sure if this MUST be here
        self.cutils = ''
        self.min_score    = 0.5
        self.clarna_opts  = "bps"
    #
    
    # key operations
    
    # when there is only one structure key recorded
    def get_RNAstruct_key(self,pdbdt):
        zz = list(pdbdt.keys())
        key = ''
        for zzi in zz:
            if (not zzi == "bp"):
                key = zzi
            #
        #
        return key
    #

    # when there are more than two structure keys recorded
    def get_other_keys(self,pdbdt):
        zz = list(pdbdt.keys())
        keys = []
        for zzi in zz:
            if (not zzi == "bp"):
                keys += [zzi]
            #
        #
        return keys
    #
    
    
    
    def make_structure_dict(self, u, v, d, pdbdict):
        debug_make_structure_dict = False
        # print u, v
        xc = u[0]
        yc = v[0]
        xr = int(u[1:])
        yr = int(v[1:])
        xn = d['n_type'][0]
        yn = d['n_type'][1]
        t = d['desc']
        w = d['weight']
        
        uc = ''; ux = 0.0; un = ''
        vc = ''; vx = 0.0; vn = ''
        sgn = 1
        if not (self.clarna_opts == "PS" or self.clarna_opts == "other"):
            # don't reorder if it is PS or other
            # check the order and switch if vc < uc
            if yr < xr:
                # must switch the order
                uc = yc; ux = yr; un = yn
                vc = xc; vx = xr; vn = xn
                sgn = -1
            else:
                uc = xc; ux = xr; un = xn
                vc = yc; vx = yr; vn = yn
                sgn = 1
            #
        else: 
            uc = xc; ux = xr; un = xn
            vc = yc; vx = yr; vn = yn
            sgn = 1
        #
        
        tag = "++"
        tt = t
        # should we change the orientation of the base-base interaction
        if sgn < 0:
            if self.cutils.clarna_opts == 'bps' or self.cutils.clarna_opts == 'stack' or self.cutils.clarna_opts == 'bp+stack':
                if t[0] == '<' or t[0] == '>':
                    if t == '<<' or  t == '>>':
                        if t == '<<':
                            tt = '>>'
                        else:
                            tt = '<<'
                    else:
                        tt = t
                    #
                else:
                    tt=t[1] + t[0] + t[2:] # we have to reverse the type
                #
                t = tt
            #
            tag = "--"
            #
            
        # Does (ux,vx) already exist
        if (uc,ux,vc,vx) in pdbdict:
            # Does pdbdict[(ux,vx)] already have configuration tt?
            if t in pdbdict[(uc,ux,vc,vx)]:
                pdbdict[(uc,ux,vc,vx)][t] += w
                pdbdict[(uc,ux,vc,vx)][t] *= 0.5
            else:
                pdbdict[(uc,ux,vc,vx)].update({t : w})
        else:
            # no key exists for (uc, ux, vc, vx)
            pdbdict[(uc, ux, vc, vx)] = {'bp': (un,vn), t:w}
        #
        key = self.get_RNAstruct_key(pdbdict[(uc, ux, vc, vx)])
        if debug_make_structure_dict:
            # The main purpose here is to manage bps and stacking, but
            # this is able to store results of options -PS and
            # -other. I do not recommend doing analysis with PS or
            # other, but for the sake of printing something, I allow
            # this. Presently, the real analysis and comparisons are
            # restricted to bps and stacking.
            print("%s (%3d,%3d), bb %s-%s, %15s, (%8.5f -> %8.5f)" % \
                (tag, ux,vx, \
                 pdbdict[(uc, ux, vc, vx)]['bp'][0], \
                 pdbdict[(uc, ux, vc, vx)]['bp'][1], \
                 t, w, pdbdict[(uc, ux, vc, vx)][t]))
                
        #else:
        #print "%s %s:%s(%3d) -- %s:%s(%3d)    %20s   %8.5f (%s)"\
        #        % (tag, uc, un, ux, vc, vn, vx, t, w, self.cutils.clarna_opts)
        #
        return pdbdict
    #
    
    
    def eval_PDB(self, cl):
        """eval_PDB"""
        debug_eval_PDB_with_clarna = False
        flnm = cl.flnm_pdb 
        min_score = cl.min_score
        analPDB_opts = cl.analPDB_opts        
        self.min_score    = cl.min_score
        self.clarna_opts  = cl.clarna_opts

        
        self.cutils = Clarna_utils()
        self.cutils.set_min_score(cl.min_score)
        self.cutils.set_clarna_opts(cl.clarna_opts)
        self.cutils.set_library(cl.clarna_opts)
        
        self.pdbclasslist={}
        
        if analPDB_opts == "Clarna":
            res_graph = self.cutils.start_clarna(flnm)
            # res_graph does not produce any displayable output.
        
            for (u,v,d) in res_graph.edges(data=True):
                # print u, v, d
                # continue
                if d['type']=='contact' and d['weight'] > min_score: 
                    self.pdbclasslist = self.make_structure_dict(u,v,d, self.pdbclasslist)
                # print u,v, d['n_type'], d['desc'], d['weight']
            #
        #
        elif analPDB_opts == "rnaview":
            # calculate using RNAVIEW
            self.pdbclasslist = self.cutils.start_rnaview(flnm)
        elif analPDB_opts == "mc_annotate":
            # calculate using MC-Annotate
            self.pdbclasslist = self.cutils.start_mc_annotate(flnm)
        elif analPDB_opts == "fr3d":
            # calculate using FR3D
            self.pdbclasslist = self.cutils.start_fr3d(flnm)
        else:
            print("ERROR: selected option (%s) not available" % analPDB_opts)
            sys.exit(1)
        #
        print(self.display_bbstack(self.pdbclasslist))
    #

    #
    def bubbleSort(self, key_order):
        debug_bubbleSort = False
        if debug_bubbleSort:
            print('length: ', len(key_order))
            print('input:  ', key_order)
        #
        for passnum in range(len(key_order)-1, 0, -1):
            for i in range(passnum):
                if key_order[i][1]>key_order[i+1][1]:
                    temp = key_order[i]
                    key_order[i] = key_order[i+1]
                    key_order[i+1] = temp
                #
            #
        #
        if debug_bubbleSort:
            print('result: ', key_order)
        return key_order
    #

    def display_bbstack(self, pdbanal):
        debug_display_bbstack = False
        bb_list = self.bubbleSort(list(pdbanal.keys()))
        if debug_display_bbstack:
            print(pdbanal)
        result = ''
        for bb_k in bb_list:
            strng = '%s %4d   %s %4d' % (bb_k[0], bb_k[1], bb_k[2], bb_k[3])
            strng += '          bp %s %s' % (pdbanal[bb_k]["bp"][0], pdbanal[bb_k]["bp"][1])
            if len(list(pdbanal[bb_k].keys())) == 2:
                key = self.get_RNAstruct_key(pdbanal[bb_k])
                #
                strng += '    %20s %8.4f' % (key, pdbanal[bb_k][key])
                if debug_display_bbstack:
                    print(key)
                    print(strng)
                    # print bb_k, pdbanal[bb_k], pdbanal[bb_k].keys()
            else:
                keys = self.get_other_keys(pdbanal[bb_k])
                for key in keys:
                    strng += '    %20s %8.4f' % (key, pdbanal[bb_k][key])
                if debug_display_bbstack:
                    print("MORE THAN ONE!")
                    print(strng)
                    # print bb_k, pdbanal[bb_k], pdbanal[bb_k].keys()
            #
            result += strng + '\n'
        #
        return result
    #

class testClarna:
    """testClarna"""
    def __init__(self):
        self.cutils = ''
    #
    
    def test_clarna(pdbid="1ehz", min_score=0.5):

        url = urllib.request.urlopen(url="http://www.rcsb.org/pdb/files/%s.pdb"%pdbid.upper())

        self.cutils = Clarna_utils()
        res_graph = self.cutils.start(url)
        for (u,v,d) in res_graph.edges(data=True):
            if d['type']=='contact' and d['weight']>min_score: print(u,v, d['n_type'], d['desc'], d['weight'])
        #
    #
#main
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("ERROR: %s requires at least one argument!" % PROGRAM)
        a = Usage_Clarna()
        print('%s' % a.get_Usage_Clarna())
        sys.exit(1)
    
    cl = CommandLine()
    cl.parse_Clarna_CL(sys.argv)
    print('chains: ', StrucFile(cl.flnm_pdb).get_info_chains())
    seeClarna = SeeClarna()

    seeClarna.eval_PDB(cl)
