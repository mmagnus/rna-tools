#!/usr/bin/env python

"""used for reconstructing the dictonary from a run_clarna.sh output
file of structure built with the flag -ipdb"""

import sys
import re
import math

# ###############  GLOBAL DEFINTIONS  ################
PROGRAM = "clarna_compare.py"  # program name
ANAL_EXT='.pdb.outCR$'      # the default extension is for CLARNA files
# ####################################################

# Provides a definition of the command line arguments used in this program.
class Usage_AnalClarna:
    """Usage_AnalClarna"""
    def __init__(self, v_extension):
        # define a usage statement
        ac_ext = v_extension[:len(v_extension)-1] # remove the '$' character that is used for class re
    
        USAGE  = PROGRAM + ' -iref <refpdb>%s -ichk <SimRNApdb>%s [ options ]\n' % (ac_ext, ac_ext)
        USAGE +='\n'
        USAGE +='flag             status        function                                 typical values\n'
        USAGE +=' -iref  fr      required     input reference clarna file               \'17ra%s\'\n'  % (ac_ext) 
        USAGE +=' -ichk  fc      required     input SimRNA result clarna file           \'17ra_clust%s\'\n'  % (ac_ext) 
        USAGE +=' -----------  ----------   ------  options  ------------------------   ------  \n'
        USAGE +=' option unique to ClaRNA:    (ext \"outCR\")                                   \n'
        USAGE +=' -thresh p         -         accept only values where threshold > p         0.6\n'
        USAGE +=' -----------  ----------   -----------------------------------------   ------  \n'
        USAGE +=' Display ClaRNA __WITH__ other classifiers:                                    \n'
        USAGE +=' -w_rnaview        -         include RNAVIEW results                           \n'
        USAGE +=' -w_mc_annotate    -         include MC-Annotate results                       \n'
        USAGE +=' -w_fr3d           -         include FR3D results                              \n'
        USAGE +=' -w_all_methods    -         include both RNAVIEW and MC-Annotate              \n'
        USAGE +=' ###########       #         ## NOTE: input extension is \"outCR\"             \n'
        USAGE +='                                but \"outRV\", and/or \"outMC\"                \n'
        USAGE +='                                 and/or \"outFR\" _must_ exist!!!              \n'
        USAGE +=' -----------  ----------   -----------------------------------------   ------  \n'
        USAGE +=' Display other classifiers (without ClaRNA):                                   \n'
        USAGE +=' -fr3d             -         analyze FR3D results (ext \"outFR\")              \n'
        USAGE +=' -mc_annotate      -         analyze MC-annotate results (ext \"outRV\")       \n'
        USAGE +=' -rnaview          -         analyze RNAVIEW results (ext \"outMC\")           \n'
        USAGE +=' -----------  ----------   ------  other options  ------------------   ------  \n'
        USAGE +=' Printing options:                                                             \n'
        USAGE +=' -CSV              -         print summary as Comma Separated Values           \n'
        USAGE +=' -details,-v       -         show all categories of TP, FP and FN              \n'
        USAGE +=' -verbose,-v       -         ditto                                             \n'
        USAGE +=' -----------  ----------   --------  help  -------------------------   ------  \n'
        USAGE +=' -h/-help/--help   -         prints this message                               \n' 
        
        self.u = USAGE

    def get_Usage_AnalClarna(self):
        return 'Usage: %s\n' % self.u

class CommandLine:

    """CommandLine"""
    def __init__(self):
        self.command_line = ''
        self.option_list  = ['-iref', '-ichk', '-w_fr3d', '-w_rnaview', '-w_mc_annotate', '-w_all_methods', '-rnaview', '-mc_annotate', '-fr3d', '-CSV', '-v', '-verbose', '-details']
        self.chkout       = ''
        self.refout       = ''
        self.classifiers  = 'Clarna' # default
        self.v_ext        = ANAL_EXT
        self.min_score    = 0.6
        self.details      = False
        self.CSV          = False

    def parse_AnalClarna_CL(self, command_line):
        """parses the command line arguments for Entropy calls"""

        self.command_line = command_line
        debug_AnalClarna_CL_Opt = False
        n = len(self.command_line)
        if (n == 0):
            print('ERROR: command line undefined. ')
            sys.exit(1)
        #
        flag_refPDB = False
        flag_chkPDB = False
        
        # scan through the command line arguments for classifier types
        i = 1
        self.v_ext = ANAL_EXT # default extension is Clarna
        
        # read a ClaRNAfied file???
        while i < n:
          if self.command_line[i] == '-rnaview':
            self.v_ext = '.pdb.outRV$'        # analyze only RNAVIEW results (ext: pdb.outRV)
            self.classifier = "rnaview"
          elif self.command_line[i] == '-mc_annotate':
            self.v_ext = '.pdb.outMC$'        # analyze only MC-Annotate results (ext: pdb.outMC)
            self.classifier = "mc_annotate"
          if self.command_line[i] == '-fr3d':
            self.v_ext = '.pdb.outFR$'        # analyze only RNAVIEW results (ext: pdb.outRV)
            self.classifier = "fr3d"
          i += 1
        #
        ac_ext = self.v_ext[:len(self.v_ext)-1] # remove the '$' character that is used for class re

        
        i = 1
        while i < n:
            if debug_AnalClarna_CL_Opt:
                print("index %d: \'%s\'" % (i, self.command_line[i]))
            # is it the input file?
            #  SimRNA: enter reference PDB structure
            elif self.command_line[i] == '-iref':
                flag_refPDB = True
                
                #  input file
                i = self.check_end(i, n)
                self.refout = self.command_line[i]
                
                #
                if (self.refout == ''):
                    emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    emsg += "ERROR: input file name undefined.\n"
                    emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    self.error(emsg)
                #
                #
                if debug_AnalClarna_CL_Opt:
                    print("reference PDB bb structure: %s" % self.refout) 
                p = re.compile(self.v_ext)
                a = p.findall(self.refout)
                # print len(a)
                if (len(a) == 0):
                    emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    emsg += "ERROR: input file (%s) must have the extension \'*%s\'.\n" % (self.refout, ac_ext)
                    emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    self.error(emsg)
                    #
                #
            #  SimRNA: enter SimRNA calculated structure (or cluster) in PDB format
            elif self.command_line[i] == '-ichk':
                flag_chkPDB = True
                #  input file
                i = self.check_end(i, n)
                self.chkout = self.command_line[i]
                
                if (self.chkout == ''):
                    emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    emsg += "ERROR: input file name undefined.\n"
                    emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    self.error(emsg)
                #
                #
                if debug_AnalClarna_CL_Opt:
                    print("SimRNA PDB bb structure:    %s" % self.chkout) 
                p = re.compile(self.v_ext)
                a = p.findall(self.chkout)
                # print len(a)
                if (len(a) == 0):
                    emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    emsg += "ERROR: input file (%s) must have the extension \'*%s\'.\n" \
                    %   (self.chkout, ac_ext)
                    emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                    self.error(emsg)
            elif self.command_line[i] == '-thresh':
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
                
                if debug_AnalClarna_CL_Opt:
                    print("threshold to be displayed ", self.min_score)
            elif (self.command_line[i] == '-details') or \
                 (self.command_line[i] == '-verbose') or \
                 (self.command_line[i] == '-v'):
                self.details = True
                if self.CSV:
                    print("WARNING: \'%s\' option not supported with option \'-CSV\'" % self.command_line[i])
                    self.CSV = False
            elif (self.command_line[i] == '-CSV'):
                self.CSV = True
                if self.details:
                    print("WARNING: \'CSV\' option not supported with option \'-details,-v,-verbose\'")
                    self.CSV = False
            # include results from RNAVIEW classifier
            elif (self.command_line[i] == '-w_rnaview'):
                self.classifiers = "w_rnaview"
            # include results from MC-Annotate classifier
            elif (self.command_line[i] == '-w_mc_annotate'):
                self.classifiers = "w_mc_annotate"
            # include results from FR3D classifier
            elif (self.command_line[i] == '-w_fr3d'):
                self.classifiers = "w_fr3d"
            # include results from RNAVIEW, MC-Annotate and FR3D classifiers
            elif (self.command_line[i] == '-w_all_methods'):
                self.classifiers = "w_all_methods"
            # analyze RNAVIEW results:
            elif (self.command_line[i] == '-rnaview'):
                self.classifiers = "rnaview"
            # analyze MC-Annotate results:
            elif (self.command_line[i] == '-mc_annotate'):
                self.classifiers = "mc_annotate"
            # analyze FR3D classifier
            elif (self.command_line[i] == '-fr3d'):
                self.classifiers = "fr3d"
            # is it a plea for help? (in one of its many different versions ...)
            elif (self.command_line[i] == '-h') or \
                 (self.command_line[i] == '--help') or \
                 (self.command_line[i] == '-help'):
                # default
                a = Usage_AnalClarna(self.v_ext)
                print('%s\n' % a.get_Usage_AnalClarna())
                sys.exit(0)
            else:
                emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                emsg += 'ERROR: command line is corrupted:\n'
                emsg += '\noffending command \'%s\'\n' % self.command_line[i]
                emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                self.error(emsg)
            i += 1
        #
        # error checking
        if (not flag_refPDB) or (not flag_chkPDB):
            if (not flag_refPDB):
                emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                emsg += 'ERROR: REFERENCE pdb.out file not specified in SimRNA analysis\n'
                emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                emsg += '\n'
                self.error(emsg)
            else:
                emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                emsg += 'ERROR: SimRNA RESULT pdb.out file not specified in SimRNA analysis\n'
                emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
                emsg += '\n'
                self.error(emsg)

        if debug_AnalClarna_CL_Opt:
            print('input reference PDB file: *.pdb = %s' % self.refout)
            print('input SimRNA PDB file:    *.pdb = %s' % self.chkout)

        return 0
    
    def check_end(self, i, n):
        if (i+1 < n):
            i += 1
        else:
            emsg  = '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
            emsg += "ERROR, insufficient number of command line arguments.\n"
            emsg += '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
            self.error(emsg)
            # terminates
        return i

    # delivers error/exit messages for this class
    def error(self, e_msg):
        print('\n')
        print(e_msg)
        ss = 'input command: ' + self.show_CL() + '\n'
        print(ss)
        a = Usage_AnalClarna(self.v_ext)
        print('%s' % a.get_Usage_AnalClarna())
        sys.exit(1)

    # separates the command line arguments into specific strings
    def show_CL(self):
        n = len(self.command_line)
        ss = ''
        for i in range(n):
            if (i < n-1):
                ss += '%s ' % self.command_line[i]
            else:
                ss += '%s' % self.command_line[i]
        return ss

class Data:
    """Data"""
    def __init__(self):
        self.flnm = ''
        self.pdbdt = {}
        self.min_score = 0.6
        self.data_type = "Clarna"
        self.verbose = False
    
    def readData(self, flnm):
        """readData, read ClaRNA output

        for debugging you can use `debug_readData = False` inside this function.
        
        :params flnm: file name"""
        debug_readData = False
        self.flnm = flnm
        try: 
            fldt = open(self.flnm, 'r')
        except IOError:
            print("ERROR: Cannot open file %s" % flnm)
            sys.exit(1)
        dlist = fldt.readlines()
        nlines = len(dlist)
        self.pdbdt = {}
        title = dlist[0].split()
        classifier = ''
        if len(title) > 0:
            if not (title[0][0:8] == "Classifi") or not (len(title) == 2):
                print("ERROR: Cannot understand the type of classifier")
                print("       valid files should contain \"Classifier: <type>\"")
                print("       in the first line of the file")
                sys.exit(1)
            else:
                classifier = title[1]
        if debug_readData:
            print("classifier: ", classifier)
        for i in range(1,nlines):
            # print len(dlist[i])
            if dlist[i].startswith('chains'): # skip line starting with chains
                continue
            if len(dlist[i]) > 25:
                # print dlist[i]
                v = self.reconstructDict(dlist[i])
                if not (v == {}):
                    self.pdbdt.update(v)
        if debug_readData:
            print('pdbdt:          ', self.pdbdt)
            print('pdbdt.keys():   ', list(self.pdbdt.keys()))
            print('Data: verbose = ', self.verbose)
        if self.verbose:
            pd = PrintData()
            print(pd.display_bbstack(self.pdbdt))
        return self.pdbdt
    
    def reconstructDict(self, strng):
        """150511wkd: I had some problems early on here where there was
        some trouble constructing v; however, it seems to have been
        fixed.  At least when I tried to do a full scan, I didn't
        encounter any issues any longer.  

        I think the reason was that RNAVIEW (at the time) did not
        override errors in "Clarnafied" types.  Most types are
        identified by Clarna, but there seem to be some rare
        varieties in MC-annotate and RNAVIEW (and probably FR3D)
        that are either not defined properly in Clarna or are bugs
        in these existing programs.  When I ask for verbose with
        run_clarna.sh, these errors are also expressed, but
        originally, I assumed that RNAVIEW and clarna had a full
        listing.  

        At any rate, if there are assignment errors generated by the
        python program, I would first guess that they are most
        likely to be occurring here.

        (this was between v and ^)
        
        :param: strng - a line of ClaRNA output
        """
        
        debug_reconstructDict = False
        # strng = "A    1   A   21          bp G C                  WW_cis   0.8637                  WW_cis   0.8637"
        # ss = strng.split()
        # ss
        # ['A', '1', 'A', '21', 'bp', 'G', 'C', 'WW_cis', '0.8637', 'WW_cis', '0.8637']
        
        sv = []
        v = {}
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        sv = strng.split()
        if len(sv) >= 7:
            v = {(sv[0], int(sv[1]), sv[2], int(sv[3])) : { sv[4]: (sv[5], sv[6])} }
        else:
            print("problems with reconstruction (%s):" % self.flnm) 
            print("sv: ", sv)
            sys.exit(1)
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Nstructs = (len(sv) - 7)/2
        key = list(v.keys())
        Nused = 0 # normally, these should equal Nstruct after the for loop.
        if debug_reconstructDict:
            print('new key: ', key)
            print('nstrcts: ', Nstructs)

        for i in range(0, int(Nstructs)):
            if (float(sv[8 + 2*i]) >= self.min_score):
                # if does not meet the minimum score, should be rejected
                struct = { sv[7 + 2*i] : float(sv[8 + 2*i]) }
                v[key[0]].update(struct)
                Nused += 1
                if debug_reconstructDict:
                    print('struct:  ', struct)
            else:
                # print user information if requested
                if self.verbose:
                    print("data below threshold (rejected sv): ", sv)
        if Nused == 0:
          v = {} # all entries were rejected, so v is canceled
        return v

class PrintData:
    """PrintData"""
    def __init__(self):
        self.debug = False
    # print out base-base stack
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
                # print 'pdbanal: ', pdbanal[bb_k]
                xx = self.get_other_keys(pdbanal[bb_k])
                key = xx[0]
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
            result += strng + '\n'
        return result

    def bubbleSort(self, key_order):
        """bubbleSort: I took the easy way....  This sorts the order of the
        keys so the output is a little more ordered than what is stored
        in the dictionary. The mixture of stacking and base pairing data
        is a bit annoying.  It looks a little more readable if you can
        separate the base pairing from the stacking.  At this point, I
        am not sure I want to waste time on the appearance of the
        output."""
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
        if debug_bubbleSort:
            print('result: ', key_order)
        return key_order
    
    # key operations
    # when there are more than two structure keys recorded
    def get_other_keys(self,bb_struct):
        zz = list(bb_struct.keys())
        keys = []
        for zzi in zz:
            if (not zzi == "bp"):
                keys += [zzi]
        return keys

class Analyze:
    """Analyze"""
    def __init__(self):
        """
        .stat_details = False # prints out all TP, FP, FN information
        .CSV = False # self.stat_details shown with CSV
        """
        self.refflnm = ''
        self.refdata = {}
        self.rWCbpKeys = []
        self.rnWCbpKeys = []
        self.rstackKeys = []
        self.refflnm = ''
        self.chkdata = {}
        self.cWCbpKeys = []
        self.cnWCbpKeys = []
        self.cstackKeys = []

        self.verbose = True
        self.stat_details = False 
        self.CSV          = False
        self.flag_pass_setup = False
        self.print_with_file_list = True  # output display
        
        self.pr = PrintData() # for scanning critical keys
        
        self.min_score    = 0.6
        # TP and FN related information (reference vs prediction)
        self.r_vs_p_WC    = [] # WC
        self.r_vs_p_nWC   = [] # nWC
        self.r_vs_p_stack = [] # stacking
        self.nref         = 0 # reference structure
        self.nchk_match   = 0 # predicted pairs that match the reference
        self.nFN          = 0 # nref - nchk_match (should be zero if perfect)
        # FP related terms information (prediction vs reference)
        self.p_vs_r_WC    = [] # WC
        self.p_vs_r_nWC   = [] # nWC
        self.p_vs_r_stack = [] # stacking
        self.nchk         = 0 # predicted structure
        self.nref_match   = 0 # reference pairs that match the prediction
        self.nFP          = 0 # nchk - nref_match (should be zero if perfect)
    
    def setup_Analysis(self, refout, chkout, cl):
        self.refflnm      = refout
        self.chkflnm      = chkout
        self.min_score    = cl.min_score
        self.verbose      = cl.details # for Data files
        self.stat_details = cl.details # for internal files
        self.CSV          = cl.CSV
    
    def set_verbose(self, b):
        self.verbose = b

    def set_print_with_file_list(self, b):
        self.print_with_file_list = b
    
    def get_outdata(self):
        # Reference input file (*.pdb.out)
        if self.verbose:
            print(self.refflnm)
        # transfer settings from CommandLine
        rdt = Data()
        rdt.verbose   = self.verbose
        rdt.min_score = self.min_score
        if self.verbose:
            print(self.refflnm)
        
        self.refdata = rdt.readData(self.refflnm)        
        self.rWCbpKeys  = self.get_bbkeys(self.refdata, "WCbp")
        self.rnWCbpKeys = self.get_bbkeys(self.refdata, "nWCbp")
        self.rstackKeys = self.get_bbkeys(self.refdata, "stack")
        # SimRNA RESULT file (*.pdb.out)
        if self.verbose:
            print(self.chkflnm)
        cdt = Data()
        cdt.verbose =   self.verbose
        cdt.min_score = self.min_score
        self.chkdata = cdt.readData(self.chkflnm)        
        self.cWCbpKeys  = self.get_bbkeys(self.chkdata, "WCbp")
        self.cnWCbpKeys = self.get_bbkeys(self.chkdata, "nWCbp")
        self.cstackKeys = self.get_bbkeys(self.chkdata, "stack")
        # print rdt.refdata 
        # print cdt.chkdata
        self.flag_pass_setup = True
    
    def get_bbkeys(self, ddata, bbtype):
        debug_get_bbkeys = False
        if bbtype == "stack":
            if debug_get_bbkeys:
                print("stacking keys")
        elif bbtype == "WCbp": 
            if debug_get_bbkeys:
                print("only canonical Watson-Crick base pairs")
        elif bbtype == "nWCbp":
            if debug_get_bbkeys:
                print("non-Watson-Crick base pairs")

        keys = list(ddata.keys())
        bbkeys = []
        self.pr = PrintData()
        for k in range(0,len(keys)):
            kstr_k = list(ddata[keys[k]].keys())
            xx =  self.pr.get_other_keys(ddata[keys[k]])
            kstr = xx[0]
            if debug_get_bbkeys:
                print('kstr = ', kstr)
            if bbtype == "stack":
                if kstr[0] == '>' or kstr[0] == '<':
                    bbkeys += [keys[k]]
            elif bbtype == "WCbp" or bbtype == "nWCbp":
                if not kstr[0] == '>' and not kstr[0] == '<': # not a stack
                    bp = ddata[keys[k]]['bp'][0] + ddata[keys[k]]['bp'][1]
                    canonical = (bp == "GC" or bp == "CG" or bp == "AU" or bp == "UA")
                    if debug_get_bbkeys:
                        print('bp: ', bp)
                    if bbtype == "WCbp":
                        if canonical and (kstr == "WW_cis"):
                            bbkeys += [keys[k]]
                            # print "flag_WC = True, bbtype = ", bbtype 
                    else:
                        # non WC pairing interactions
                        if not canonical or not (kstr == "WW_cis"):
                            bbkeys += [keys[k]]
                            # print "flag_WC = False, bbtype = ", bbtype 
            if debug_get_bbkeys:
                print(kstr)
                print(kstr_k)
        if debug_get_bbkeys:
            print('%s keys:    ' % (bbtype), bbkeys)
            for bpk in bbkeys:
                print(list(ddata[bpk].keys()))
        return bbkeys
    
    def analyze_data(self):
        """Analyze data"""
        if not self.flag_pass_setup:
            print("ERROR: input data is not properly define")
            print("       need to run get_outdata or get_pdbdata.")
            sys.exit(1)
        debug_analyze_data = False
        if debug_analyze_data:
            print("organized listing of structures found")
            print('refdata  WCbp: ')  
            for k in range(0, len(self.rWCbpKeys)):
                print(self.refdata[self.rWCbpKeys[k]])
            print('refdata nWCbp: ')  
            for k in range(0, len(self.rnWCbpKeys)):
                print(self.refdata[self.rnWCbpKeys[k]])
            print('refdata stack: ')  # , self.refdata.keys()
            for k in range(0, len(self.rstackKeys)):
                print(self.refdata[self.rstackKeys[k]])
            
            print('chkdata  WCbp: ')  # , self.chkdata.keys()
            for k in range(0, len(self.cWCbpKeys)):
                print(self.chkdata[self.cWCbpKeys[k]])
            print('chkdata nWCbp: ')  # , self.chkdata.keys()
            for k in range(0, len(self.cnWCbpKeys)):
                print(self.chkdata[self.cnWCbpKeys[k]])
            print('chkdata stack: ')  # , self.chkdata.keys()
            for k in range(0, len(self.cstackKeys)):
                print(self.chkdata[self.cstackKeys[k]])
            #
        #
        
        # match reference canonical WC base pairs
        self.r_vs_p_WC    = self.test_against_reference_motifs(self.rWCbpKeys, "WC")
        # match reference non-canonical WC base pairs
        self.r_vs_p_nWC   = self.test_against_reference_motifs(self.rnWCbpKeys, "nWC")
        # match reference stacking
        self.r_vs_p_stack = self.test_against_reference_motifs(self.rstackKeys, "stack")
        
        # reference structure
        self.nref       = self.r_vs_p_WC[0] + self.r_vs_p_nWC[0] + self.r_vs_p_stack[0]
        self.nchk_match = self.r_vs_p_WC[1] + self.r_vs_p_nWC[1] + self.r_vs_p_stack[1]
        self.nFN        = self.nref - self.nchk_match

        # match predicted canonical WC base pairs
        self.p_vs_r_WC    = self.test_against_predicted_motifs(self.cWCbpKeys, "WC")
        # match predicted non-canonical WC base pairs
        self.p_vs_r_nWC   = self.test_against_predicted_motifs(self.cnWCbpKeys, "nWC")
        # match predicted stacking
        self.p_vs_r_stack = self.test_against_predicted_motifs(self.cstackKeys, "stack")
        
        # predicted structure
        self.nchk       = self.p_vs_r_WC[0] + self.p_vs_r_nWC[0] + self.p_vs_r_stack[0]
        self.nref_match = self.p_vs_r_WC[1] + self.p_vs_r_nWC[1] + self.p_vs_r_stack[1]
        self.nFP        = self.nchk - self.nref_match
        rslt = self.summarize_analysis()
        return rslt
    
    def test_against_reference_motifs(self, refkeys, mtype):
        """count how many instances of the given types of motifs of the
        reference structure (refkeys) are matched in the predicted
        structure. This should correspond to the number of true
        positives (motifs matching the reference) and false
        negatives (reference motifs missing in the predition).
        
        Hence, 
         
        TP = n_match_ref/n_ref_motifs
        FN = (n_ref_motifs - n_match_ref)/n_ref_motifs
        
        and 
        
        1 = TP + FN"""
        
        debug_this = False
        n_match_ref = 0
        n_ref_motifs = len(refkeys)
        for k in range(0, n_ref_motifs):
            if refkeys[k] in self.chkdata:
                refmotif = self.pr.get_other_keys(self.refdata[refkeys[k]])
                chkmotif = self.pr.get_other_keys(self.chkdata[refkeys[k]])
                if len(refmotif) > 1:
                    print("WARNING: %s has more than one values %s: " % (mtype, self.refflnm), refmotif)
                if len(chkmotif) > 1:
                    print("WARNING: %s has more than one values %s: " % (mtype, self.chkflnm), chkmotif)
                if debug_this:
                    print(refmotif, " <-> ", chkmotif)
                n_match_ref += 1
        if debug_this:
            print("%5s    ref: %5d   pred match: %5d" % (mtype, n_ref_motifs, n_match_ref))
        return [n_ref_motifs, n_match_ref]
    
    def test_against_predicted_motifs(self, chkkeys, mtype):
        """Count how many instances of the given types of motifs of the
        predicted structure (chkkeys) are matched in the calculated
        structure. The DIFFERENCE between the counted motifs of the
        predicted structure and the matching motifs of the reference
        structure should correspond to the number of false positive
        motifs (FP).
        
        Hence,
        
        FP = (n_chk_motifs - n_match_chk)/n_chk_motifs"""
        
        debug_this = False
        n_match_chk = 0
        n_chk_motifs = len(chkkeys)
        for k in range(0, n_chk_motifs):
            if chkkeys[k] in self.refdata:
                refmotif = self.pr.get_other_keys(self.refdata[chkkeys[k]])
                chkmotif = self.pr.get_other_keys(self.chkdata[chkkeys[k]])
                if len(refmotif) > 1:
                    print("WARNING: %s has more than one values %s: " % (mtype, self.refflnm), refmotif)
                if len(chkmotif) > 1:
                    print("WARNING: %s has more than one values %s: " % (mtype, self.chkflnm), chkmotif)
                if debug_this:
                    print(refmotif, " <-> ", chkmotif)
                n_match_chk += 1
            #
        #
        if debug_this:
            print("%5s   pred: %5d    ref match: %5d" % (mtype, n_chk_motifs, n_match_chk))
        return [n_chk_motifs, n_match_chk]

    def calc_inf(self, TP, FP, FN, verbose=False):
        """Interaction Network Fidelity (essentially Matthew's
        correlation coefficient when TN --> infinity)
        
        -- magnus::
        
         if x == 0 and y == 0:
            inf = 0"""
        inf = 0 # -999.999
        x = TP + FP
        y = TP + FN
        if verbose: print('TP:', TP, 'FP:', FP, 'FN:', FN)
        if x > 0.0 and y > 0.0:
            inf = math.sqrt(TP**2/((TP + FP)*(TP + FN)))
        if x == 0 and y == 0:
            inf = 0
        return inf
    
    def summarize_analysis(self): 
        """
        The definitions are taken from https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
        
        -- magnus change::

         if TP_WC == 0 and FP_WC == 0.0:
             PPV_WC = 0
        
        true positive (TP)
             eqv. with hit
        true negative (TN)
             eqv. with correct rejection
        false positive (FP)
             eqv. with false alarm, Type I error
        false negative (FN)
             eqv. with miss, Type II error         
        
        P and N:
         P = TP + FN
         N = FP + TN
        
        sensitivity or true positive rate (TPR)
            eqv. with hit rate, recall
            TPR = TP / P = TP / (TP+FN)

        specificity (SPC) or True Negative Rate
            SPC = TN / N = TN / (FP + TN) 
             
        precision or positive predictive value (PPV)
            PPV = TP / (TP + FP)
        negative predictive value (NPV)
             NPV = TN / (TN + FN)
        fall-out or false positive rate (FPR)
             FPR = FP / N = FP / (FP + TN)
        false discovery rate (FDR)
             FDR = FP / (FP + TP) = 1 - PPV 
        Miss Rate or False Negative Rate (FNR)
             FNR = FN / (FN + TP) 
        
        accuracy (ACC)
             ACC = (TP + TN) / (P + N)

        F1 score
             is the harmonic mean of precision and sensitivity
             F1 = 2 TP / (2 TP + FP + FN)

        Matthews correlation coefficient (MCC)
             (TP*TN - FP*FN)/[(TP+FP)(TP+FN)(TN+FP)(TN+FN)]^(1/2)  
        
        Informedness
            TPR + SPC - 1
        Markedness
            PPV + NPV - 1 
        
        #####################################################
        
        Unfortunately, converting this information into something
        that is easy to understand is still a bit complicated!!!
        Let's work with an example to illustrate how the terms are
        used::
        
           motifs         motifs
          reference      prediction
             A               
             B               B
             C               C
             D               
             E               E
                             F
                             G
                             H
        
        where we assume that the reference motifs are the "truth
        incarnate" (this is not necessarily so, but we assume it)::
        
         nTP = 3: B,C,E
         nFP = 3: F,G,H
         nFN = 2: A,D
        
        so this much is easy.
        
        Now, recall that "ref" is the reference structure and "chk"
        is the predicted structure (to be checked)::
        
         nref       = r_vs_p_WC[0] + r_vs_p_nWC[0] + r_vs_p_stack[0]
         nchk_match = r_vs_p_WC[1] + r_vs_p_nWC[1] + r_vs_p_stack[1]
         nFN        = nref - nchk_match
         nTP        = nchk_match
        
        clearly ``nref = nFP + nTP, so ppv is nTP/nref``::
        
         nchk       = p_vs_r_WC[0] + p_vs_r_nWC[0] + p_vs_r_stack[0]
         nref_match = p_vs_r_WC[1] + p_vs_r_nWC[1] + p_vs_r_stack[1]
         nFP        = nchk - nref_match
        
        so for individual cases such as WC pairing::
        
         nTP = r_vs_p_WC[1]
         nFN = r_vs_p_WC[0] - r_vs_p_WC[1]
         nFP = p_vs_r_WC[0] - p_vs_r_WC[1]
        """
        debug_summarize_analysis = False
        
        rf = self.refflnm.split("/")
        cf = self.chkflnm.split("/")
        
        p = re.compile('/')
        mcf = p.findall(self.chkflnm)
        if debug_summarize_analysis:
            print('# matchings of ref file with / path', len(mcf))
        mrf = p.findall(self.refflnm)
        if debug_summarize_analysis:
            print('# matchings of ref file with / path', len(mrf))
        refflnm = self.refflnm
        if (len(mcf) > 0):
            rf = self.refflnm.split("/")
            refflnm = rf[len(rf)-1]
        chkflnm = self.chkflnm
        if (len(mcf) > 0):
            cf = self.chkflnm.split("/")
            chkflnm = cf[len(cf)-1]
        #
        
        sdt = '' # string data
        if  self.stat_details:
            sdt += "compare:\n"
            sdt +=  "%s     %s\n\n" % (refflnm, chkflnm)
            # TP and FN related information (reference vs prediction)
            sdt +=  "reference vs prediction:\n"
            sdt +=  "%5s    ref %5d   pred-match %5d\n" % ("WC",    self.r_vs_p_WC[0],    self.r_vs_p_WC[1])
            sdt +=  "%5s    ref %5d   pred-match %5d\n" % ("nWC",   self.r_vs_p_nWC[0],   self.r_vs_p_nWC[1])
            sdt +=  "%5s    ref %5d   pred-match %5d\n" % ("stack", self.r_vs_p_stack[0], self.r_vs_p_stack[1])
            sdt +=  "\n"
            # FP related terms information (prediction vs reference)
            sdt +=  "prediction vs reference:"
            sdt +=  "%5s   pred %5d    ref-match %5d\n" % ("WC",    self.p_vs_r_WC[0],    self.p_vs_r_WC[1])
            sdt +=  "%5s   pred %5d    ref-match %5d\n" % ("nWC",   self.p_vs_r_nWC[0],   self.p_vs_r_nWC[1])
            sdt +=  "%5s   pred %5d    ref-match %5d\n" % ("stack", self.p_vs_r_stack[0], self.p_vs_r_stack[1])
            sdt +=  "\n"
            sdt +=  "for all motifs:"
            sdt +=  "nref %4d,  nchkTP %4d,  nFN %4d\n" % (self.nref, self.nchk_match, self.nFN)
            sdt +=  "nchk %4d,  nrefTP %4d,  nFP %4d\n" % (self.nchk, self.nref_match, self.nFP)
            sdt +=  "\n\n"
        #
        
        # TP = n_match_ref
        TP_all     = float(self.nchk_match)
        TP_WC      = float(self.r_vs_p_WC[1])
        TP_nWC     = float(self.r_vs_p_nWC[1])
        TP_stack   = float(self.r_vs_p_stack[1])
        
        # FN      number_of_reference_motifs  - number_of_predicted_motifs
        FN_all     = float(self.nref - self.nchk_match)
        if debug_summarize_analysis:
            print("compare: ", self.nref - self.nchk_match, self.nFN)
        FN_WC      = float(self.r_vs_p_WC[0]    - self.r_vs_p_WC[1])
        FN_nWC     = float(self.r_vs_p_nWC[0]   - self.r_vs_p_nWC[1])
        FN_stack   = float(self.r_vs_p_stack[0] - self.r_vs_p_stack[1])
        
        # FP = (n_chk_motifs - n_match_chk)/n_chk_motifs
        FP_all     = float(self.nchk - self.nref_match)
        if debug_summarize_analysis:
            print("compare: ", self.nchk - self.nref_match, self.nFP)
        FP_WC      = float(self.p_vs_r_WC[0]    - self.p_vs_r_WC[1])
        FP_nWC     = float(self.p_vs_r_nWC[0]   - self.p_vs_r_nWC[1])
        FP_stack   = float(self.p_vs_r_stack[0] - self.p_vs_r_stack[1])
        
        inf_all   = self.calc_inf(TP_all,   FP_all,   FN_all)
        inf_stack = self.calc_inf(TP_stack, FP_stack, FN_stack)
        inf_WC    = self.calc_inf(TP_WC,    FP_WC,    FN_WC)
        inf_nWC   = self.calc_inf(TP_nWC,   FP_nWC,   FN_nWC)
        
        # SNS_all = TP_all/float(self.nchk) # still not sure about this one
        # PPV_all = TP_all/(TP_all + FP_all)
        #

        SNS_WC = 0 #-999.999
        if self.r_vs_p_WC[0] > 0:
            SNS_WC = TP_WC/float(self.r_vs_p_WC[0])

        PPV_WC = 0 # -999.999
        if (TP_WC + FP_WC) > 0.0:
            PPV_WC = TP_WC/(TP_WC + FP_WC)
        # --magnus 
        if TP_WC == 0 and FP_WC == 0:
            PPV_WC = 0
        # --end
        SNS_nWC = 0 #-999.999
        if self.r_vs_p_nWC[0] > 0:
            SNS_nWC = TP_nWC/float(self.r_vs_p_nWC[0])
        PPV_nWC = 0 # -999.999
        if (TP_nWC + FP_nWC) > 0.0:
            PPV_nWC = TP_nWC/(TP_nWC + FP_nWC)
        # --magnus 
        if TP_nWC == 0 and FP_nWC == 0:
            PPV_nWC = 0
        # --end
        if  self.stat_details:
            sdt +=  "summary:\n"
            sdt +=  "inf_all   %8.3f\n" % inf_all
            sdt +=  "inf_stack %8.3f\n" % inf_stack
            sdt +=  "inf_WC    %8.3f\n" % inf_WC
            sdt +=  "inf_nWC   %8.3f\n" % inf_nWC
            sdt +=  "SNS_WC    %8.3f\n" % SNS_WC
            sdt +=  "PPV_WC    %8.3f\n" % PPV_WC
            if SNS_nWC < 0.0:
              sdt +=  "SNS_nWC     NA\n" 
            else:
              sdt +=  "SNS_nWC   %8.3f\n" % SNS_nWC
            if PPV_nWC < 0.0:
              sdt +=  "PPV_nWC     NA\n" 
            else:
              sdt +=  "PPV_nWC   %8.3f\n" % PPV_nWC
              
        else:
            if self.CSV:
                # prints a single line
                if self.print_with_file_list:
                    sdt += "%s, %s," % (refflnm, chkflnm)
                # print inf_nWC if it exists. [Note that there is no
                # case where inf_nWC < 0 and both (SNS_nWC > 0 and
                # PPV_nWC > 0)
                if inf_all < 0.0:
                    sdt += " NA" 
                else:
                    sdt += "%8.3f" % inf_all
                #
                if inf_stack < 0.0:
                    sdt += ", NA" 
                else:
                    sdt += ", %8.3f" % inf_stack
                #
                if inf_WC < 0.0:
                    sdt += ", NA" 
                else:
                    sdt += ", %8.3f" % inf_WC
                #
                if inf_nWC < 0.0:
                    sdt += ", NA"
                else:
                    sdt += ", %8.3f" % inf_nWC
                #
                # print SNS_nWC and/or PPV_nWC if it exists
                if SNS_WC < 0.0:
                    sdt += ", NA"
                else:
                    sdt += ", %8.3f" % SNS_WC
                #
                if PPV_WC < 0.0:
                    sdt += ", NA"
                else:
                    sdt += ", %8.3f" % PPV_WC
                #
                if SNS_nWC < 0.0:
                    sdt += ", NA"
                else:
                    sdt += ", %8.3f" % SNS_nWC
                #
                if PPV_nWC < 0.0:
                    sdt += ", NA"
                else:
                    sdt += ", %8.3f" % PPV_nWC
                #
            #
            else:
                # prints a single line
                if self.print_with_file_list:
                    sdt += "%s     %40s" % (refflnm, chkflnm)
                # print inf_nWC if it exists. [Note that there is no
                # case where inf_nWC < 0 and both (SNS_nWC > 0 and
                # PPV_nWC > 0)

                if inf_all < 0.0:
                    sdt += "      NA   " 
                else:
                    sdt += "   %8.3f" % inf_all

                if inf_stack < 0.0:
                    sdt += "      NA   " 
                else:
                    sdt += "   %8.3f" % inf_stack

                if inf_WC < 0.0:
                    sdt += "      NA   " 
                else:
                    sdt += "   %8.3f" % inf_WC

                if inf_nWC < 0.0:
                    sdt += "      NA   "
                else:
                    sdt += "   %8.3f" % inf_nWC

                # print SNS_nWC and/or PPV_nWC if it exists
                if SNS_WC < 0.0:
                    sdt += "      NA   "
                else:
                    sdt += "   %8.3f" % SNS_WC

                if PPV_WC < 0.0:
                    sdt += "      NA   "
                else:
                    sdt += "   %8.3f" % PPV_WC

                if SNS_nWC < 0.0:
                    sdt += "      NA   "
                else:
                    sdt += "   %8.3f" % SNS_nWC

                if PPV_nWC < 0.0:
                    sdt += "      NA   "
                else:
                    sdt += "   %8.3f" % PPV_nWC

        return sdt

# main
if __name__ == '__main__':
    cl = CommandLine()
    cl.parse_AnalClarna_CL(sys.argv)
    s_clarna = ''
    if not ((cl.classifiers == "mc_annotate") or (cl.classifiers == "rnaview")):
        # either only ClaRNA or with other applications
        refoutCR = cl.refout
        chkoutCR = cl.chkout
        anal_clarna = Analyze()
        # default settings for output
        anal_clarna.set_verbose(False)
        anal_clarna.set_print_with_file_list(True)
        # -----
        anal_clarna.setup_Analysis(refoutCR, chkoutCR, cl)
        anal_clarna.get_outdata()
        s_clarna +=  anal_clarna.analyze_data()
        if cl.classifiers == "w_mc_annotate" or cl.classifiers == "w_all_methods": 
            refoutMC = cl.refout[:-2] + "MC"
            chkoutMC = cl.chkout[:-2] + "MC"
            # print "refoutMC,chkoutMC: ", refoutMC, chkoutMC
            anal_mc_annotate = Analyze()
            # default settings for output
            anal_mc_annotate.set_verbose(False)
            anal_mc_annotate.set_print_with_file_list(False)
            # -----

            anal_mc_annotate.setup_Analysis(refoutMC, chkoutMC, cl)
            anal_mc_annotate.get_outdata()
            s_clarna += anal_mc_annotate.analyze_data()
        if cl.classifiers == "w_rnaview" or cl.classifiers == "w_all_methods": 
            refoutRV = cl.refout[:-2] + "RV"
            chkoutRV = cl.chkout[:-2] + "RV"
            # print "refoutRV,chkoutRV: ", refoutRV, chkoutRV
            anal_rnaview = Analyze()
            # default settings for output
            anal_rnaview.set_verbose(False)
            anal_rnaview.set_print_with_file_list(False)
            # -----
            anal_rnaview.setup_Analysis(refoutRV, chkoutRV, cl)
            anal_rnaview.get_outdata()
            s_clarna +=  anal_rnaview.analyze_data()
        if cl.classifiers == "w_fr3d" or cl.classifiers == "w_all_methods": 
            refoutFR = cl.refout[:-2] + "FR"
            chkoutFR = cl.chkout[:-2] + "FR"
            # print "refoutFR,chkoutFR: ", refoutFR, chkoutFR
            anal_fr3d = Analyze()
            # default settings for output
            anal_fr3d.set_verbose(False)
            anal_fr3d.set_print_with_file_list(False)
            # -----
            anal_fr3d.setup_Analysis(refoutFR, chkoutFR, cl)
            anal_fr3d.get_outdata()
            s_clarna +=  anal_fr3d.analyze_data()
    else:
        if cl.classifiers == "mc_annotate":
            refoutMC = cl.refout
            chkoutMC = cl.chkout
            anal_mc_annotate = Analyze()
            # default settings for output
            anal_mc_annotate.set_verbose(False)
            anal_mc_annotate.set_print_with_file_list(True)
            # -----
            anal_mc_annotate.setup_Analysis(refoutMC, chkoutMC, cl)
            anal_mc_annotate.get_outdata()
            s_clarna += anal_mc_annotate.analyze_data()
        elif cl.classifiers == "rnaview":
            refoutRV = cl.refout
            chkoutRV = cl.chkout
            anal_rnaview = Analyze()
            # default settings for output
            anal_rnaview.set_verbose(False)
            anal_rnaview.set_print_with_file_list(True)
            # -----
            anal_rnaview.setup_Analysis(refoutRV, chkoutRV, cl)
            anal_rnaview.get_outdata()
            s_clarna +=  anal_rnaview.analyze_data()
        elif cl.classifiers == "rnaview":
            refoutFR = cl.refout
            chkoutFR = cl.chkout
            anal_rnaview = Analyze()
            # default settings for output
            anal_fr3d.set_verbose(False)
            anal_fr3d.set_print_with_file_list(True)
            # -----
            anal_fr3d.setup_Analysis(refoutFR, chkoutFR, cl)
            anal_fr3d.get_outdata()
            s_clarna +=  anal_fr3d.analyze_data()
        else:
            print("no information!")

    print(s_clarna)
