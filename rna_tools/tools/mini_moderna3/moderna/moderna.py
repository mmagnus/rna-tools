#!/usr/bin/env python
#
# moderna_cli.py
#
# Command-line interface
# 
# http://iimcb.genesilico.pl/moderna/
#
__author__ = "Magdalena Rother, Tomasz Puton, Kristian Rother"
__copyright__ = "Copyright 2008, The Moderna Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "1.7.1"
__maintainer__ = "Magdalena Rother"
__email__ = "mmusiel@genesilico.pl"
__status__ = "Production"

"""
 ModeRNA

 A program for comparative RNA structure modeling

 Authors: 
    Magdalena Rother
    Kristian Rother
    Tomasz Puton
    Janusz Bujnicki

 (c) 2008
 
 www.genesilico.pl 

"""

import sys, os
from rna_tools.tools.mini_moderna3.moderna.Constants import MODULE_PATH

#KR: manipulating the import order in sys.path 
# to avoid a conflict between moderna/ and moderna.py
# It's ugly but tested.
# DON'T TRY THIS AT HOME.
while MODULE_PATH in sys.path:
    sys.path.remove(MODULE_PATH)
sys.path.append(MODULE_PATH)

from rna_tools.tools.mini_moderna3.moderna import *
from optparse import OptionParser
from commands import *

########################################  MAIN  ##########################################

def prepare_parser():
    """
    ModeRNA (C) 2008 by Magdalena Musielak, Kristian Rother et al.
    
    Prepares simple script base interface.
    """
    usage = """usage: %prog [options] arg
    
Examples:
python moderna.py -t my_template.pdb -c my_template_chain -a my_alignment.fasta -o model_output.pdb
python moderna.py -s my_template.pdb -c my_template_chain -m m22G -p 3 -o modification_output.pdb"""
    parser = OptionParser(usage, version="%prog "+__version__)

    parser.add_option("-a", "--alignment", dest="a_file", help="Reads alignment from A_FILE")
    parser.add_option("-t", "--template", dest="t_file", help="Reads a template structure from T_FILE (.pdb or .ent).")
    parser.add_option("-e", "--examine", action="store_true", dest="examine_template", help="Checks whether template structure may cause any problems during modeling")
    parser.add_option("-l", "--clean", action="store_true", dest="clean_template", help="Cleans up template structure")
    parser.add_option("-s", "--structure", dest="s_file", help="Reads a single chain from the structure in S_FILE (.pdb 0r .ent).")
    parser.add_option("-c", "--chain", dest="chain_name", default="A", help="Chain you want to read (-s or -t) (default is A)")
    parser.add_option("-o", "--output", dest="o_file", default="ModeRNA.out.pdb", help="Writes the output model structure to O_FILE")
    parser.add_option("-m", "--modification", dest="modification", help="Modification name to be inserted (abbreviation like m2G, mnm5U, Y)")
    parser.add_option("-p", "--position", dest="position", help="Modification position (residue number)")
    #parser.add_option("-l" "--logfile", dest="l_file", default="ModeRNA.out.log", help="Writes log message to L_FILE")

    (options, args) = parser.parse_args()
    return parser, options


def parser_usage(parser, options):
    """
    Defines the script behaviour after calling program with some options.
    """
    model = None

    if options.t_file:
        try:
            t = load_template(options.t_file, options.chain_name)
            if options.examine_template: 
                print(examine_structure(t))
            if options.clean_template: 
                clean_structure(t)
            print('TEMPLATE SEQUENCE:\n%s'%t.get_sequence())
           
        except KeyError:
            parser.error("Bad chain name. Add your template chain name as -c option")

    if options.a_file and options.t_file:
        a = load_alignment(options.a_file)
        model = create_model(t,a)
     
    if options.s_file:
        try:
            model=load_model(options.s_file, options.chain_name)
        except KeyError:
            parser.error("Bad chain name. Add your template chain name as -c option")

    if options.position and options.modification and model:
        add_modification(model[options.position], options.modification, model)

    if options.a_file and not options.t_file:
        parser.error("Both alignment and template are required to create model")

    if (options.position or options.modification) and ( not (options.position and options.modification) or not model):
        parser.error("Modification name, modification position (residue number) and structure required for adding modification") 
    
    if options.a_file and options.t_file and options.s_file:
        parser.error("Too many options")
           
    if model: 
        write_model(model,options.o_file)


if __name__ == "__main__":
    p,o=prepare_parser()
    parser_usage(p,o)

