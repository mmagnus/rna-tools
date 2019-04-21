#!/usr/bin/env python

"""
Example::

    ./rna_add_chain.py -c X ../../input/1msy_rnakbmd_decoy999_clx_noChain.pdb \
    > ../../output/1msy_rnakbmd_decoy999_clx_noChain_Xchain.pdb

From::

    ATOM      1  O5'   U     1      42.778  25.208  46.287  1.00  0.00
    ATOM      2  C5'   U     1      42.780  26.630  45.876  1.00  0.00
    ATOM      3  C4'   U     1      42.080  27.526  46.956  1.00  0.00
    ATOM      4  O4'   U     1      43.013  28.044  47.963  1.00  0.00
    ATOM      5  C1'   U     1      42.706  29.395  48.257  1.00  0.00
    ATOM      6  N1    U     1      43.857  30.305  47.703  1.00  0.00
    ATOM      7  C6    U     1      45.057  29.857  47.308  1.00  0.00
    ATOM      8  C5    U     1      46.025  30.676  46.763  1.00  0.00
    ATOM      9  C4    U     1      45.720  32.110  46.702  1.00  0.00
    ATOM     10  O4    U     1      46.444  32.975  46.256  1.00  0.00

to::

    ATOM      1  O5'   U X   1      42.778  25.208  46.287  1.00  0.00
    ATOM      2  C5'   U X   1      42.780  26.630  45.876  1.00  0.00
    ATOM      3  C4'   U X   1      42.080  27.526  46.956  1.00  0.00
    ATOM      4  O4'   U X   1      43.013  28.044  47.963  1.00  0.00
    ATOM      5  C1'   U X   1      42.706  29.395  48.257  1.00  0.00
    ATOM      6  N1    U X   1      43.857  30.305  47.703  1.00  0.00
    ATOM      7  C6    U X   1      45.057  29.857  47.308  1.00  0.00
    ATOM      8  C5    U X   1      46.025  30.676  46.763  1.00  0.00
    ATOM      9  C4    U X   1      45.720  32.110  46.702  1.00  0.00
    ATOM     10  O4    U X   1      46.444  32.975  46.256  1.00  0.00"""
from rna_tools.rna_tools_lib import *
import argparse

def get_parser():
    parser =  argparse.ArgumentParser()#usage="%prog [<options>] <pdb files (test_data/*)>")
    parser.add_argument('file', help="file")
    parser.add_argument('-c',"--chain",
                         dest="chain",
                        help='a new chain, e.g. A')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print((parser.print_help()))
        sys.exit(1)

    s = RNAStructure(args.file)
    for l in s.lines:
        if l.startswith('ATOM'):
            nl = list(l)
            nl[21] =  args.chain # new chain
            nl = ''.join(nl)
            print(nl)
        else:
            print(l)
