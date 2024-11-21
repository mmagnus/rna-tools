#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import requests


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("family", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.family):
        args.family = [args.family]

    for family in args.family:
            
            # URL of the Rfam seed alignment
            url = f"https://rfam.org/family/{family}/alignment/fastau"
            
            # Output file name
            output_file = f"{family}.seed.sto"

            # Download the file
            response = requests.get(url)
            if response.status_code == 200:
                with open(output_file, "wb") as f:
                    f.write(response.content)
                print(f"Downloaded seed alignment saved as {output_file}")
            else:
                print(f"Failed to download. Status code: {response.status_code}")
