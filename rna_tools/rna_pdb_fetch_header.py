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


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument('-', "--", help="", default="")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser

import requests

def fetch_pdb_header(pdb_id):
    """
    Fetches the header information for a given PDB ID using the RCSB PDB GraphQL API.

    :param pdb_id: The PDB ID to query.
    :return: A dictionary containing the header information or an error message.
    """
    # Define the GraphQL query
    query = """
    query getPDBHeader($id: String!) {
      entry(entry_id: $id) {
        rcsb_id
        struct {
          title
        }
        exptl {
          method
        }
        rcsb_accession_info {
          deposit_date
          initial_release_date
        }
        citation {
          title
          rcsb_authors
          year
          pdbx_database_id_DOI
        }
      }
    }
    """
    
    # Define the variables
    variables = {"id": pdb_id}
    
    # API endpoint
    url = "https://data.rcsb.org/graphql"
    headers = {"Content-Type": "application/json"}
    
    try:
        # Send the POST request
        response = requests.post(url, json={"query": query, "variables": variables}, headers=headers)
        response.raise_for_status()  # Raise HTTPError for bad responses (4xx and 5xx)
        data = response.json()
        
        # Return the data or error
        if "data" in data:
            return data["data"]
        else:
            return {"error": "No data returned. Check the PDB ID or API response."}
    
    except requests.exceptions.RequestException as e:
        return {"error": str(e)}
    

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = args.file
 
    for pdb_id in [f.replace('.pdb', '').replace('.cif', '') for f in args.file]:
        if '_' in pdb_id:
            pdb_id, chain = pdb_id.split('_')
        header_info = fetch_pdb_header(pdb_id)
        title = header_info.get('entry', {}).get('struct', {}).get('title', 'Title not found')
        # Print the extracted title
        print(pdb_id, title)
