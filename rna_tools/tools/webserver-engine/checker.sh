#!/bin/bash
ssh aws 'cd /home/ubuntu/rna-tools/rna_tools/tools/webserver-engine/media/jobs && /home/ubuntu/miniconda3/bin/python check.py' > check.md
emacsclient check.md
