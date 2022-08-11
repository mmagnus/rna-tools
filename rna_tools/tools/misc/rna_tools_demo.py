#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from zipfile import ZipFile
try:
    import wget
except:
    print('wget missing, please "pip install wget"')
    sys.exit(1)
             
url = "http://rna-tools.online/static/app/demo/rp21.zip"
f = wget.download(url)
fn = 'rp21.zip'
with ZipFile(fn, 'r') as zip:
    zip.printdir()
    zip.extractall()
