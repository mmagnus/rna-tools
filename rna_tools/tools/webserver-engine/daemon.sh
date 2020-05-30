#!/bin/bash
cd /home/rnamasonry/rnamasonryweb_env/rnamasonry-web
source ../bin/activate && python daemon.py --prod --verbose $@
