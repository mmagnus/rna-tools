##!/bin/sh
# curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py

#ABSPATH=$(cd "$(dirname "$0")"; pwd)
#. $ABSPATH/../bin/activate
#source rna_tools_env/bin/activate
kill -9 `lsof -t -i:8080`
python3 manage.py runserver --settings web.settings  localhost:8080
