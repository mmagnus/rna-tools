#!/bin/sh

# curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py

ABSPATH=$(cd "$(dirname "$0")"; pwd)

. $ABSPATH/../bin/activate

python manage.py runserver --settings web.settings  0:8080
