#!/usr/bin/env python
# -*- coding: utf-8 -*-
import fileinput

for line in fileinput.input():
    line = reversed(line.strip())
    print(''.join(line))
