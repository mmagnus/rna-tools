#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

SECOND_PK_CHAR = ['{', '}']

def get_one_line(lines):
    if len(lines) == 2:
            ss = ''
            for i,j in zip(lines[0], lines[1]):
                if j == '(':
                    i = '['
                if j == ')':
                    i = ']'
                ss += i
    if len(lines) == 3:
            ss = ''
            for i,j,k in zip(lines[0], lines[1], lines[2]):
                if j == '(':
                    i = '['
                if j == ')':
                    i = ']'
                if k == '(':
                    i = SECOND_PK_CHAR[0]
                if k == ')':
                    i = SECOND_PK_CHAR[1]
                ss += i
    return ss

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='file')
    args = parser.parse_args()
    if args.file:
        ss = open(args.file)
        print((ss.readline().strip()))
        lines = ss.readlines()
        print((get_one_line(lines)))
        
