#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='file')
    args = parser.parse_args()
    if args.file:
        ss = open(args.file)
        print ss.readline().strip()
        lines = ss.readlines()
        if len(lines) == 2:
            ss = ''
            for i,j in zip(lines[0], lines[1]):
                if j == '(':
                    i = '['
                if j == ')':
                    i = ']'
                ss += i
            print ss.strip()
        if len(lines) == 3:
            ss = ''
            for i,j,k in zip(lines[0], lines[1], lines[2]):
                if j == '(':
                    i = '['
                if j == ')':
                    i = ']'
                if k == '(':
                    i = '<'
                if k == ')':
                    i = '>'
                ss += i
            print ss.strip()
            
