#!/usr/bin/python

# Copyright (C)  2011 Zhang Initative Research Unit
# fileio.py -description

# system 
import sys, os

# biopython

# user-defined modules


def write_2_file( out_file, content ):
  f = open(out_file, "w")
  f.write(content)
  f.close()
