#!/usr/bin/python

# Copyright (C)  2011 Zhang Initative Research Unit
# potential_region.py -description

# system 
import sys, os

# biopython

# user-defined modules
import fileio


def potintial_wrong_region( model_pdb, out_dir, avg_rms, mutual_rmsds):
  line_to_be_printed = ""
  print model_pdb + " " + str(avg_rms)
  for count in range(len(mutual_rmsds)):
    line_to_be_printed += str(count + 1) + " "  + "%0.2f" %mutual_rmsds[count] + "\n" 
    count += 1
  out_file = os.path.basename(model_pdb).replace(".pdb", ".dat")
  fileio.write_2_file(os.path.join(out_dir,out_file), line_to_be_printed)

def potential_wrong_regions_native_rmsds(model_pdb, out_dir, native_rmsds, mutual_rmsds):
  line_to_be_printed = ""
  for count in range(len(mutual_rmsds)):
    line_to_be_printed += str(count + 1) + " "  + "%0.2f" %mutual_rmsds[count] + " " + \
                                                  "%0.2f" %native_rmsds[count] + "\n" 
    count += 1
  out_file = os.path.basename(model_pdb).replace(".pdb", ".dat")
  fileio.write_2_file(os.path.join(out_dir, out_file), line_to_be_printed)

  
