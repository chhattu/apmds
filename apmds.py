#!/usr/bin/python 

# Copyright (C)  2011 Zhang Initative Research Unit
# mutual_structural_align.py -description

# system 

import sys, os, commands, math, numpy

# mulitprocesing 
from multiprocessing import Process, Queue, current_process, freeze_support

# biopython

# user-defined modules
import compare, fileio, potential_region

out_dir = "" 

def open_pdb(input_pdb):
  line = open(input_pdb, "r")
  lines = line.read() 


# back string
def matchable_residues( ref_pdb, model_pdb ):
  # ref info
  cmd = "awk '{ if ($0 ~ /^ATOM/) print $4, $5, $6 }' " +  ref_pdb + " | uniq | head -1"
  ref_top_info = commands.getoutput(cmd).split()
  cmd = "awk '{ if ($0 ~ /^ATOM/) print $4, $5, $6 }' " +  ref_pdb + " | uniq | tail -1"
  ref_bottom_info = commands.getoutput(cmd).split()
  # model info
  cmd = "awk '{ if ($0 ~ /^ATOM/) print $4, $5, $6 }' " + model_pdb  + " | uniq | head -1"
  model_top_info = commands.getoutput(cmd).split()
  cmd = "awk '{ if ($0 ~ /^ATOM/) print $4, $5, $6 }' " + model_pdb + " | uniq | tail -1"
  model_bottom_info = commands.getoutput(cmd).split()
  # parse
  ref_top_res_name     = ref_top_info[0] 
  ref_top_chain        = ref_top_info[1] 
  ref_top_res_id       = ref_top_info[2] 
  ref_bottom_res_na    = ref_bottom_info[0]
  ref_bottom_chain     = ref_bottom_info[1]  
  ref_bottom_res_id    = ref_bottom_info[2]
  model_top_res_na     = model_top_info[0] 
  model_top_chain      = model_top_info[1] 
  model_top_res_id     = model_top_info[2] 
  model_bottom_res_na  = model_bottom_info[0]
  model_bottom_chain   = model_bottom_info[1]  
  model_bottom_res_id  = model_bottom_info[2]
  # check the flow of residues
  # make syntax like "[('A:5', 'A:10'), ('B:3', 'B:19')]"
  segment1 = "[('" + ref_top_chain    + ":" + ref_top_res_id    + "', "  + \
              "'" + ref_bottom_chain + ":" + ref_bottom_res_id + "')]" 
  segment2 = "[('" + model_top_chain + ":" + model_top_res_id  + "', "  + \
              "'" + model_bottom_chain + ":" + model_bottom_res_id  + "')]"
  return segment1 , segment2 

def make_writeable( model_name, overall_rmsd, rms2native, resi_rmsds):
  in_line = os.path.basename(model_name)
  in_line += " " + " %0.2f" %rms2native
  in_line += " " + " %0.2f" %overall_rmsd 
  for resi_rmsd in resi_rmsds:
    in_line += " %.2f" %(resi_rmsd)
  return in_line + "\n"

def make_writeable2native( model_name, rms2native, rmsds2native):
  in_line = os.path.basename(model_name)
  in_line += " " + " %0.2f" %rms2native
  for rmsd2native in rmsds2native:
    in_line += " %0.2f" %(rmsd2native)
  return in_line
  
def write_cc(model_id, cc, rms2native, overallrms):
  correlation = "%0.2f" %cc
  correlation += " " + "%0.2f" %rms2native
  correlation += " " + "%0.2f" %overallrms
  return model_id + " " + correlation + "\n"

'''
  xs rmsd2native, ys rmsd2models and its avearge
'''
def cal_corr(xs, ys):
  try:
    N = len(xs)
    if N != len(ys):
      print "Fatal error: residue mismatch in native and model" 
      sys.exit(1)
    if N == 0:
      print "Warning: there is no data to compute correlation."
      return 
    xbar = compute_mean ( xs, N)
    ybar = compute_mean ( ys, N)
    covXY = cal_covariance(xs, ys, xbar, ybar, N)
    sdx   = compute_sq_diff(xs, xbar, N)
    sdy   = compute_sq_diff(ys, ybar, N)
    cc    = covXY / (sdx * sdy)
    return round(cc, 2)
  except ZeroDivisionError:
    print "fatal: zero division error."
    sys.exit(1)
  except:
    print "fatal: unexpected error:", sys.exc_info()[0]
    sys.exit(1)

def cal_covariance(xs, ys, xbar, ybar, N):
  try:
    sumXY = 0.0
    for i in range(N):
      sumXY += ((xs[i] - xbar) * (ys[i] - ybar))
    return sumXY 
  except ZeroDivisionError:
    print "fatal: zero division error."
    sys.exit(1)
  except:
    print "fatal: unexpected error:", sys.exc_info()[0]
    sys.exit(1)

def compute_mean(xs, N):
  try:
    sumX = 0.0
    for i in range(N):
      sumX += xs[i]
    return sumX / N
  except ZeroDivisionError:
    print "fatal: zero division error."
    sys.exit(1)
  except:
    print "fatal: unexpected error:", sys.exc_info()[0]
    sys.exit(1)

def compute_sq_diff(xs, xbar , N):
  try:
    xdiff, xdiffsum = 0.0, 0.0
    for i in range(N):
      xdiff = (xs[i] - xbar)
      xdiffsum += (xdiff * xdiff)
    sd = math.sqrt(xdiffsum)
    return sd
  except ZeroDivisionError:
    print "fatal: zero division error."
    sys.exit(1)
  except:
    print "fatal: unexpected error:", sys.exc_info()[0]
    sys.exit(1)

def compute_stddev(xs):
  try:
    N = len(xs)
    mean = compute_mean(xs, N)
    xdiff, xdiffsum = 0.0, 0.0
    for i in range(N):
      xdiff = (xs[i] - mean)
      xdiffsum += math.pow(xdiff,2)
    sd = math.sqrt(xdiffsum/(N))
    return sd
  except ZeroDivisionError: # as (errno, strerror):
    print "fatal: divide by zero error." # ({0}): {1}" format(errno, strerror)
    sys.exit(1)
  except:
    print "fatal: unexpected error:", sys.exc_info()[0]
    sys.exit(1)


def compare2native( native_pdb, model_pdb):
  segment1, segment2 = matchable_residues(native_pdb, model_pdb)
  allresidues, overall = compare.get_best_alignment_with_residues( \
                     native_pdb, model_pdb, eval(segment1), eval(segment2), ['CA'])
  return allresidues, overall 

def compute_apmds(path2models, num_of_proc):
  args_results = {}
  segment1, segment2  = matchable_residues(path2models[0], path2models[1])

  if num_of_proc > 1:
    print "Using " + str(num_of_proc) + " cpus."
    args_results = parallel_process(path2models, segment1, segment2, num_of_proc)
    compute(args_results, path2models)
  else:
    print "Single CPU mode..."
    model_rmsds, protein_size, rmsd_matrix = fill_matrix(path2models, segment1, segment2)
    compute_apmds_in_one(model_rmsds, path2models, protein_size, rmsd_matrix)

def get_key(fullpath1, fullpath2):
  filename1 = os.path.splitext(os.path.basename(fullpath1))[0]
  filename2 = os.path.splitext(os.path.basename(fullpath2))[0]
  return str(filename1) + "-" + str(filename2)

def parallel_process(model_pdbs, segment1, segment2, nprocs):

  def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
      arg_result = {} 
      result     = func(*args)
      key = get_key(args[0], args[1])
      arg_result[key] = result 
      output.put(arg_result)

  NUMBER_OF_PROCESSES = nprocs 
  TASKS = []
  j = 0
  matrix_size = len(model_pdbs)
  t_queue = Queue()
  d_queue = Queue()

  for i in range(0, len(model_pdbs)):
    j       = i + 1
    ref_pdb = model_pdbs[i]
    while( j < matrix_size and len(model_pdbs[j]) > 0):
      model_pdb = model_pdbs[j]
      TASKS.append((compare.get_best_alignment_with_residues, 
                  (ref_pdb, model_pdb,  eval(segment1), eval(segment2),  ['CA'])))
      j = j + 1
    i = i + 1
 
  results = {}
  for task in TASKS: 
    t_queue.put(task)
  for i in range(NUMBER_OF_PROCESSES):
    Process(target=worker, args=(t_queue, d_queue)).start()
  for i in range(len(TASKS)):
    temp_results = d_queue.get() 
    for key in temp_results.keys(): 
      results[key] = temp_results[key] 
  for i in range(NUMBER_OF_PROCESSES):
    t_queue.put('STOP')

  return results  

def write2file(out_file, content):
  try:
    f = open(out_file, "w")
    f.write(content)
    f.close()
  except IOError:
    print "error: cannot produced output file."
    sys.exit(1)

def print2file(model_pdb, rsd_rmsd_avg, rsd_rmsd_std):
  line_to_be_printed = ""
  if(len(rsd_rmsd_avg) != len(rsd_rmsd_std)):
    print "Fatal: avg. and stdev. errors."
    sys.exit(1)

  for count in range(0, len(rsd_rmsd_avg)):
    line_to_be_printed += str(count + 1) + " "  + \
                          "%0.2f" %rsd_rmsd_avg[count] + " " \
                          "%0.2f" %rsd_rmsd_std[count] + "\n" 

  out_file = os.path.basename(model_pdb).replace(".pdb", ".dat")
  write2file(os.path.join(out_file), line_to_be_printed)

# def print_rsd2native(model_pdb, rsd_rms_avg_all, rsd_rms_std_all, rsd_rms_avg_native):
#   line_to_be_printed = ""
#   assert(len(rsd_rms_avg_all) == len(rsd_rms_avg_native))
#   assert(len(rsd_rms_avg_all) == len(rsd_rms_std_all))
# 
#   for count in range(0, len(rsd_rms_avg_all)):
#     line_to_be_printed += str(count + 1) + " "  + \
#                           "%0.2f" %rsd_rms_avg_all[count] + " " \
#                           "%0.2f" %rsd_rms_std_all[count] + " " \
#                           "%0.2f" %rsd_rms_avg_native[count] + "\n" 
# 
#   out_file = os.path.basename(model_pdb).replace(".pdb", ".dat")
#   write2file(os.path.join(out_file), line_to_be_printed)

# calculate avearge
def compute(args_result, path2models): 
  row = 0
  num_of_models = len(path2models) 
  apmds = {}

  # get keys
  first_key      = get_key(path2models[0], path2models[1]) 
  first_results  = args_result[first_key]
  protein_length = len(first_results[0])
  
  while (row < num_of_models): 
    column, ca_apmds  = 0, 0.0
    rsd_rmsd_sum = numpy.zeros((protein_length, 1), float) 
    rsd_rmsd_avg = numpy.zeros((protein_length, 1), float) 
    rsd_rmsd_std = numpy.zeros((protein_length, 1), float) 
    rsd_rmsd_std_matrix = numpy.zeros((protein_length, num_of_models), float)
     
    while (column < num_of_models):
      if (row != column):
        if row<column:
          key = get_key(path2models[row], path2models[column])
        elif column < row:
          key = get_key(path2models[column], path2models[row])
        results = args_result[key]
        rsd_rmsd = results[0]
        ca_apmds += results[1]
        for m in range(0, protein_length):
          rsd_rmsd_sum[m] += rsd_rmsd[m] 
          rsd_rmsd_std_matrix[m, column] = rsd_rmsd[m]
      column += 1

    for n in range(0, protein_length):
      rsd_rmsd_avg[n] = rsd_rmsd_sum[n] / (num_of_models-1) 
      rsd_rmsd_std[n] = numpy.std(numpy.delete(rsd_rmsd_std_matrix[n], row, 0))
    
    apmds[path2models[row]] = ca_apmds / (num_of_models-1)
    print2file(path2models[row], rsd_rmsd_avg, rsd_rmsd_std)
    row += 1 
  print_apmds2file(apmds)

def fill_matrix(model_pdbs, segment1, segment2):
  matrix_size = len(model_pdbs)
  rmsd_matrix = numpy.zeros((matrix_size, matrix_size), float)
  model_rmsds = {}
  j = 0
  
  for i in range(0, len(model_pdbs)):
    j       = i + 1
    ref_pdb = model_pdbs[i]
    while( j < matrix_size and len(model_pdbs[j]) > 0):
      model_pdb = model_pdbs[j]
      residue_rmsds, rmsd_matrix[i, j] = compare.get_best_alignment_with_residues( \
      ref_pdb, model_pdb, eval(segment1), \
      eval(segment2),  ['CA'])
      key              = str(i) + "-" + str(j)
      model_rmsds[key] = residue_rmsds
      protein_size     = len(residue_rmsds)
      j = j + 1
    i = i + 1
  # return model_rmsds, row, col , protein_size, rmsd_matrix
  return model_rmsds, protein_size, rmsd_matrix

def compute_apmds_in_one(model_rmsds, model_pdbs, protein_size, rmsd_matrix):
  avg_sum, row = 0.0, 0
  num_of_models = len(model_pdbs)
  apmds = {}
  while (row < num_of_models):
    rsd_rmsd_sum = numpy.zeros((protein_size, 1), float)
    rsd_rmsd_avg = numpy.zeros((protein_size, 1), float)
    
    rsd_rmsd_std = numpy.zeros((protein_size, 1), float)
    rsd_rmsd_std_matrix = numpy.zeros((protein_size, num_of_models), float)
    
    column, avg_sum = 0, 0.0
    cur_row, cur_col = 0, 0
    
    while (column < num_of_models):
      if (row != column):
        if(row < column):
          key = str(row) + "-" + str(column)
          cur_row = row
          cur_col = column
        elif (row > column):
          key = str(column) + "-" + str(row)
          cur_row = column
          cur_col = row
        rsd_rmsd = model_rmsds[key]
        for m in range(0, protein_size):
          rsd_rmsd_sum[m] += rsd_rmsd[m]
          rsd_rmsd_std_matrix[m, column] = rsd_rmsd[m] 
        avg_sum += rmsd_matrix[cur_row, cur_col]
      column = column + 1

    for n in range(0, protein_size):
      rsd_rmsd_avg[n]  = rsd_rmsd_sum[n] / (num_of_models - 1) 
      # since it does not compare with itself 
      rsd_rmsd_std[n] = numpy.std(numpy.delete(rsd_rmsd_std_matrix[n], row, 0))

    apmds[model_pdbs[row]] = avg_sum / (num_of_models - 1)
    print2file(model_pdbs[row], rsd_rmsd_avg, rsd_rmsd_std)
    row = row + 1
  print_apmds2file(apmds)


def print_apmds2file(apmds):
  try:
    content = ""
    keys    = apmds.keys()
    for key in keys:
      content += key + " " + "%0.2f" %apmds[key] + "\n"

    f = open("apmds.out", "w")
    f.write(content)
    f.close()
  except IOError:
    print "error: cannot produced output file."
    sys.exit(1)

# def compute_with_native(args_result, path2models, path2native):
#   row = 0
#   num_of_models = len(path2models) 
#   N = num_of_models - 1
# 
#   # get keys
#   first_key      = get_key(path2models[0], path2models[1]) 
#   first_results  = args_result[first_key]
#   protein_length = len(first_results[0])
#   
#   correlation, rms_native = "", "" 
#   while (row < num_of_models): 
#     column, rms2all = 0, 0.0;
#     rsd_rms_sum = numpy.zeros((protein_length, 1), float) 
#     rsd_rms_avg = numpy.zeros((protein_length, 1), float) 
#     rsd_rms_std = numpy.zeros((protein_length, 1), float) 
#     rsd_rms_std_mat = numpy.zeros((protein_length, num_of_models), float)
# 
#     while (column < num_of_models):
#       if (row != column):
#         if row<column:
#           key = get_key(path2models[row], path2models[column])
#         elif column < row:
#           key = get_key(path2models[column], path2models[row])
#         results = args_result[key]
#         rsd_rms = results[0]
#         rms2all += results[1]
#         for m in range(0, protein_length):
#           rsd_rms_sum[m] += rsd_rms[m] 
#           rsd_rms_std_mat[m, column] = rsd_rms[m] 
#       column += 1
# 
#     for n in range(0, protein_length):
#       rsd_rms_avg[n] = rsd_rms_sum[n] / N 
#       rsd_rms_std[n] = numpy.std(rsd_rms_std_mat[n]) 
# 
#     rsd_rms2native , rms2native = compare2native( path2native, path2models[row])
#     p_correlation  = cal_corr(rsd_rms2native, rsd_rms_avg)
#     correlation += write_cc(path2models[row], p_correlation, rms2native, rms2all/N) 
#     rms_native += make_writeable2native(path2models[row], rms2native, rsd_rms2native) 
#     print_rsd2native(path2models[row], rsd_rms_avg, rsd_rms_std, rsd_rms2native)
#     row += 1 
# 
#   # print rmsd_matrix
#   write2file(os.path.join("correlation.out"), correlation)
#   write2file(os.path.join("rmsd2native.out"), rms_native)

def read_model_list( in_file_list ):
  files = []
  try:
    f_lists    = open(in_file_list, "r")
    temp_files = f_lists.read().split("\n")
    for temp_file in temp_files:
      if len(temp_file) > 0:
        files.append(temp_file)
    f_lists.close()
    if len(files) <= 2:
      print "Fatal: not enough input models."
      sys.exit(0)
  except IOError:
    print "error: no file pdb list."
    sys.exit(0)
  return files


if __name__=="__main__":

  if len(sys.argv) < 2 or len(sys.argv) > 9:
    print "Usage:" 
    print "apmds.py pdb_list(decoys) [ -p 3 ]"
    print " -p number of processor."  
    sys.exit(0)

  path2models = sys.argv[1]

  opt_arg = 2
  argu    = {}
  while opt_arg < len(sys.argv):
    argu[sys.argv[opt_arg]] = sys.argv[opt_arg+1]
    opt_arg += 2

  num_proc  = 1 
  for key in argu.keys():
    if key == "-p": num_proc    = int(argu["-p"]) 

  files = read_model_list(path2models)
  print "Please wait..."
  compute_apmds(files, num_proc)
  print "Finished."
