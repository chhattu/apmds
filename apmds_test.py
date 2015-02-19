#!/usr/bin/python

# Copyright (C)  2011 Zhang Initative Research Unit
# mutual_test.py -description

# system 
import compare, apmds
import math
import unittest

class rmsdtest(unittest.TestCase):

  def rmsdToSelf( self ):
    pdb1 = './test/model1.pdb' 
    pdb2 = './test/model1.pdb'
    segment1 = eval("[('A:1', 'A:128')]")
    segment2 = segment1 
    residue_rmsds, overall_rmsd = \
                      compare.get_best_alignment_with_residues( 
                      pdb1, pdb2, segment1, segment2, ['CA'])
    self.assertEqual(overall_rmsd, 0)

  def compare2native( self ):
    native_pdb = './test/native.pdb'
    model_pdb = './test/model1.pdb'
    residues, overall = apmds.compare2native( native_pdb, model_pdb )
    self.assertEqual(overall, 3.11)

  def stddevTest(self):
    xs = [10, 20, 50, 100, 150, 200, 230]
    stdev = apmds.compute_stddev(xs)
    self.assertEqual(stdev, 87.45067) 
 
  def corTest(self):
    xs = [10, 20, 50, 100, 150, 200, 230]
    ys = [10, 20, 50, 100, 150, 200, 230]
    corr = apmds.cal_corr(xs, ys)
    self.assertEqual(corr, 1.0)

    xs = [10, 20, 50, 100, 150, 200, 230]
    ys = [230, 200, 150, 100, 50, 20, 10 ]
    corr = round(apmds.cal_corr(xs, ys), 2)
    self.assertEqual(corr, -0.98)

  def matrix(self):
    native_pdb = "./test/native.pdb"
    models = []
    models.append("./test/model3.pdb")
    models.append("./test/model4.pdb")
    models.append("./test/model5.pdb")
    models.append("./test/model1.pdb")
    segment1, segment2 = apmds.matchable_residues(models[0], models[1])  
    results = apmds.parallel_process(models, segment1, segment2, 1)
    keys =  results.keys()
    # order has meaning
    known_keys     = ['model5-model1', 'model4-model1', 'model3-model5', 
                      'model3-model4', 'model3-model1', 'model4-model5']
    ca_rmsds_known = [6.56, 6.53, 4.71, 4.58, 6.47, 4.62]

    for i in range(0, len(known_keys)):
      self.assertEqual(keys[i], known_keys[i])
      rmsds = results[known_keys[i]]
      self.assertEqual( round(rmsds[1], 2), ca_rmsds_known[i]) 
   
    apmds.compute(results, models)

  def statistics(self):

    def read_data(filepath, delimiter):
      try:
        mean = {} 
        stdev = {}
        f = open(filepath)
        lines = f.read().split("\n")
        for line in lines:
          if line:
            values = line.split(delimiter)
            mean[values[0]] = values[1]
            stdev[values[0]] = values[2]
        f.close()
        return mean, stdev
      except IOError:
        print "fatal: file could not be opened."

    model1_mean, model1_sd = read_data("./model1.dat", " ")
    model1_mean_known, model1_sd_known = read_data("./test/model1_valid.csv",",")
    self.assertEqual(len(model1_mean), len(model1_mean_known)) 
    self.assertEqual(len(model1_sd), len(model1_sd_known)) 

    for key in model1_mean.keys():
      self.assertEqual(math.ceil(float(model1_mean[key])), math.ceil(float(model1_mean_known[key])))
      self.assertEqual(math.ceil(float(model1_sd[key])), math.ceil(float(model1_sd_known[key]))) 

    model3_mean, model3_sd = read_data("./model3.dat", " ")
    model3_mean_known, model3_sd_known = read_data("./test/model3_valid.csv",",")
    self.assertEqual(len(model3_mean), len(model3_mean_known)) 
    self.assertEqual(len(model3_sd), len(model3_sd_known)) 

    for key in model3_mean.keys():
      self.assertEqual(math.ceil(float(model3_mean[key])), math.ceil(float(model3_mean_known[key])))
      self.assertEqual(math.ceil(float(model3_sd[key])), math.ceil(float(model3_sd_known[key]))) 

  def matrix_native(self):
    native_pdb = "./test/native.pdb"
    models = []
    models.append("./test/model3.pdb")
    models.append("./test/model4.pdb")
    models.append("./test/model5.pdb")
    models.append("./test/model1.pdb")
    segment1, segment2 = apmds.matchable_residues(models[0], models[1])  
    results = apmds.parallel_process(models, segment1, segment2, 1)
    keys =  results.keys()
    # order has meaning
    known_keys     = ['model5-model1', 'model4-model1', 'model3-model5', 
                      'model3-model4', 'model3-model1', 'model4-model5']
    ca_rmsds_known = [6.56, 6.53, 4.71, 4.58, 6.47, 4.62]

    for i in range(0, len(known_keys)):
      self.assertEqual(keys[i], known_keys[i])
      rmsds = results[known_keys[i]]
      self.assertEqual( round(rmsds[1], 2), ca_rmsds_known[i]) 
   
    apmds.compute_with_native(results, models, native_pdb)

  def statistics_native(self):

    def read_data(filepath, delimiter):
      try:
        mean = {} 
        stdev = {}
        f = open(filepath)
        lines = f.read().split("\n")
        for line in lines:
          if line:
            values = line.split(delimiter)
            mean[values[0]] = values[1]
            stdev[values[0]] = values[2]
        f.close()
        return mean, stdev
      except IOError:
        print "fatal: file could not be opened."

    model1_mean, model1_sd = read_data("./model1.dat", " ")
    model1_mean_known, model1_sd_known = read_data("./test/model1_valid.csv",",")
    self.assertEqual(len(model1_mean), len(model1_mean_known)) 
    self.assertEqual(len(model1_sd), len(model1_sd_known)) 

    for key in model1_mean.keys():
      self.assertEqual(math.ceil(float(model1_mean[key])), math.ceil(float(model1_mean_known[key])))
      self.assertEqual(math.ceil(float(model1_sd[key])), math.ceil(float(model1_sd_known[key]))) 

    model3_mean, model3_sd = read_data("./model3.dat", " ")
    model3_mean_known, model3_sd_known = read_data("./test/model3_valid.csv",",")
    self.assertEqual(len(model3_mean), len(model3_mean_known)) 
    self.assertEqual(len(model3_sd), len(model3_sd_known)) 

    for key in model3_mean.keys():
      self.assertEqual(math.ceil(float(model3_mean[key])), math.ceil(float(model3_mean_known[key])))
      self.assertEqual(math.ceil(float(model3_sd[key])), math.ceil(float(model3_sd_known[key]))) 

if __name__=="__main__":
  testSuite = unittest.TestSuite()
  testSuite.addTest(rmsdtest("compare2native"))
  testSuite.addTest(rmsdtest("rmsdToSelf"))
  testSuite.addTest(rmsdtest("corTest"))
  testSuite.addTest(rmsdtest("matrix"))
  testSuite.addTest(rmsdtest("statistics"))
  # has to test
  # testSuite.addTest(rmsdtest("matrix_native"))
  # testSuite.addTest(rmsdtest("statistics_native"))

  unittest.TextTestRunner(verbosity=2).run(testSuite)
