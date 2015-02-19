#!/usr/bin/python

"""
Set of routines to calculate the RMSD between two molecular structures.
The module can be run from the command line using PDB files as input.
"""

import math
import numpy
import vector3d, util, molecule, polymer


def rmsd(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array"""

  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)

  n_vec = numpy.shape(crds1)[0]
  correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
  v, s, w = numpy.linalg.svd(correlation_matrix)
  is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  return numpy.sqrt(rmsd_sq)


def optimal_superposition(crds1, crds2):
  """Returns best-fit rotation matrix as [3x3] numpy matrix"""
  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)
  correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
  v, s, w = numpy.linalg.svd(correlation_matrix)
  is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0
  if is_reflection:
    v[-1,:] = -v[-1,:]
  return numpy.dot(v, w)


def get_i_residue(residues, tag):

  def get_tag(residue):
    tag = ""
    if residue.chain_id != " " and residue.chain_id != "":
      tag += residue.chain_id + ":"
    tag += str(residue.num)
    if residue.insert:
      tag += residue.insert
    return tag  

  # clean up tag
  tag = tag.strip()
  if tag[0] == ":":
    tag = tag[1:]
  if not tag[0].isdigit() and tag[1].isdigit():
    tag = tag[0] + ":" + tag[1:]

  for i, residue in enumerate(residues):
    if tag.lower() == get_tag(residue).lower():
      return i
  raise "Can't find residue", tag

def get_superposable_all_atoms(polymer, segments): 
  atom_types=['CA', 'N', 'C', 'CB',
              'CD', 'CD1', 'CD2', 'CE', 'CE1',
              'CE2', 'CE3', 'CG', 'CG1', 'CG2',
              'CH2', 'CZ', 'CZ2', 'CZ3', 'ND1', 
              'ND2', 'NE', 'NE1', 'NE2', 'NH1', 
              'NH2', 'NZ', 'O', 'OD1', 'OD2', 
              'OE1', 'OE2', 'OG', 'OG1', 'OH', 
              'OXT', 'SD', 'SG', 'CEN'] 
  result = []
  allowed_i = []
  residues = polymer.residues()
  for res_num_i, res_num_j in segments:
    i = get_i_residue(residues, str(res_num_i))
    j = get_i_residue(residues, str(res_num_j))
    allowed_i.extend(range(i,j+1))
  for i, residue in enumerate(residues):
    if i in allowed_i:
      result.extend([a for a in residue.atoms()
                     if a.type in atom_types])
  return result


def get_superposable_atoms(polymer, segments, 
           atom_types=['CA', 'N', 'C', 'CB','O']):
  result = []
  allowed_i = []
  residues = polymer.residues()
  for res_num_i, res_num_j in segments:
    i = get_i_residue(residues, str(res_num_i))
    j = get_i_residue(residues, str(res_num_j))
    allowed_i.extend(range(i,j+1))
  for i, residue in enumerate(residues):
    if i in allowed_i:
      result.extend([a for a in residue.atoms()
                     if a.type in atom_types])
  return result


def get_crds(atoms):
  crds = numpy.zeros((len(atoms), 3), float)
  for i, a in enumerate(atoms):
    crds[i,0] = a.pos.x
    crds[i,1] = a.pos.y
    crds[i,2] = a.pos.z
  return crds


def calculate_superposition_matrix(atoms1, atoms2):

  def convert_to_matrix3d(numpy_matrix3d):
    result = vector3d.Matrix3d()
    for i in range(3):
      for j in range(3):
        result.setElem(i, j, numpy_rotation[j, i])
    return result

  numpy_rotation = optimal_superposition(get_crds(atoms1), get_crds(atoms2))
  return convert_to_matrix3d(numpy_rotation)
    

def sum_rmsd(atoms1, atoms2):
  sum_squared = 0.0
  for atom1, atom2 in zip(atoms1, atoms2):
    sum_squared += vector3d.pos_distance(atom1.pos, atom2.pos)**2
  return math.sqrt(sum_squared/float(len(atoms1)))
  
def residues_rmsd(atoms1, atoms2):
  sum_squared = [] 
  for atom1, atom2 in zip(atoms1, atoms2):
    each_diff = math.sqrt(vector3d.pos_distance(atom1.pos, atom2.pos)**2)
    sum_squared.append(round(each_diff, 2))
  return sum_squared, round(sum_rmsd(atoms1, atoms2), 2) 

def get_raw_rmsd(pdb1, pdb2, segments1, segments2, atom_types):
  polymer1 = polymer.Polymer(pdb1)
  polymer2 = polymer.Polymer(pdb2)
  atoms1 = get_superposable_atoms(polymer1, segments1, atom_types)
  atoms2 = get_superposable_atoms(polymer2, segments2, atom_types)
  
  return sum_rmsd(atoms1, atoms2)


def get_best_alignment(pdb1, pdb2, segments1, segments2, atom_types):
  """Returns rmsd and filename of transformed pdb2."""
  polymer1 = polymer.Polymer(pdb1)
  atoms1   = get_superposable_atoms(polymer1, segments1, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2   = get_superposable_atoms(polymer2, segments2, atom_types)
  center1  = molecule.get_center(atoms1)
  polymer1.transform(vector3d.Translation(-center1))
  polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))
  polymer2.transform(calculate_superposition_matrix(atoms1, atoms2))

  rmsd = sum_rmsd(atoms1, atoms2)
  return rmsd

def get_best_alignment_with_residues(pdb1, pdb2, segments1, segments2, atom_types):
  polymer1 = polymer.Polymer(pdb1)
  atoms1   = get_superposable_atoms(polymer1, segments1, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2   = get_superposable_atoms(polymer2, segments2, atom_types)
  center1  = molecule.get_center(atoms1)
  polymer1.transform(vector3d.Translation(-center1))
  polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))
  polymer2.transform(calculate_superposition_matrix(atoms1, atoms2))

  residue_rmsd, overall_rmsd = residues_rmsd(atoms1, atoms2)
  return residue_rmsd, overall_rmsd  

def get_best_alignment_with_all_residues(pdb1, pdb2, segments1, segments2, atom_types):
  polymer1 = polymer.Polymer(pdb1)
  atoms1   = get_superposable_atoms(polymer1, segments1, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2   = get_superposable_atoms(polymer2, segments2, atom_types)
  center1  = molecule.get_center(atoms1)
  polymer1.transform(vector3d.Translation(-center1))
  polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))
  polymer2.transform(calculate_superposition_matrix(atoms1, atoms2))

  allatoms1   = get_superposable_all_atoms(polymer1, segments1)
  allatoms2   = get_superposable_all_atoms(polymer2, segments1)

  residue_rmsd, overall_rmsd = residues_rmsd(allatoms1, allatoms2)
  return residue_rmsd, overall_rmsd  

def get_total_num_atoms(pdb1, segment1):
  polymer1 = polymer.Polymer(pdb1)
  atoms1   = get_superposable_all_atoms(polymer1, segment1)
  return len(atoms1)

def get_rmsd(pdb1, pdb2, segments1, segments2, atom_types):
  polymer1 = polymer.Polymer(pdb1)
  atoms1 = get_superposable_atoms(polymer1, segments1, atom_types)
  polymer2 = polymer.Polymer(pdb2)
  atoms2 = get_superposable_atoms(polymer2, segments2, atom_types)

  center1 = molecule.get_center(atoms1)
  polymer1.transform(vector3d.Translation(-center1))
  polymer2.transform(vector3d.Translation(-molecule.get_center(atoms2)))

  crds1 = get_crds(atoms1)
  crds2 = get_crds(atoms2)
  return rmsd(crds1, crds2)


def segments_str(segments):
  residues = []
  for i, j in segments:
    if i == j:
      residues.append(str(i))
    else:
      residues.append("%s-%s" % (i,j))
  return ', '.join(residues)
  

