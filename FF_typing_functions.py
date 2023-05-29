#!/usr/bin/python

###################################
# 
# forcefield typing module Eric Chan Jan 2105
# use open babel bindings to facilitie assignmnet of forcefield types
# built on desision trees 
#
# These software tools are currently
# Licensed under the Academic Free License version 3.0
# The authors hold no responsibility whatsoever for any software 
# or hardware related damages that can occur as a result of the use of these codes.
#
# See https://opensource.org/licenses/AFL-3.0 for the details.
#
#####################################

import re
import sys
import os

import numpy as np
import openbabel

# send in molecule and atom, return FF_type and atom_symbol

def get_FFtype_DREIDING(atom):

# The dreiding FF is reletivly easy to make assignments but the current method is not full prof 
# right now and needs to be checked
# often the hybridiation of the atom in babel is not correct and the wrong FFtype can be assigned  
# so you have to put in extra statements to help fix for this 
# if the hybridisation assignment is influenced by the geometry then the hybridization is wrong
# babel may rely exclusivly on geometry to figure out connectivity and atom typing  

# basically does the typing how I want it to go 

    hybr = atom.GetHyb()  # get the hybridization   
                      
    atomicnum = atom.GetAtomicNum() # get the atomic number 

    if atom.IsHydrogen():
     atom_symbol = "H"  
     if hybr == 0 : FFtype = "H_"
     if atom.IsHbondDonorH(): FFtype = "H___A"
     # check how many bonds are connected to the hydrogen. if greater than 1 then it is bridging 
     if np.size(openbabel.OBAtomAtomIter(atom)) > 1: FFtype = "H___b"

     # print hybr
     # for jatom in openbabel.OBAtomAtomIter(atom):
     # bond = atom.GetBond(jatom)
     # print bond.GetBondOrder()
     
     if hybr != 0 : print "check your hydrogen conectivities"

    if atomicnum == 5:
     atom_symbol = "B"  
     if hybr == 2 : FFtype = "B_2"
     if hybr == 3 : FFtype = "B_3"
     if not FFtype:  print "check your hydrogen conectivities"  
      
    if atom.IsCarbon():   
     atom_symbol = "C"  
     if hybr == 1 : FFtype = "C_1"
     if hybr == 2 : FFtype = "C_2"
     if hybr == 3 : FFtype = "C_3"
     if atom.IsAromatic(): FFtype = "C_R"

    if atom.IsNitrogen():   
     atom_symbol = "N"  
     if hybr == 1 : FFtype = "N_1"
     if hybr == 2 : FFtype = "N_2"
     if hybr == 3 : FFtype = "N_3"
     if atom.IsAromatic(): FFtype = "N_R"
     if atom.IsAmideNitrogen(): FFtype = "N_R"

    if atom.IsOxygen():  
     atom_symbol = "O"  
    # print hybr
     if hybr == 1 : FFtype = "O_1"
     if hybr == 2 : 
      FFtype = "O_2"
     # not currently sure how to deal with this it may be best to just treat all Carboxl oxygens as SP2
     # I think babel automatically makes the check for Carboxlate H-donors and assigns SP2  
     # fore esters the oxygen is supposedly SP2  -- see https://classesv2.yale.edu/wiki/site/chem124_f08/ester.html 
     # should try do the theoretical calc
      if not atom.HasDoubleBond() : 
       FFtype = "O_3" # in most case SP2 oxygens should have double bonds or are aromatic
       if atom.IsCarboxylOxygen(): FFtype = "O_R" # if it is carboxylic acid hydroxyl group then make it have aromatic character
     if hybr == 3 : FFtype = "O_3"
     if atom.IsAromatic(): FFtype = "O_R"
     if atom.IsNitroOxygen(): FFtype = "O_R"
     if atom.IsSulfateOxygen(): FFtype = "O_2"
     if atom.IsPhosphateOxygen(): FFtype = "O_2"
     # if atom.IsEster(): print "ester" 
     # you still need to check if the oxygen is an ester bond then assign it SP3  

    if atomicnum == 9:
     atom_symbol = "F"  
     FFtype = "F_"

    if atomicnum == 13:
     atom_symbol = "Al"  
     FFtype = "Al3"

    if atomicnum == 14:
     atom_symbol = "Si"  
     FFtype = "Si3"

    if atomicnum == 15:
     atom_symbol = "P"  
     FFtype = "P_3"

    if atomicnum == 16:
     atom_symbol = "S"  
     FFtype = "S_3"

    if atomicnum == 17:
     atom_symbol = "Cl"  
     FFtype = "Cl"

    if atomicnum == 31:
     atom_symbol = "Ga"  
     FFtype = "Ga3"

    if atomicnum == 32:
     atom_symbol = "Ge"  
     FFtype = "Ge3"

    if atomicnum == 33:
     atom_symbol = "As"  
     FFtype = "As3"

    if atomicnum == 34:
     atom_symbol = "Se"  
     FFtype = "Se3"

    if atomicnum == 35:
     atom_symbol = "Br"  
     FFtype = "Br"

    if atomicnum == 49:
     atom_symbol = "In"  
     FFtype = "In3"

    if atomicnum == 50:
     atom_symbol = "Sn"  
     FFtype = "Sn3"

    if atomicnum == 51:
     atom_symbol = "Sb"  
     FFtype = "Sb3"

    if atomicnum == 52:
     atom_symbol = "Te"  
     FFtype = "Te3"

    if atomicnum == 53:
     atom_symbol = "I"  
     FFtype = "I_"

    if atomicnum == 11:
     atom_symbol = "Na"  
     FFtype = "Na"

    if atomicnum == 20:
     atom_symbol = "Ca"  
     FFtype = "Ca"

    if atomicnum == 22:
     atom_symbol = "Ti"  
     FFtype = "Ti"

    if atomicnum == 26:
     atom_symbol = "Fe"  
     FFtype = "Fe"

    if atomicnum == 30:
     atom_symbol = "Zn"  
     FFtype = "Zn"

    if atomicnum == 43:
     atom_symbol = "Tc"  
     FFtype = "Tc"

    if atomicnum == 44:
     atom_symbol = "Ru"  
     FFtype = "Ru"

#
# this is use of the "try" and exception error handling components of python 
# see http://stackoverflow.com/questions/1592565/determine-if-variable-is-defined-in-python 
#

    try:
      FFtype
    except NameError:
      print "there is no FF type described for this atom using this FF"
      print "check the input for errors "
 #   else:
 #     print "FFtype is defined."

    return FFtype, atom_symbol


def main():

 from sys import argv
 script, fname = argv

 #  use babel to get the data we want
 obConv = openbabel.OBConversion()
 obConv.SetInFormat("arc")

 mol = openbabel.OBMol()       
 obConv.ReadFile(mol, fname)   
 print mol.GetMolWt()

# get a single atom from the molecule by index (starts at 1)
# atom = mol.GetAtom(1)       

# or loop through all atoms 
 for atom in openbabel.OBMolAtomIter(mol): 
 # atomidx = obatom.GetIdx() 

  FFtype, atom_symbol = get_FFtype_DREIDING(atom)
  print "%s %s" %(FFtype, atom_symbol)

# atomnames,resnames,resnos,FF_types,atomtypes = readmol_from_arcfile(fname,0) 
# print atomnames,resnames,resnos,FF_types,atomtypes 

if __name__ == "__main__":
    main()




