#!/opt/local/bin/python2.7

# 
#
# These software tools are currently
# Licensed under the Academic Free License version 3.0
# The authors hold no responsibility whatsoever for any software 
# or hardware related damages that can occur as a result of the use of these codes.
#
# See https://opensource.org/licenses/AFL-3.0 for the details.
#
# Eric Chan Jan 2016
#--------------------
#
# october 2016 - modified to run antechamber to perform FF assignments as well and copy over the user assigned charges   
# you can now place a 1 or 0 if you want to mol2car to run antechamber  
# 
# for 2016 version 
# added in the option to do 'drd' = dreiding FF
#                        or 'gaff' = genralised amber  FF  
#
# for output using gaff you need two mol2 files. the original and one that comes with the amber fftypes. 
# you need to make sure this file is in the same formating as the main input mol2  
#
# if we read just the one for amber then obabel function will crash.
# so we need to do typing in antechamber and then just use both mol2 files to create the .car and .mdf 
#
# the name of this extra mol2 file should always be rootname_gaff.mol2, the program will look for this file.
#
# so prior to setting up the gaff implemetation in lammps you will need to run 
#
# eg. 
# antechamber -fi mol2 -fo mol2 -i rootname.mol2 -o rootname_gaff.mol2 -c gas  
# 
# the running of this program is then exactly the same 
#
# for more complex structures you might have to do some further scripting that so that the input .mol2 files  
# are formatted  correctly so they can be converted by the program. 
#
# the program mol2car is not limmited in the number of molecules that can be read in but does not describe  
# residues in the manner that protein structures are often utilised. 
#
# the strategy is to read in the .mol2 file as text as well as read the data into babelclasses 
# so that we can do interesting thing=s by using babel functions 
# it is based on the original script for reading .arc files that are output from GDIS only work for .mol2 
# format 
#
# you can concatonate  multiple molecules in the same .mol2 file and it will still work 
# you should only have to specify the unit cell at the end if there is one.
# there should be the header @<TRIPOS>MOLECULE  
#
# this means the standard  .mol2 output form mercury and openbabel should work 
# and in many cases babel can be used to convert to .mol2 
# 
# note that the purpose of this is to create the input .car and .mdf files for MSI2LMP.exe
# to create the LAMMPS input files  
# note that LAMMPS does not differentiate residue numbers so reading in residue numbers is not actually nessesary  
# 
# run by using python2.7 mol2car.py rootname 
#
# if you want seperate molecule numbers for the molecules in your unit cell use eg.
#  the separate residues will be placed into the format for doing this 
# 
# this provides some versatility with the organisation of molecule numbers    
# 
#  babel --separate Acetylsalicylic_acid_CSD.mol2 test.mol2
# 
# it is the opposite of babel -join 
#
 
import re
import sys
import os

import numpy as np
import openbabel
import read_moldata as rdat
import FF_typing_functions as gft
import time

# import matplotlib.pyplot as plt
# import matplotlib

# import scipy
# import matplotlib  


def read_write_mol_data(mol,carfile_handle,mdffile_handle,molno,ffname,mol2file):

 # in a seperate function lets read in the .arc file and get the atom data we want as lists 
 # atomnames,resnames,resnos,FF_types,atomtypes = rdat.readmol_from_arcfile(arcfile,molno) 
 
 mollabel=molno+1
 #write the molecule number to the .mdf file 
 mdffile_handle.write("\n@molecule mol %d\n\n" %(mollabel) )
  
 # mol.AddHydrogens()

# for gaff, get the FFtypes and symobls for this molecule from the mol2 file
 if ffname == 'gaff' : fftypes,symbs,gaffcharges = rdat.get_fftype_data_from_mol2file(mol2file,mollabel)
 
 # print mol.NumAtoms()
 # print mol.NumBonds()
 # print mol.NumResidues()
 
 mol.UnsetPartialChargesPerceived();
 for obatom in openbabel.OBMolAtomIter(mol): 
  atomidx = obatom.GetIdx() 
  atomid = obatom.GetId() 
  atomNo1 = obatom.GetAtomicNum() 
  charge = obatom.GetPartialCharge()           
  formal_charge = obatom.GetFormalCharge()           
  isotope = obatom.GetIsotope()           
  atype = obatom.GetType()           
  xpos =  obatom.GetX() 
  ypos =  obatom.GetY() 
  zpos =  obatom.GetZ() 
 # print "%3i %3i %s %6f %6f %6f %6f " %(atomid, atomNo1, atype, charge, xpos, ypos, zpos)
 #  this will automatically recalculate charges using gast method 
 # I dont really know how to change from gastieger charges
 # but you are supposed to use the OBChargeModel class which is not straightforward to implement
  
 # example of how to iterate over atoms and bond types 
 # to construct the .MDF file we can use property requests on each bond connecting each atom  
 #
 # it turns out that msi2lmp if very sensitive to the order of the  conectivity 
 #
 # in some cases the force field types can be read directly from the .mol2 file eg. gaff forcefield
 # in these cases it is best to just do the typing directly from labels in the .mol2 file  
 # we can do this on a per molecule basis sicne this is the way that mol2car has been designed to run 

  if ffname == 'drd' : FFtype1, atom_symbol1 = gft.get_FFtype_DREIDING(obatom)
  if ffname == 'gaff' : 
    FFtype1 =  fftypes[atomid]  
    atom_symbol1 = symbs[atomid]
    charge = float(gaffcharges[atomid])           

  na_string_list = []   # dump all the neighbor atom strings into a list 
  na_index_list = []   # dump all the neighbor atom strings into a list 
  for neighbour_atom in openbabel.OBAtomAtomIter(obatom):
     atomNo2 = neighbour_atom.GetAtomicNum()
     atomid2 = neighbour_atom.GetId()
     atomidx2 = neighbour_atom.GetIdx()
     bond = obatom.GetBond(neighbour_atom)
     bond_order = bond.GetBondOrder()
     if bond.IsAromatic(): bondtype = "/1.5" # bondtype = "aromatic"
     if bond.IsSingle():   bondtype = "" # bondtype = "single"
     if bond.IsDouble():   bondtype = "/2.0" # bondtype = "double"
     if bond.IsTriple():   bondtype = "/3.0" # bondtype = "triple"  
 
 #   i dont think we need this for the .mdf format 
 #    if bond.IsInRing(): print "ring"


     try:
      bondtype
     except NameError:
      print "cannot determine the bondtype - should assign manually " 
      bondtype = "/1.0"

 
     atomname = str(atom_symbol1)+str(atomidx)

     if ffname == 'drd' : FFtype2, atom_symbol2 = gft.get_FFtype_DREIDING(neighbour_atom)
     if ffname == 'gaff' : 
      FFtype2 =  fftypes[atomid2]  
      atom_symbol2 = symbs[atomid2]

 
 #    print atomNo2, bond_order, bondtype  
     neighbor_atom_string = str(atom_symbol2)+str(atomidx2)+str(bondtype)  
 
 #    print neighbor_atom_string
     na_string_list.append(neighbor_atom_string)
     na_index_list.append(atomidx2)
 
 # sort out the fftypes - this is the part where we can substitute FF typing syntax   
 # myFFtype = FF_types[atomid]
 # if myFFtype == "H_n": myFFtype = "H_ "
 # if myFFtype == "H_HB": myFFtype = "H___A "
 
 # we need to properly sort the na string prior to printing 
  na_string_list=np.array(na_string_list)
  na_index_list=np.array(na_index_list)
 
 # tricky way to sort out the array against the other
  sorted_index_list=np.sort(na_index_list)
  sorted_string_list = []
  print sorted_index_list 
  for i in range(np.size(na_string_list)):
   xx = na_string_list[na_index_list  ==  sorted_index_list[i]] 
   sorted_string_list.append(xx.tolist()) 
  new_list = np.hstack(sorted_string_list)  # notice the use of np.hstack here to collapse the array before printing out to the string 
  
 # print na_string_list
  delimiter = ' '
  na_string = delimiter.join(new_list)    
 
 # sort out the atomprelabel for the mdf 
  resname = "CORE"
  resno = str(mollabel)
  atomprelabel = resname+"_"+resno+":"+atomname 
  charge_group = "?" 
  occupancy = 1.000000
  uiso = 0.000000 
 
  mdffile_handle.write("%s  %s  %s  %s %3i %3i %12f  0  0  8 %5f %5f %s \n" %(atomprelabel, str(atom_symbol1), FFtype1, charge_group, isotope, formal_charge, charge, occupancy, uiso, na_string))
  carfile_handle.write("%s   %12f  %12f  %12f  %s  %s  %s  %s  %12f \n" %(atomname, xpos, ypos, zpos, resname, resno, FFtype1, str(atom_symbol1), charge))
 

def main():

 from sys import argv
 script, rootname, ffname, run_amber = argv

 print "run by typing \n" 
 print "python mol2car.py rootname ffname run_amber \n"

 if not (ffname == 'drd'  or  ffname == 'gaff') : 
  print "please choose appropriote force field\n"
  exit()


 mol2file = rootname+".mol2"
 gaff_mol2file = rootname+"_gaff.mol2" # this is the file output be antechamber where the gaff types are stored 
 carfile = rootname+".car"
 mdffile = rootname+".mdf"
 print mol2file

# run antechamber here
 run_amber = int(run_amber) 
 if (run_amber == 1) :  
 # in the mac version we dont need to specify user charges  
 # amber_cmd = 'antechamber -fi mol2 -fo mol2 -i '+mol2file+' -o '+gaff_mol2file+' -c user' 
  amber_cmd = 'antechamber -fi mol2 -fo mol2 -i '+mol2file+' -o '+gaff_mol2file  
  os.system(amber_cmd)

# in a seperate function lets read in the .mol2 file and get any unit cell data as well as the number of molecules  
 cellparams,nmol = rdat.getheader_data_from_mol2file(mol2file)

 localtime = time.asctime( time.localtime(time.time()) )


# print out the header for our car file 
 carfile_handle = open(carfile, 'w')


####################################################
# print out the header for the .car file 
#

 carfile_handle.write("!BIOSYM archive 3\n")
 if cellparams[0] == 'none': 
  carfile_handle.write("PBC=OFF\n")
 else: carfile_handle.write("PBC=ON\n")
 carfile_handle.write("Created by ERIC CHAN version 1.0.0 - sorry sir but there was no other way... \n")
 carfile_handle.write("!DATE : %s \n" % (localtime))
 if not cellparams[0] == 'none': 
  carfile_handle.write("PBC   %7f   %7f   %7f   %7f   %7f   %7f   (P1) \n" % (float(cellparams[0]), float(cellparams[1]), float(cellparams[2]), float(cellparams[3]), float(cellparams[4]), float(cellparams[5])))



####################################################
# print out the header for the .mdf file 
#

 mdffile_handle = open(mdffile, 'w')

 mdfheadertxt = """!BIOSYM molecular_data 4
 
!ERIC CHAN's OBabel Generated MDF file
 
#topology

@column 1 element
@column 2 atom_type
@column 3 charge_group
@column 4 isotope
@column 5 formal_charge
@column 6 charge
@column 7 switching_atom
@column 8 oop_flag
@column 9 chirality_flag
@column 10 occupancy
@column 11 xray_temp_factor
@column 12 connections

"""

 mdffile_handle.write(mdfheadertxt)

####################################################

# this is a way to loop through multiple molecules in the .arc file  
# just keep re-initialising the molecule object by using the obConv.Read() function   
# example is below.
# it should be most reliable if you are using a for loop  

# mol = openbabel.OBMol()
# obConv.ReadFile(mol,arcfile)
# print mol.GetMolWt()
# mol = openbabel.OBMol()
# obConv.Read(mol)
# print mol.GetMolWt()


 # now we use babel to get the data we want
 obConv = openbabel.OBConversion()
 obConv.SetInFormat("mol2")
 # obConversion.OpenInAndOutFiles("ethane.mol2","ethane.cif")

 for molnum in range(nmol):
  mol = openbabel.OBMol()          # reinitialise the molecule object 
  if molnum == 0:
   obConv.ReadFile(mol, mol2file)   # open babel reads the arc file the first time  
  else:
   obConv.Read(mol)                # it does this each consecutive time  
  print molnum
  print mol.GetMolWt()


  read_write_mol_data(mol,carfile_handle,mdffile_handle,molnum,ffname,gaff_mol2file)
  carfile_handle.write("end\n")


 mdffile_handle.write("\n#end\n")
 carfile_handle.write("end\n")

 carfile_handle.close()
 mdffile_handle.close()


if __name__ == "__main__":
    main()


####################################################
# in case you need to output as a different format 

# obConv.SetOutFormat("mol2")
# obConv.WriteFile(mol, 'ethane.mol2')
# obConv.SetOutFormat("pdb")
# obConv.WriteFile(mol, 'ethane.pdb')
 
###############################################################
# 
#
