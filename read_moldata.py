#!/usr/bin/python

# run by typing
# python regex.py result_t1.out 

# These software tools are currently
# Licensed under the Academic Free License version 3.0
# The authors hold no responsibility whatsoever for any software 
# or hardware related damages that can occur as a result of the use of these codes.
#
# See https://opensource.org/licenses/AFL-3.0 for the details.


import re
import sys
import os


# This is a standard way to loop through a block of data  
# using 'break'  
# atomnum = 1;

def getheader_data_from_mol2file(fname):

 endcounter = 0
 with open(fname) as file:
        for line in file:
         molstart = re.search('@<TRIPOS>MOLECULE', line)
         pbc_flag = re.search('@<TRIPOS>CRYSIN', line)
      #   create_flag = re.search('Created', line)
         if molstart: endcounter += 1
         if (pbc_flag):  
          line = file.next() 
          line = line.replace("\n", " ") # get rid of newline 
          cellparams = line.split()

# here nmol should not be -1 because there is always an extra end in an arcfile
 nmol = endcounter          

 try: 
  cellparams
 except NameError:
  # not assigned
  cellparams = ["none"]
#  else:
#   print "FFtype is defined."

 return cellparams,nmol

def getheader_data_from_arcfile(fname):

 endcounter = 0
 headerlines = []
 with open(fname) as file:
        for line in file:
         comment = re.search('!', line)
         pbc_flag = re.search('PBC', line)
         create_flag = re.search('Created', line)
         if re.search('^\end', line): endcounter += 1
         if (comment or pbc_flag or create_flag):  
          headerlines.append(line)

 nmol = endcounter-1          
 return headerlines,nmol


def readmol_from_arcfile(fname,molnumber):

 endcounter = 0
 atomnames = []
 resnames = []
 resnos = []
 FF_types = []
 atomtypes = []

 with open(fname) as file:
        for line in file:
        # print line 
        # if re.search('^\s\w', line):
        # if re.search('!', line): line = file.next()  # this is the handy trick to jump the next line
         comment = re.search('!', line)
         pbc_flag = re.search('PBC', line)
         create_flag = re.search('Created', line)
         if re.search('^\end', line): endcounter += 1
         if not (comment or pbc_flag or create_flag):  
          if endcounter == molnumber: 
           if not re.search('^\end', line):
            words = line.split()
            atomnames.append(str(words[0]))               
            resnames.append(str(words[4]))               
            resnos.append(str(words[5]))               
            FF_types.append(str(words[6]))               
            atomtypes.append(str(words[7]))               

 return atomnames,resnames,resnos,FF_types,atomtypes


def get_fftype_data_from_mol2file(fname,molnum):

 fftypes = []
 symbs = []
 charges = []
 molcount = 0 
 with open(fname) as file:
     for line in file:
      molstart = re.search('@<TRIPOS>MOLECULE', line)
      if molstart: molcount += 1
      if molcount == molnum: 

       for line in file:
        if re.search('@<TRIPOS>ATOM', line): break
       for line in file:
        if re.search('@<TRIPOS>BOND', line): break
        fty = line.split()[5]
        fftypes.append(fty)  
        charge = line.split()[8]
        charges.append(charge)  
        atom_symb = line.split()[1]
        atom_symb = re.sub("\d","",atom_symb)
        symbs.append(atom_symb)    
       break
      
 return fftypes,symbs,charges


def main():

 from sys import argv
 script, fname = argv

# headerlines, nmol = getheader_data_from_arcfile(fname)

 cellparams,nmol = getheader_data_from_mol2file(fname)
 fftypes,symbs = get_fftype_data_from_mol2file(fname,1)

 print cellparams
 print nmol 
 print fftypes
 print symbs

#
# atomnames,resnames,resnos,FF_types,atomtypes = readmol_from_arcfile(fname,0) 
# print atomnames,resnames,resnos,FF_types,atomtypes 
#

if __name__ == "__main__":
    main()




