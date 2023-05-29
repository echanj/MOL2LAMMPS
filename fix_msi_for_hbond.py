#!/opt/local/bin/python2.7

################################################
#
# this is the script that is adapted ffrom the previous fix_lmp_hbond_input.py that will overwrite 
# lammps.in and .data files that are output from fix_msi2lmp
#
#
# sept 2015 -eric chan  - dev version
#
# create modified lammps input files for hydrogen bonding  
#
# These software tools are currently
# Licensed under the Academic Free License version 3.0
# The authors hold no responsibility whatsoever for any software 
# or hardware related damages that can occur as a result of the use of these codes.
#
# See https://opensource.org/licenses/AFL-3.0 for the details.

#
#########################################################
#
#

import re
import sys
import os
import numpy as np

###########################################

def get_pair_coeffs_from_datafile(datafile,ntyp):

 n = 0
 eps = []
 sigma = []
 atomtyping = []
 with open(datafile) as file:
        for line in file:
         if re.search('Pair Coeffs', line): break
        line = file.next()  # this is the handy trick to jump the next line
        for line in file:
         n += 1
         # print line
         words = line.split()[4]
         print words
         eps.append(float(line.split()[1]))
         sigma.append(float(line.split()[2]))
         atomtyping.append(line.split()[4])
         if n >= ntyp: break  

 return eps,sigma,atomtyping

def get_num_bond_types_from_datafile(datafile):

 with open(datafile) as file:
        for line in file:
         if re.search('bond types', line): 
            words = line.split()
            numbondtypes = int(words[0])

 return numbondtypes

def get_num_atom_types_from_datafile(datafile):

 with open(datafile) as file:
        for line in file:
         if re.search('atom types', line): 
            words = line.split()
            numatomtypes = int(words[0])

 return numatomtypes

def get_hbond_donor_acceptor(dataname,atomtyping,outfh):

 ntyp = get_num_bond_types_from_datafile(dataname)

 # get possible acceptors
 acceptor_ids = []
 for a in range(np.size(atomtyping)):
  if re.match('H___A',atomtyping[a]): hbonding_hyd_type_num = a+1
  alabel = re.sub('_\w+','',atomtyping[a])  
  if alabel == 'O' or alabel == 'N' or alabel == 'S' or alabel == 'Cl' or alabel == 'F'  :    # do add to the list of possible H bonding aceptors when needed 
   acceptor_ids.append(a+1)                    # store the acceptor id  
 
 print hbonding_hyd_type_num 
 print acceptor_ids 

 with open(dataname) as file:
        for line in file:
         if re.search('Bond Coeffs', line):  
          line = file.next()  # this is the handy trick to jump the next line
          for i in range(ntyp): 
           line = file.next()  # this is the handy trick to jump the next line
           words=line.split()
           print line            
           if re.search('H___A', line): 
            print "h-bond donor identified" 
            bond_string = words[4]
            donorlabel = re.sub('H___A','',bond_string)   # w is for word characters, W is for non-word chracters 
            donorlabel = re.sub('-','',donorlabel)   # w is for word characters, W is for non-word chracters 

         # get the atom type id for the donor 
            for a in range(np.size(atomtyping)): 
             if re.match(donorlabel,atomtyping[a]): hbonding_donor_type_num = a+1
            print donorlabel, hbonding_donor_type_num

         # create lines for every possible donor acceptors combination where the donor number <= aceptor number  

            for ai in acceptor_ids:
             if ai >= hbonding_donor_type_num:
              print "pair_coeff   %i  %i  hbond/dreiding/lj  %i  i   4.0  2.7500   4\n" %(hbonding_donor_type_num,ai,hbonding_hyd_type_num ) 
              outfh.write("pair_coeff   %i  %i  hbond/dreiding/lj  %i  i   4.0  2.7500   4\n" %(hbonding_donor_type_num,ai,hbonding_hyd_type_num )) 
             else:
              print "pair_coeff   %i  %i  hbond/dreiding/lj  %i  j   4.0  2.7500   4\n" %(ai,hbonding_donor_type_num,hbonding_hyd_type_num ) 
              outfh.write("pair_coeff   %i  %i  hbond/dreiding/lj  %i  j   4.0  2.7500   4\n" %(ai,hbonding_donor_type_num,hbonding_hyd_type_num )) 

          break


def fix_lammps_input(rootname):

 dataname = rootname+".data"
 newdataname = rootname+"_hbond.data"
 inpname = rootname+".in"
 outname = rootname+"_hbond.in"

 numatomtypes = get_num_atom_types_from_datafile(dataname)
 eps,sigma,atomtyping = get_pair_coeffs_from_datafile(dataname,numatomtypes)

 eps=np.array(eps)
 sigma=np.array(sigma)
 
 outfh = open(outname, 'wb')

 with open(inpname) as file:
        for line in file:
         if not re.search('pair_style', line): 
          if re.search('boundary', line): 
           outfh.write("boundary p p p \n")
           outfh.write("pair_style      hybrid/overlay hbond/dreiding/lj 4 6 6.5 90 lj/charmm/coul/long 8.5000  12.5 \n")
           outfh.write("kspace_style pppm 0.0001 \n")

          elif re.search('# read_data', line): 
           re.sub('read_data','',line)   # w is for word characters, W is for non-word chracters 
          elif re.search('read_data', line): 
            outfh.write("read_data "+newdataname+"\n\n")
            X = np.arange(numatomtypes)
            for i in X:
             for j in X:
               if(i<=j):
                 atom_i = i+1
                 atom_j = j+1
                 epsilon_ij = np.sqrt(eps[i] * eps[j])
                 sigma_ij = (sigma[i] + sigma[j]) / 2.0 
                 outfh.write("pair_coeff   %i  %i  lj/charmm/coul/long   %.6f   %.6f" %(atom_i,atom_j, epsilon_ij, sigma_ij)) 
                 outfh.write("\n")

# this is the part where you need to figure out a general way to tell lammps what acceptors and donors there are  
# should be able to automate this using the bonding profile but for now can be adhoc 
            outfh.write("# please edit the below h-bonding information for each system \n")
            outfh.write("# or edit this part in the fix_msi_hbond_script directly \n")
           #  outfh.write("pair_coeff   2  3  hbond/dreiding/lj  4  j   4.0  2.7500   4\n\n")

# this is where we construct the donor aceptor line
            get_hbond_donor_acceptor(dataname,atomtyping,outfh)
            outfh.write("\n # end_pair_coeff #  \n\n")
            outfh.write("# uncomment the below lines if you require output of hbond info\n\n")  
            outfh.write("# group          mobile union all\n")
            outfh.write("# compute   hb all pair hbond/dreiding/lj\n")
            outfh.write("# variable    C_hbond equal c_hb[1] #number hbonds\n")
            outfh.write("# variable    E_hbond equal c_hb[2] #hbond energy\n\n")
            outfh.write("# thermo_style 	custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong v_E_hbond v_C_hbond press vol\n")
            outfh.write("# thermo_modify	line multi format float %14.6f\n\n")

          else:
           outfh.write(line)
       # line = file.next()  # this is the handy trick to jump the next line

 outfh.close()

def fix_data_input(rootname):

 dataname = rootname+".data"
 newdataname = rootname+"_hbond.data"

 numatomtypes = get_num_atom_types_from_datafile(dataname)
 
 outfh = open(newdataname, 'wb')
  
 with open(dataname) as file:
        for line in file:
         if re.search('Pair Coeffs', line):
          for i in range(numatomtypes+3):             # the +3 here is to get rid of unwanted lines
           line = file.next()  # this is the handy trick to jump the next line
         outfh.write(line)

#########################################
####################################

def main():

 from sys import argv
 script, rootname = argv

 fix_lammps_input(rootname)
 fix_data_input(rootname)

 

if __name__ == "__main__":
    main()

#####################################
# spare code 
##############################################################

 # timesteps,com = get_time_ave_data(fname) 
 # ftimesteps,natoms,fx,fy,fz = get_force_data(fname) 

 # atom_pos1 = get_atom_pos_from_datafile(fname,1) 
 # atom_pos2 = get_atom_pos_from_datafile(fname,2) 

#  dh_atoms_index = [3,1,2,6]

# ftimesteps,dangles = get_dihedral_data(fname,dh_atoms_index) 

 # atompos = get_atom_pos_from_last_frame(fname,1)

# print atom_pos1
# print atom_pos2

 # energies = get_energies(fname) 

#  print energies


#          if line.strip() == '_atom_site_label':  # Or whatever test is needed
#           break
#        for line in file:
#          line = line.replace("\n", " ") # get rid of newline 
#          match_underscore = re.search('^\_', line) 
#          if not match_underscore :   # this is the block we want 
#           words = line.split()
#           type_symbol = str(words[1])   # in most .cif formats the symbol comes after the initial label 
#           if type_symbol:
#            new_label = type_symbol+str(atomnum)+" "
#            atomnum += 1
#            newline = re.sub('^\w+',new_label,line)   # w is for word characters, W is for non-word chracters 
#            print newline               
#          else: print line


        #  if not match_whitespace: print line
        #  else: break

        #  if line.strip() == '\#END':  
        #   break
        #  if line.strip() == 'loop_':  
        #   break
# print the rest of the cif out 
       # for line in cif:
       #   line = line.replace("\n", " ") # get rid of newline 
       #   print line               



