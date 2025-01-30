#!/opt/apes/miniconda3/bin/python

import numpy as np
import os

##############################
#          mbe library       #
##############################
# contains stuff
##############################
#
# A) Scripts used to generate n-mers in many-body expansion
#
# B) Scripts to read and process output files with energies 
#
# C) Helper functions, these can need ASE toolkit
#
##############################


#############################
#   get_distance(nm1,nm2)   #
#############################
#
# calculate average distance between two molecules
# given by their xyz structures in Angstrom
#
#############################
# nm1 are the xyz coordinates of the first molecule
# nm2 are the xyz coordinates of the second molecule
#############################
def get_distance(nm1,nm2):
   dist=0.0
   for iat1 in range(nm1.shape[0]):
      for iat2 in range(nm2.shape[0]):
         dist+=np.sqrt((nm1[iat1,0]-nm2[iat2,0])**2+(nm1[iat1,1]-nm2[iat2,1])**2+(nm1[iat1,2]-nm2[iat2,2])**2)
   #print dist/nm1.shape[0]/nm2.shape[0]
   return dist/float(nm1.shape[0])/float(nm2.shape[0])

#####################################################
#   get_distance_from_vec(nm1,vec1,nm2,vec2,cell)   #
#####################################################
#
# calculate a distance between two molecules which reside in unit cells
# defined by the vec lists
#
###############
# nm1 and nm2 are the coordinates of two unit cell monomers stored as [natoms,xyz]
# vec1 and vec2 are shifts [0 1 0]
# cell is a 3x3 matrix of the unit cell
###############
def get_distance_from_vec(nm1,vec1,nm2,vec2,cell):
   nm1_act=nm1.copy()
   nm2_act=nm2.copy()
   #get new structures in the correct positions
   for iat in range(nm1.shape[0]):
      for ix in range(3):
         nm1_act[iat,ix]=nm1[iat,ix]+vec1[0]*cell[0,ix]+vec1[1]*cell[1,ix]+vec1[2]*cell[2,ix]
         #print(vec1[0]*cell[0,ix])
         #print(vec1[1]*cell[1,ix])
         #print(vec1[2]*cell[2,ix])
   for iat in range(nm2.shape[0]):
      for ix in range(3):
         nm2_act[iat,ix]=nm2[iat,ix]+vec2[0]*cell[0,ix]+vec2[1]*cell[1,ix]+vec2[2]*cell[2,ix]
   #print(vec1)
   #print(nm1_act[0,:])
   #print(vec2)
   #print(nm2_act[0,:])
   #dist=0.0
   #for iat1 in range(nm1.shape[0]):
   #   for iat2 in range(nm2.shape[0]):
   #      dist+=np.sqrt((nm1_act[iat1,0]-nm2_act[iat2,0])**2+(nm1_act[iat1,1]-nm2_act[iat2,1])**2+(nm1_act[iat1,2]-nm2_act[iat2,2])**2)
   ##print dist/nm1.shape[0]/nm2.shape[0]
   #call get_distance to get the distances between the two molecules
   dist=get_distance(nm1_act,nm2_act)
   return dist

##########################
#   get_vec(idn,mlist)   #
##########################
#
# Returns the vec
# TODO expand
def get_vec(idn,mlist):
   return mlist[idn,2:6]

##############################
#   get_vec_dir(idx,mlist)   #
##############################
#
# returns a string containing the directory name
# corresponding to the index of molecule according 
# to monomer list mlist
# we want to get this format shell0_X0_Y0_Z0_mono1_fix1
##############################
#
# id0 index of molecule in the list of monomers
# mlist list of monomers
#
def get_vec_dir(idx,mlist):
#                mlist[idx,0]=idx
#                mlist[idx,1]=shell
#                mlist[idx,2]=x
#                mlist[idx,3]=y
#                mlist[idx,4]=z
#                mlist[idx,5]=mono
   string='shell'+str(mlist[idx,1])
   string+='_X'+str(mlist[idx,2])
   string+='_Y'+str(mlist[idx,3])
   string+='_Z'+str(mlist[idx,4])
   string+='_mono'+str(mlist[idx,5])
   string+='_fix1/'
   return string

#get the index into mlist knowing the shell, XYZ and monomer data stored in a list as ints
def get_idx_vecall(mlist,vecall):
   for i in range(mlist.shape[0]):
     if mlist[i,1]==vecall[0]:
       if mlist[i,2]==vecall[1]:
         if mlist[i,3]==vecall[2]:
           if mlist[i,4]==vecall[3]:
              if mlist[i,5]==vecall[4]:
                return i

# get the numbers from the monomer string
def get_idx_str(mlist,string):
   #split shell0_X0_Y0_Z0_mono1_fix1 
   spl=string.split('_')
   #print(spl)
   #store all including shell and monomer info
   #print(spl)
   vecall=[int(spl[0][5:]),int(spl[1][1:]),int(spl[2][1:]),int(spl[3][1:]),int(spl[4][4:])]
   #print(vecall)
   idx=get_idx_vecall(mlist,vecall)
   return idx


#####################################
#   mono_shift_vec(nm1,vec1,cell)   #
#####################################
#
# shift a monomer geometry according to the vector
#
###############
# nm1 are the coordinates of unit cell monomer stored as [natoms,xyz]
# vec1 is shift [0 1 0]
# cell is a 3x3 matrix of the unit cell
###############
def mono_shift_vec(nm1,vec1,cell):
   nm1_act=nm1.copy()
   for iat in range(nm1.shape[0]):
      for ix in range(3):
         nm1_act[iat,ix]=nm1[iat,ix]+vec1[0]*cell[0,ix]+vec1[1]*cell[1,ix]+vec1[2]*cell[2,ix]
   return nm1_act

###################################################
#   get_mono_struc_from_idx(idx,mlist,coord,cell)   #
###################################################
#
# get a monomer structure when we know its index into mlist 
#
###################################################
#
# idx index into mlist
# mlist list of monomers 
# coord coordinates of monomers in the unit cell
# cell unit cell vectors
def get_mono_struc_from_idx(idx,mlist,coord,cell):
   nm1=coord[mlist[idx,5],:,:].copy()
   vec1=mlist[idx,2:5]
   for iat in range(nm1.shape[0]):
      for ix in range(3):
         nm1[iat,ix]=nm1[iat,ix]+vec1[0]*cell[0,ix]+vec1[1]*cell[1,ix]+vec1[2]*cell[2,ix]
   return nm1


##############################
#   get_dir_vec(mdir)   #
##############################
#
# on input reads the monomer identifier in format shell0_X0_Y0_Z0_mono1_fix1
# returns a vector containing the xyz shift and monomer index
##############################
#
# mdir is the directory
#
def get_dir_vec(mdir):
   spl=mdir.split('_')
   x=int(spl[1][1:])
   y=int(spl[2][1:])
   z=int(spl[3][1:])
   mono=int(spl[4][4:])
   vec=[x,y,z,mono]
   return vec



#########################################################
#   get_distance_from_idx(idx1,idx2,mlist,coord,cell)   #
#########################################################
#
# Calculates the distance between two molecules
# The molecules are defined by their indices into the list of monomers
# The routine simply finds the shifts and monomer numbers in the list
# and calls get_distance_from_vec
####################################################
#
# idx1 is the index of the first monomer
# idx2 is the index of the second monomer
# mlist is the list between monomer ID and xyz shift and monomer number in unit cell
# coord are the coordinates of monomers in the unit cell
####################################################
def get_distance_from_idx(idx1,idx2,mlist,coord,cell):
#                mlist[idx,0]=idx
#                mlist[idx,1]=shell
#                mlist[idx,2]=x
#                mlist[idx,3]=y
#                mlist[idx,4]=z
#                mlist[idx,5]=mono
   #get the xyz coordinates of monomers 1 and 2
   vec1=mlist[idx1,2:5]
   vec2=mlist[idx2,2:5]
   nm1=coord[mlist[idx1,5],:,:]
   nm2=coord[mlist[idx2,5],:,:]
   dist=get_distance_from_vec(nm1,vec1,nm2,vec2,cell)
   return dist


#######################
#   read_cell(path)   #
#######################
#
#  read the unit cell
#  the cell object (numpy array) is created by the function
#
#######################
#
# path is the (relative) directory path to place where 
# the file "cell" is located
#######################
def read_cell(path):
   cell=np.zeros((3,3),dtype=np.float64)
   if os.path.isfile(path) == True:
      #print('file is')
      with open(path) as cell_file:
         cell_lin=cell_file.readlines()
         idx=0
         #read three lines from the unit cell and store them in the cell array
         #in case we have one trailing empty line
         for line in cell_lin:
            spl=line.split()
            #print (spl)
            for j in range(3):
               cell[idx,j]=float(spl[j])
            idx+=1
            #three lines ought be enough for everyone
            if idx==3:
              break
      print('Unit cell found and read')
   else:
      print('No unit cell found')
  
   return cell


###############################
#   read_monomers_uni(path)   #
###############################
#
# Reads the monomer structures and stores them in an array
# The array is created on the fly and is returned to the
# calling programme
# At the beginning we don't know how many monomers there are
# and we expect them to be called as 0.tail ... 
# the tail could be a function parameter
# This works for monomers which are identical, hence the _uni
# suffix
#
#####################
#
# str path: location of the monomers
# ###################
def read_monomers_uni(path):
   #OK, 128 is pretty arbitrary
   for nmono in range(128):
      #test how many monomers there are
      if os.path.isfile(path+'/'+str(nmono)+'.tail') != True:
         break
   print('Found '+str(nmono)+' monomers in '+path)

   #count how many items there are
   natoms = 0
   for line in open(path+'/0.tail'):
      natoms += 1
   
   #create array with monomer coordinates
   #mind the C-style indexing!
   coord=np.zeros((nmono,natoms,3))
   elems=[]

   #loop over all monomers 
   for idx_mono in range(nmono):
      elem=[]
      if os.path.isfile(path+'/'+str(idx_mono)+'.tail') == True:
         #open one monomer file with the structure
         with open(path+'/'+str(idx_mono)+'.tail') as mono_file:
            mono_lin=mono_file.readlines()
            idx=0
            #go through the lines in the structure file and store them 
            for line in mono_lin:
               spl=line.split()
               #print (spl)
               #possible TODO? do we assume that line with four entries is data? really? fcplm
               if len(spl)==4:
                  for j in range(3):
                     coord[idx_mono,idx,j]=float(spl[j+1])
                  elem.append(spl[0])
               idx+=1
         elems.append(elem)
      else:
         print('No monomer structures found, I give up')
         exit()

   #return the coordinates on exit 
   return coord, elems

def get_Z(elem):
   if elem=='H':
      return 1
   elif elem=='C':
      return 6
   elif elem=='N':
      return 7
   elif elem=='O':
      return 8

#################################################
#################################################
#                                               #
#    generate, read, write a list of monomers   #

########################
#    read_list(path)   #
########################
#
# reads the monomer list from a given file
#
########################
#
# str path is the directory path to list file
########################
def read_list(path):
   if os.path.isfile(path) == True:
      #count how many items there are
      nlines = 0
      for line in open(path):
         nlines += 1
      #for each monomer we store its idx, shell, x, y, z, mono
      mlist=np.zeros((nlines,6),dtype=int)
      idx=0
      for line in open(path):
         spl=line.split()
         for item in range(6):
             mlist[idx,item]=int(spl[item])
         idx+=1
      print('List found and read with '+str(nlines)+' items')
   else:
      print('List not found')
   return mlist

########################
#    write_list(path)   #
########################
#
# writes the monomer list to a given file
#
########################
#
# str path is the name of the list file
# mlist is the monomer list
########################
def write_list(path, mlist):
   with open(path,'w+') as f:
      for idx in range(mlist.shape[0]):
         #this is ok, but too many empty spaces
         #f.write(str(mlist[idx,:])[1:-1].lstrip(' ')+'\n')
         #f.write(''.join(map(str,mlist[idx,:])))
         #f.write(np.array2string(mlist[idx,:])[1:-1]+'\n')
         #f.write(' '.join(str(mlist[idx,n]) for n in mlist[idx,0:4]))
         #f.write(*mlist[idx,:].flatten()+'\n')
         #this creates a list from the row, changes the elements to string
         #adds them together, adds to end line and prints it
         f.write(" ".join(map(str,mlist[idx,:].tolist()))+'\n')
   
   #if os.path.isfile(path) == True:
   #   #count how many items there are
   #   nlines = 0
   #   for line in open(path):
   #      nlines += 1
   #   #for each monomer we store its idx, shell, x, y, z, mono
   #   mlist=np.zeros((nlines,6),dtype=int)
   #   idx=0
   #   for line in open(path):
   #      spl=line.split()
   #      for item in range(6):
   #          mlist[idx,item]=int(spl[item])
   #      idx+=1
   #   print('List found and read with '+str(nlines)+' items')
   #else:
   #   print('List not found')
   #return mlist

##################################
#   gen_list(ncell,nmono)   #
###############################################################################
#
# generate a monomer list for a given number of monomers and number of shells 
###############################################################################
#
# ncell is the number of shells
# nmono is the number of monomers in the unit cell
##################################
def gen_list(ncell,nmono):

   #nmono=4
   #nmono is the numer of monomers in a unit cell
   #ncell=6
   #ncell is the number of shells 
   #(shell 0 is the unit cell, shell 1 is around it)

   nitems=nmono*(2*ncell+1)**3
   #for each monomer we store its idx, shell, x, y, z, mono
   mlist=np.zeros((nitems,6),dtype=int)
   #ass_list=np.full((2*ncell+1,2*ncell+1,2*ncell+1,nmono),-1,dtype=int)
   #ass_list=np.zeros((2*ncell+1,2*ncell+1,2*ncell+1,nmono),dtype=int)
   #so if we want x indices (-1,0,1)
   #we have python indices (0,1,2)
   #and need to use a shift of ncell

   idx=0
   for shell in range(ncell+1):
     #print(shell)
     for z in [-shell,shell]:
       for x in range(-shell,shell+1):
          for y in range(-shell,shell+1):
             #the monomers are indexed from zero
             for mono in range(nmono):
                mlist[idx,0]=idx
                mlist[idx,1]=shell
                mlist[idx,2]=x
                mlist[idx,3]=y
                mlist[idx,4]=z
                mlist[idx,5]=mono
                #ass_list[x+ncell,y+ncell,z+ncell,mono]=idx
                #print(' '.join(map(str,mlist[idx-1])))
                idx+=1
       #break from the loop in case shell is zero so that it does not appear twice
       if -shell == shell:
          break
   
     if not shell == 0:
        for x in [-shell,shell]:
          for z in range(-shell+1,shell):
             for y in range(-shell,shell+1):
                for mono in range(nmono):
                   mlist[idx,0]=idx
                   mlist[idx,1]=shell
                   mlist[idx,2]=x
                   mlist[idx,3]=y
                   mlist[idx,4]=z
                   mlist[idx,5]=mono
                   #ass_list[x+ncell,y+ncell,z+ncell,mono]=idx
                   #print(list[idx-1])
                   #print(' '.join(map(str,mlist[idx-1])))
                   idx+=1
          #break from the loop in case shell is zero so that it does not appear twice
          if -shell == shell:
             break
   
        for y in [-shell,shell]:
          for x in range(-shell+1,shell):
             for z in range(-shell+1,shell):
                for mono in range(nmono):
                   mlist[idx,0]=idx
                   mlist[idx,1]=shell
                   mlist[idx,2]=x
                   mlist[idx,3]=y
                   mlist[idx,4]=z
                   mlist[idx,5]=mono
                   #ass_list[x+ncell,y+ncell,z+ncell,mono]=idx
                   #print(list[idx-1])
                   #print(' '.join(map(str,mlist[idx-1])))
                   idx+=1
          #break from the loop in case shell is zero so that it does not appear twice
          if -shell == shell:
             break
   return mlist


#############################################
#############################################
#                                           #
#   generate mers from a list of monomers   #

#######################################################
#   gen_dimers_old(ncell,coord,cell,file_out)   #
#######################################################
#
# create a list of dimers stored in a file
# in contrast to the directory names used for trimers,
# the dimers are stored shellN/XN_YN_ZN_monoN_fix1
#
# superseded by gen_dimers which uses the same convention as 
# gen_trimers and higher routines, this makes it easier
# to work with more reference molecules in the unit cell
#######################################################
#
# ncell is the number of shells to use 
def gen_dimers_old(nref,ncell,coord,cell,file_out):
   #create the list first, for this we need to know the number of monomers
   mlist=gen_list(ncell,coord.shape[0])

   #all we have to do now is turn the list into strings to print
   with open(file_out,'w+') as f:
      for idx in range(mlist.shape[0]):
         if idx==nref:
            continue
         row=mlist[idx,:].tolist()
         row=[str(i) for i in row]
         #row=list(map(str,mlist[idx,:].tolist()))
         #row=mlist[idx,:].tolist()
         #shell2/X0_Y1_Z2_mono0_fix1
         f.write('shell'+row[1]+'/X'+row[2]+'_Y'+row[3]+'_Z'+row[4]+'_mono'+row[5]+'_fix'+str(nref)+'\n')
   
#######################################################
#   gen_dimers(ncell,coord,cell,file_out)   #
#######################################################
#
# create a list of dimers stored in a file
# the same convention is used as for trimers and tetramers
# the dimers are stored shellN_XN_YN_ZN_monoN_fixN
#
#######################################################
#
# ncell is the number of shells to use 
def gen_dimers(nref,ncell,coord,cell,file_out):
   #create the list first, for this we need to know the number of monomers
   mlist=gen_list(ncell,coord.shape[0])

   #all we have to do now is turn the list into strings to print
   with open(file_out,'w+') as f:
      for idx in range(mlist.shape[0]):
         if idx==nref:
            continue
         row=mlist[idx,:].tolist()
         row=[str(i) for i in row]
         #row=list(map(str,mlist[idx,:].tolist()))
         #row=mlist[idx,:].tolist()
         #shell2/X0_Y1_Z2_mono0_fix1
         #f.write('shell'+row[1]+'/X'+row[2]+'_Y'+row[3]+'_Z'+row[4]+'_mono'+row[5]+'_fix'+str(nref)+'\n') 
         #use a directory structure where the indices ordering is preserved
         if nref < idx:  
            sid0=get_vec_dir(nref,mlist)
            sid1=get_vec_dir(idx,mlist)
         else:
            sid1=get_vec_dir(nref,mlist)
            sid0=get_vec_dir(idx,mlist)
         f.write(sid0+sid1+' '+str(get_distance_from_idx(idx,nref,mlist,coord,cell))+'\n')

#################################################
#   gen_tetramers(nmono_max,mlist,coord,cell)   #
#################################################
#
#
#
def gen_tetramers(nmono_max,dist_cut,mlist,coord,cell,file_out):
   dist_min=3.68-0.1
   id0=0
   list_tot=[]
   #with open(file_out,'w+') as f:
   for id1 in range(1,nmono_max):
      print('working on '+str(id1)+' out of '+str(nmono_max))
      dist01=get_distance_from_idx(id0,id1,mlist,coord,cell)
      if dist01>dist_cut-5*dist_min:
         continue
      for id2 in range(id1+1,nmono_max):
         dist02=get_distance_from_idx(id0,id2,mlist,coord,cell)
         dist12=get_distance_from_idx(id1,id2,mlist,coord,cell)
         if dist01+dist02+dist12>dist_cut-3*dist_min:
            continue
         for id3 in range(id2+1,nmono_max):
            dist03=get_distance_from_idx(id0,id3,mlist,coord,cell)
            dist13=get_distance_from_idx(id1,id3,mlist,coord,cell)
            dist23=get_distance_from_idx(id2,id3,mlist,coord,cell)
            tot_dist=dist01+dist02+dist12+dist03+dist13+dist23
            if tot_dist<=dist_cut:
               #print(tot_dist,get_vec(id1,mlist),get_vec(id2,mlist))
               #print (str(id0)+' '+str(id1)+' '+str(id2)+' '+str(id3))
               #f.write(get_vec_dir(id0,mlist)+get_vec_dir(id1,mlist)+get_vec_dir(id2,mlist)+get_vec_dir(id3,mlist)+' '+str(dist01)+' '+str(dist02)+' '+str(dist03)+' '+str(tot_dist)+'\n')
               list_tot.append([id0,id1,id2,id3,dist01,dist02,dist03,tot_dist])
   list_tot.sort(key=lambda x: x[7])
   #print(list_tot) 
   with open(file_out,'w+') as f:
      for tetr in list_tot:
         sid0=get_vec_dir(int(tetr[0]),mlist)
         sid1=get_vec_dir(int(tetr[1]),mlist)
         sid2=get_vec_dir(int(tetr[2]),mlist)
         sid3=get_vec_dir(int(tetr[3]),mlist)
         f.write(sid0+sid1+sid2+sid3+' '+' '.join(map(str,tetr[4:8]))+'\n')

# BIG TODO: make these n-independent

#########################################################################
#   gen_trimers_dist_cache(shells,dist_cut,mlist,coord,cell,file_out)   #
#########################################################################
#
# Generate the list of trimers using a distance criterion
# an ordered list of monomers needs to be put to the routine
# a nested loop over the monomers is performed 0<id1<id2
# trimers are found, stored in a list, and finally ordered 
# according to distance and written to the file
#
# the dimer distances are precalculated here
#
# modification of the tetramer routine
############################################################
#
# ##nmono_max: number of monomers to read from the monomer list mlist
# shells:    largest shell to consider  
# dist_cut:  distance criterion in Angstrom
# mlist:     monomer list as from gen_list(), each item holds: index shell x y z nmono
# coord:     monomer coordinates in the unit cell
# cell:      unit cell
# file_out:  output file to write the list to
############################################################
def gen_trimers_dist_cache(shells,dist_cut,mlist,coord,cell,file_out):

   #number of monomers in the unit cell
   cell_mono=coord.shape[0]

   #The number of monomers that we will consider is defined by the number of shells
   #which is an input parameter. Moreover we have the mlist which was also generated
   #for some number of shells. We need to make sure that these are consistent.
   #Typically, that the list contains all monomers from all the required shells.
   #Find the number of monomers that we will use
   nmono_max=-1
   #First check if the monomer list was generated for the same number of shells as we will consider
   if mlist[-1,1]==shells:
      #if the last monomer is from the largest shell, then use all monomers
      nmono_max=mlist.shape[0]-1
   #if not, check if we have a list generated for less shells than we want to use
   elif mlist[-1,1]<shells:
      print ('Error, monomer list not sufficient, rerun mbe.gen_list with a larger number of shells')
      exit()
   #if we have a longer list, find which is the highest monomer index in the largest shell we want to use
   else:
      #find the largest monomer index from the shell
      #loop over all the monomers in the monomer list
      for idx in range(mlist.shape[0]):
         #we find a monomer from a larger shell
         if mlist[idx,1]>shells:
            nmono_max=idx-1
            break
      #this should not happen if everything is fine
      if nmono_max==-1:
         print ('Can not identify monomer from your shell')   
         #That's one possible reason to fail
         print ('Check if you used the correct number of monomers to generate mlist')
         exi()

   #Now precalculate the distances between the monomers, the array needs to
   # be larger as we need to have distances between monomers in, say, cells +n and -n
   # and also -n and +n
   #dist_ar is the array which stores all the distances
   dist_ar=np.zeros((cell_mono,cell_mono,2*shells+1,4*shells+1,4*shells+1)) 
   #this is just a helper array
   shift_base=[0,0,0] 
   #dist_min is the smallest distance between any dimers, it is used to speed-up the generation of nmers
   dist_min=1000000
   #loop over the monomers in the unit cell
   for m1 in range(cell_mono):
      #first monomer coordinates in the unit cell
      cm1=coord[m1,:,:].copy()
      #loop over another monomer
      for m2 in range(cell_mono):
         #second monomer coordinates in the unit cell
         cm2=coord[m2,:,:].copy()
         # Loop over 
         #difference between the x shift of the monomers, 
         #the largest difference between X shifts for a given number of shells
         #is 2*shells, from -shells to shells
         #the smallest one can be set to zero as we can choose the monomer
         for xsh in range(2*shells+1):
            #for y and z we can't shuffle the monomers, so the difference between cells
            #is from -2*shells to 2*shells, in total 4*shells
            #in the range ysh=2*shells corresponds to the zero difference between the cells
            for ysh in range(4*shells+1):
               for zsh in range(4*shells+1):
                  shift_m2=[xsh,ysh-2*shells,zsh-2*shells]
                  dist_act=get_distance_from_vec(cm1,shift_base,cm2,shift_m2,cell)
                  dist_ar[m1,m2,xsh,ysh,zsh]=dist_act
                  if dist_act<dist_min:
                     dist_min=dist_act
                  
  
   print('Precalculated distances array size ')
   print(dist_ar.shape[:])

   id0=0
   list_tot=[]
   #with open(file_out,'w+') as f:
   for id1 in range(1,nmono_max):
      print('working on '+str(id1)+' out of '+str(nmono_max))
      #dist01=get_distance_from_idx(id0,id1,mlist,coord,cell)
      #get xyz and monomer information
      dat1=mlist[id1,2:6]
      #loop up the distance between the two monomers in the dist_ar array
      #if the xshift is 0 or larger, use the shift directly
      #need to shift the y and z index to larger values
      if dat1[0]>=0:
         dist01=dist_ar[0,dat1[3],dat1[0],dat1[1]+2*shells,dat1[2]+2*shells]
      else: 
         #if the xshift is negative, swap the monomers and shift
         dist01=dist_ar[dat1[3],0,-dat1[0],-dat1[1]+2*shells,-dat1[2]+2*shells]
      #if the current distance plus two smallest are larger than the cut-off
      #we do not need to continue 
      if dist01>dist_cut-2*dist_min:
         continue
      for id2 in range(id1+1,nmono_max):
         #dist02=get_distance_from_idx(id0,id2,mlist,coord,cell)
         dat2=mlist[id2,2:6]
         #if the xshift is 0 or larger, use the shift directly to obtain the precalculated distance
         if dat2[0]>=0:
            dist02=dist_ar[0,dat2[3],dat2[0],dat2[1]+2*shells,dat2[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift when accessing the distance array
            dist02=dist_ar[dat2[3],0,-dat2[0],-dat2[1]+2*shells,-dat2[2]+2*shells]
         #dist12=get_distance_from_idx(id1,id2,mlist,coord,cell)
         #calculate the distance to second monomer
         #if the xshift is same or larger than that of monomer1, use the shift directly
         if dat2[0]>=dat1[0]:
            dist12=dist_ar[dat1[3],dat2[3],dat2[0]-dat1[0],dat2[1]-dat1[1]+2*shells,dat2[2]-dat1[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift
            dist12=dist_ar[dat2[3],dat1[3],dat1[0]-dat2[0],dat1[1]-dat2[1]+2*shells,dat1[2]-dat2[2]+2*shells]
          
         tot_dist=dist01+dist02+dist12
         if tot_dist<=dist_cut:
            #print(tot_dist,get_vec(id1,mlist),get_vec(id2,mlist))
            #print (str(id0)+' '+str(id1)+' '+str(id2)+' '+str(id3))
            #f.write(get_vec_dir(id0,mlist)+get_vec_dir(id1,mlist)+get_vec_dir(id2,mlist)+get_vec_dir(id3,mlist)+' '+str(dist01)+' '+str(dist02)+' '+str(dist03)+' '+str(tot_dist)+'\n')
            list_tot.append([id0,id1,id2,dist01,dist02,dist12,tot_dist])
   #sort the list
   list_tot.sort(key=lambda x: x[6])
   #print(list_tot) 
   #write the list to file
   with open(file_out,'w+') as f:
      for trim in list_tot:
         sid0=get_vec_dir(int(trim[0]),mlist)
         sid1=get_vec_dir(int(trim[1]),mlist)
         sid2=get_vec_dir(int(trim[2]),mlist)
         f.write(sid0+sid1+sid2+' '+' '.join(map(str,trim[3:7]))+'\n')

##############################################################################
#   gen_trimers_dist_cache_cell(shells,dist_cut,mlist,coord,cell,file_out)   #
##############################################################################
#
# Generate the list of trimers using a distance criterion
# an ordered list of monomers needs to be put to the routine
# a nested loop over the monomers is performed 0<id1<id2
# for the last index all the monomers in the unit cell are considered
# if one of them satisfies the cut-off criterion, all are used
# trimers are found, stored in a list, and written to the output file
# ordered according to distance and not ordered
#
# the dimer distances are precalculated here
#
# modification of the tetramer routine
############################################################
#
# ##nmono_max: number of monomers to read from the monomer list mlist
# shells:    largest shell to consider  
# dist_cut:  distance criterion in Angstrom
# mlist:     monomer list as from gen_list(), each item holds: index shell x y z nmono
# coord:     monomer coordinates in the unit cell
# cell:      unit cell
# file_out:  output file to write the list to
############################################################
def gen_trimers_dist_cache_cell(shells,dist_cut,mlist,coord,cell,file_out):

   #number of monomers in the unit cell
   cell_mono=coord.shape[0]

   #The number of monomers that we will consider is defined by the number of shells
   #which is an input parameter. Moreover we have the mlist which was also generated
   #for some number of shells. We need to make sure that these are consistent.
   #Typically, that the list contains all monomers from all the required shells.
   #Find the number of monomers that we will use
   nmono_max=-1
   #First check if the monomer list was generated for the same number of shells as we will consider
   if mlist[-1,1]==shells:
      #if the last monomer is from the largest shell, then use all monomers
      nmono_max=mlist.shape[0]-1
   #if not, check if we have a list generated for less shells than we want to use
   elif mlist[-1,1]<shells:
      print ('Error, monomer list not sufficient, rerun mbe.gen_list with a larger number of shells')
      exit()
   #if we have a longer list, find which is the highest monomer index in the largest shell we want to use
   else:
      #find the largest monomer index from the shell
      #loop over all the monomers in the monomer list
      for idx in range(mlist.shape[0]):
         #we find a monomer from a larger shell
         if mlist[idx,1]>shells:
            nmono_max=idx-1
            break
      #this should not happen if everything is fine
      if nmono_max==-1:
         print ('Can not identify monomer from your shell')   
         #That's one possible reason to fail
         print ('Check if you used the correct number of monomers to generate mlist')
         exi()

   #Now precalculate the distances between the monomers, the array needs to
   # be larger as we need to have distances between monomers in, say, cells +n and -n
   # and also -n and +n
   #dist_ar is the array which stores all the distances
   dist_ar=np.zeros((cell_mono,cell_mono,2*shells+1,4*shells+1,4*shells+1)) 
   #this is just a helper array
   shift_base=[0,0,0] 
   #dist_min is the smallest distance between any dimers, it is used to speed-up the generation of nmers
   dist_min=1000000
   #loop over the monomers in the unit cell
   for m1 in range(cell_mono):
      #first monomer coordinates in the unit cell
      cm1=coord[m1,:,:].copy()
      #loop over another monomer
      for m2 in range(cell_mono):
         #second monomer coordinates in the unit cell
         cm2=coord[m2,:,:].copy()
         # Loop over 
         #difference between the x shift of the monomers, 
         #the largest difference between X shifts for a given number of shells
         #is 2*shells, from -shells to shells
         #the smallest one can be set to zero as we can choose the monomer
         for xsh in range(2*shells+1):
            #for y and z we can't shuffle the monomers, so the difference between cells
            #is from -2*shells to 2*shells, in total 4*shells
            #in the range ysh=2*shells corresponds to the zero difference between the cells
            for ysh in range(4*shells+1):
               for zsh in range(4*shells+1):
                  shift_m2=[xsh,ysh-2*shells,zsh-2*shells]
                  dist_act=get_distance_from_vec(cm1,shift_base,cm2,shift_m2,cell)
                  dist_ar[m1,m2,xsh,ysh,zsh]=dist_act
                  if dist_act<dist_min:
                     dist_min=dist_act
  
   print('Precalculated distances array size ')
   print(dist_ar.shape[:])

   id0=0
   list_tot=[]
   #with open(file_out,'w+') as f:
   for id1 in range(1,nmono_max):
      print('working on '+str(id1)+' out of '+str(nmono_max))
      #dist01=get_distance_from_idx(id0,id1,mlist,coord,cell)
      #get xyz and monomer information
      dat1=mlist[id1,2:6]
      #loop up the distance between the two monomers in the dist_ar array
      #if the xshift is 0 or larger, use the shift directly
      #need to shift the y and z index to larger values
      if dat1[0]>=0:
         dist01=dist_ar[0,dat1[3],dat1[0],dat1[1]+2*shells,dat1[2]+2*shells]
      else: 
         #if the xshift is negative, swap the monomers and shift
         dist01=dist_ar[dat1[3],0,-dat1[0],-dat1[1]+2*shells,-dat1[2]+2*shells]
      #if the current distance plus two smallest are larger than the cut-off
      #we do not need to continue 
      if dist01>dist_cut-2*dist_min:
         continue

      #make a loop over the last monomer in the trimer
      #we check if the monomer belongs to the same cell as the previous one
      #if not, we print all the trimers from the previous cell if there
      #was one with small distance
      #temporary list with trimers in one cell
      list_tmp=[]
      trim_out=False
      dat2_tmp=[-100,-100,-100,-100]
      for id2 in range(id1+1,nmono_max):
         #dist02=get_distance_from_idx(id0,id2,mlist,coord,cell)
         #get current monomer info
         #the numbers should be x, y, z coordinates of cell and monomer index
         dat2=mlist[id2,2:6]

         #check if we have moved to a new cell
         if (dat2[0]!=dat2_tmp[0]) or (dat2[1]!=dat2_tmp[1]) or (dat2[2]!=dat2_tmp[2]):
            #check if some of the previous trimers is below cut-off
            if trim_out==True:
                print(len(list_tmp),' trimers added')
                #add the trimers to the global list
                for elem in list_tmp:
                   list_tot.append(elem)
            #start making a new list for this cell
            list_tmp=[]
            trim_out=False
            dat2_tmp=dat2

         #if the xshift is 0 or larger, use the shift directly to obtain the precalculated distance
         if dat2[0]>=0:
            dist02=dist_ar[0,dat2[3],dat2[0],dat2[1]+2*shells,dat2[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift when accessing the distance array
            dist02=dist_ar[dat2[3],0,-dat2[0],-dat2[1]+2*shells,-dat2[2]+2*shells]
         #dist12=get_distance_from_idx(id1,id2,mlist,coord,cell)
         #calculate the distance to second monomer
         #if the xshift is same or larger than that of monomer1, use the shift directly
         if dat2[0]>=dat1[0]:
            dist12=dist_ar[dat1[3],dat2[3],dat2[0]-dat1[0],dat2[1]-dat1[1]+2*shells,dat2[2]-dat1[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift
            dist12=dist_ar[dat2[3],dat1[3],dat1[0]-dat2[0],dat1[1]-dat2[1]+2*shells,dat1[2]-dat2[2]+2*shells]

         #add the trimer to a temporary list
         tot_dist=dist01+dist02+dist12
         list_tmp.append([id0,id1,id2,dist01,dist02,dist12,tot_dist])
         if tot_dist<=dist_cut:
            trim_out=True
            #print(tot_dist,get_vec(id1,mlist),get_vec(id2,mlist))
            #print (str(id0)+' '+str(id1)+' '+str(id2)+' '+str(id3))
            #f.write(get_vec_dir(id0,mlist)+get_vec_dir(id1,mlist)+get_vec_dir(id2,mlist)+get_vec_dir(id3,mlist)+' '+str(dist01)+' '+str(dist02)+' '+str(dist03)+' '+str(tot_dist)+'\n')
   
   with open(file_out,'w+') as f:
      for trim in list_tot:
         sid0=get_vec_dir(int(trim[0]),mlist)
         sid1=get_vec_dir(int(trim[1]),mlist)
         sid2=get_vec_dir(int(trim[2]),mlist)
         f.write(sid0+sid1+sid2+' '+' '.join(map(str,trim[3:7]))+'\n')

   #sort the list
   list_tot.sort(key=lambda x: x[6])
   #print(list_tot) 
   #write the list to file
   with open('sort_'+file_out,'w+') as f:
      for trim in list_tot:
         sid0=get_vec_dir(int(trim[0]),mlist)
         sid1=get_vec_dir(int(trim[1]),mlist)
         sid2=get_vec_dir(int(trim[2]),mlist)
         f.write(sid0+sid1+sid2+' '+' '.join(map(str,trim[3:7]))+'\n')

###########################################################################
#   gen_tetramers_dist_cache(shells,dist_cut,mlist,coord,cell,file_out)   #
###########################################################################
#
# Generate the list of tetramers using a distance criterion
# an ordered list of monomers needs to be put to the routine
# a nested loop over the monomers is performed 0<id1<id2<id3
# tetramers are found, stored in a list, and finally ordered 
# according to distance and written to the file
#
# compared to gen_tetramers() the dimer distances are precalculated
# here
############################################################
#
# #nmono_max: number of monomers to read from the monomer list mlist
# shells:    number of shells to consider
# dist_cut:  distance criterion in Angstrom
# mlist:     monomer list as from gen_list()
# coord:     monomer coordinates in the unit cell
# cell:      unit cell
# file_out:  output file to write the list to
############################################################
def gen_tetramers_dist_cache(shells,dist_cut,mlist,coord,cell,file_out):
   #number of monomers in the unit cell
   cell_mono=coord.shape[0]
   #number of shells that we are using
   #the largest index of the distance array will have to be then 2*shells+1
   #shells=mlist[nmono_max,1] 
   nmono_max=-1
   if mlist[-1,1]==shells:
      #if the last monomer is from the largest shell, then use all monomers
      nmono_max=mlist.shape[0]-1
   elif mlist[-1,1]<shells:
      print ('Error, monomer list not sufficient, rerun mbe.gen_list with a larger number of shells')
      exit()
   else:
      #find the largest monomer index from the shell
      #loop over all the monomers 
      for idx in range(mlist.shape[0]):
         #we find a monomer from a larger shell
         if mlist[idx,1]>shells:
            nmono_max=idx-1
            break
      if nmono_max==-1:
         print ('Can not identify monomer from your shell')

   #dist_ar is the array which stores all the distances
   dist_ar=np.zeros((cell_mono,cell_mono,2*shells+1,4*shells+1,4*shells+1)) 
   #finished here
   shift_base=[0,0,0]
   dist_min=100000
   #index the first monomer
   for m1 in range(cell_mono):
      #first monomer coordinates in the unit cell
      cm1=coord[m1,:,:].copy()
      for m2 in range(cell_mono):
         #second monomer coordinates in the unit cell
         cm2=coord[m2,:,:].copy()
         #difference between the x shift of the monomers, 
         #the largest difference between X shifts for a given number of shells
         #is 2*shells, from -shells to shells
         #the smallest one can be set to zero as we can choose the monomer
         for xsh in range(2*shells+1):
            #for y and z we can't shuffle the monomers, so the difference between cells
            #is from -2*shells to 2*shells, in total 4*shells
            #in the range ysh=2*shells corresponds to the zero difference between the cells
            for ysh in range(4*shells+1):
               for zsh in range(4*shells+1):
                  shift_m2=[xsh,ysh-2*shells,zsh-2*shells]
                  #dist_ar[m1,m2,xsh,ysh,zsh]=get_distance_from_vec(cm1,shift_base,cm2,shift_m2,cell)
                  dist_act=get_distance_from_vec(cm1,shift_base,cm2,shift_m2,cell)
                  dist_ar[m1,m2,xsh,ysh,zsh]=dist_act
                  #keep track of the shortest distance between dimers
                  if dist_act<dist_min:
                     dist_min=dist_act

   id0=0
   list_tot=[]
   #with open(file_out,'w+') as f:
   for id1 in range(1,nmono_max):
      print('working on '+str(id1)+' out of '+str(nmono_max))
      #dist01=get_distance_from_idx(id0,id1,mlist,coord,cell)
      #get xyz and monomer information
      dat1=mlist[id1,2:6]
      #if the xshift is 0 or larger, use the shift directly
      #need to shift the y and z index to larger values
      if dat1[0]>=0:
         dist01=dist_ar[0,dat1[3],dat1[0],dat1[1]+2*shells,dat1[2]+2*shells]
      else: 
         #if the xshift is negative, swap the monomers and shift
         dist01=dist_ar[dat1[3],0,-dat1[0],-dat1[1]+2*shells,-dat1[2]+2*shells]
      if dist01>dist_cut-5*dist_min:
         continue
      for id2 in range(id1+1,nmono_max):
         #dist02=get_distance_from_idx(id0,id2,mlist,coord,cell)
         dat2=mlist[id2,2:6]
         #if the xshift is 0 or larger, use the shift directly
         if dat2[0]>=0:
            dist02=dist_ar[0,dat2[3],dat2[0],dat2[1]+2*shells,dat2[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift
            dist02=dist_ar[dat2[3],0,-dat2[0],-dat2[1]+2*shells,-dat2[2]+2*shells]
         #dist12=get_distance_from_idx(id1,id2,mlist,coord,cell)
         #if the xshift is same or larger than that of monomer1, use the shift directly
         if dat2[0]>=dat1[0]:
            dist12=dist_ar[dat1[3],dat2[3],dat2[0]-dat1[0],dat2[1]-dat1[1]+2*shells,dat2[2]-dat1[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift
            dist12=dist_ar[dat2[3],dat1[3],dat1[0]-dat2[0],dat1[1]-dat2[1]+2*shells,dat1[2]-dat2[2]+2*shells]
         if dist01+dist02+dist12>dist_cut-3*dist_min:
            continue
         for id3 in range(id2+1,nmono_max):
            dat3=mlist[id3,2:6]
            #dist03=get_distance_from_idx(id0,id3,mlist,coord,cell)
            #if the xshift is 0 or larger, use the shift directly
            if dat3[0]>=0:
               dist03=dist_ar[0,dat3[3],dat3[0],dat3[1]+2*shells,dat3[2]+2*shells]
            else: 
               #if the xshift is negative, swap the monomers and shift
               dist03=dist_ar[dat3[3],0,-dat3[0],-dat3[1]+2*shells,-dat3[2]+2*shells]
            #dist13=get_distance_from_idx(id1,id3,mlist,coord,cell)
            if dat3[0]>=dat1[0]:
               dist13=dist_ar[dat1[3],dat3[3],dat3[0]-dat1[0],dat3[1]-dat1[1]+2*shells,dat3[2]-dat1[2]+2*shells]
            else: 
               #if the xshift is negative, swap the monomers and shift
               dist13=dist_ar[dat3[3],dat1[3],dat1[0]-dat3[0],dat1[1]-dat3[1]+2*shells,dat1[2]-dat3[2]+2*shells]

            #dist23=get_distance_from_idx(id2,id3,mlist,coord,cell)
            if dat3[0]>=dat2[0]:
               dist23=dist_ar[dat2[3],dat3[3],dat3[0]-dat2[0],dat3[1]-dat2[1]+2*shells,dat3[2]-dat2[2]+2*shells]
            else: 
               #if the xshift is negative, swap the monomers and shift
               dist23=dist_ar[dat3[3],dat2[3],dat2[0]-dat3[0],dat2[1]-dat3[1]+2*shells,dat2[2]-dat3[2]+2*shells]

            tot_dist=dist01+dist02+dist12+dist03+dist13+dist23
            if tot_dist<=dist_cut:
               #print(tot_dist,get_vec(id1,mlist),get_vec(id2,mlist))
               #print (str(id0)+' '+str(id1)+' '+str(id2)+' '+str(id3))
               #f.write(get_vec_dir(id0,mlist)+get_vec_dir(id1,mlist)+get_vec_dir(id2,mlist)+get_vec_dir(id3,mlist)+' '+str(dist01)+' '+str(dist02)+' '+str(dist03)+' '+str(tot_dist)+'\n')
               list_tot.append([id0,id1,id2,id3,dist01,dist02,dist03,tot_dist])
   list_tot.sort(key=lambda x: x[7])
   #print(list_tot) 
   with open(file_out,'w+') as f:
      for tetr in list_tot:
         sid0=get_vec_dir(int(tetr[0]),mlist)
         sid1=get_vec_dir(int(tetr[1]),mlist)
         sid2=get_vec_dir(int(tetr[2]),mlist)
         sid3=get_vec_dir(int(tetr[3]),mlist)
         f.write(sid0+sid1+sid2+sid3+' '+' '.join(map(str,tetr[4:8]))+'\n')

################################################################################
#   gen_tetramers_dist_cache_cell(shells,dist_cut,mlist,coord,cell,file_out)   #
################################################################################
#
# Generate the list of tetramers using a distance criterion
# an ordered list of monomers needs to be put to the routine
# a nested loop over the monomers is performed 0<id1<id2<id3
# tetramers are found, stored in a list, and finally ordered 
# according to distance and written to the file
#
# compared to gen_tetramers() the dimer distances are precalculated
# here
############################################################
#
# #nmono_max: number of monomers to read from the monomer list mlist
# shells:    number of shells to consider
# dist_cut:  distance criterion in Angstrom
# mlist:     monomer list as from gen_list()
# coord:     monomer coordinates in the unit cell
# cell:      unit cell
# file_out:  output file to write the list to
############################################################
def gen_tetramers_dist_cache_cell(shells,dist_cut,mlist,coord,cell,file_out):
   #number of monomers in the unit cell
   cell_mono=coord.shape[0]
   #number of shells that we are using
   #the largest index of the distance array will have to be then 2*shells+1
   #shells=mlist[nmono_max,1] 
   nmono_max=-1
   if mlist[-1,1]==shells:
      #if the last monomer is from the largest shell, then use all monomers
      nmono_max=mlist.shape[0]-1
   elif mlist[-1,1]<shells:
      print ('Error, monomer list not sufficient, rerun mbe.gen_list with a larger number of shells')
      exit()
   else:
      #find the largest monomer index from the shell
      #loop over all the monomers 
      for idx in range(mlist.shape[0]):
         #we find a monomer from a larger shell
         if mlist[idx,1]>shells:
            nmono_max=idx-1
            break
      if nmono_max==-1:
         print ('Can not identify monomer from your shell')

   #dist_ar is the array which stores all the distances
   dist_ar=np.zeros((cell_mono,cell_mono,2*shells+1,4*shells+1,4*shells+1)) 
   #finished here
   shift_base=[0,0,0]
   dist_min=100000
   #index the first monomer
   for m1 in range(cell_mono):
      #first monomer coordinates in the unit cell
      cm1=coord[m1,:,:].copy()
      for m2 in range(cell_mono):
         #second monomer coordinates in the unit cell
         cm2=coord[m2,:,:].copy()
         #difference between the x shift of the monomers, 
         #the largest difference between X shifts for a given number of shells
         #is 2*shells, from -shells to shells
         #the smallest one can be set to zero as we can choose the monomer
         for xsh in range(2*shells+1):
            #for y and z we can't shuffle the monomers, so the difference between cells
            #is from -2*shells to 2*shells, in total 4*shells
            #in the range ysh=2*shells corresponds to the zero difference between the cells
            for ysh in range(4*shells+1):
               for zsh in range(4*shells+1):
                  shift_m2=[xsh,ysh-2*shells,zsh-2*shells]
                  #dist_ar[m1,m2,xsh,ysh,zsh]=get_distance_from_vec(cm1,shift_base,cm2,shift_m2,cell)
                  dist_act=get_distance_from_vec(cm1,shift_base,cm2,shift_m2,cell)
                  dist_ar[m1,m2,xsh,ysh,zsh]=dist_act
                  #keep track of the shortest distance between dimers
                  if dist_act<dist_min:
                     dist_min=dist_act

   id0=0
   list_tot=[]
   #with open(file_out,'w+') as f:
   for id1 in range(1,nmono_max):
      print('working on '+str(id1)+' out of '+str(nmono_max))
      #dist01=get_distance_from_idx(id0,id1,mlist,coord,cell)
      #get xyz and monomer information
      dat1=mlist[id1,2:6]
      #if the xshift is 0 or larger, use the shift directly
      #need to shift the y and z index to larger values
      if dat1[0]>=0:
         dist01=dist_ar[0,dat1[3],dat1[0],dat1[1]+2*shells,dat1[2]+2*shells]
      else: 
         #if the xshift is negative, swap the monomers and shift
         dist01=dist_ar[dat1[3],0,-dat1[0],-dat1[1]+2*shells,-dat1[2]+2*shells]
      if dist01>dist_cut-5*dist_min:
         continue
      for id2 in range(id1+1,nmono_max):
         #dist02=get_distance_from_idx(id0,id2,mlist,coord,cell)
         dat2=mlist[id2,2:6]
         #if the xshift is 0 or larger, use the shift directly
         if dat2[0]>=0:
            dist02=dist_ar[0,dat2[3],dat2[0],dat2[1]+2*shells,dat2[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift
            dist02=dist_ar[dat2[3],0,-dat2[0],-dat2[1]+2*shells,-dat2[2]+2*shells]
         #dist12=get_distance_from_idx(id1,id2,mlist,coord,cell)
         #if the xshift is same or larger than that of monomer1, use the shift directly
         if dat2[0]>=dat1[0]:
            dist12=dist_ar[dat1[3],dat2[3],dat2[0]-dat1[0],dat2[1]-dat1[1]+2*shells,dat2[2]-dat1[2]+2*shells]
         else: 
            #if the xshift is negative, swap the monomers and shift
            dist12=dist_ar[dat2[3],dat1[3],dat1[0]-dat2[0],dat1[1]-dat2[1]+2*shells,dat1[2]-dat2[2]+2*shells]
         if dist01+dist02+dist12>dist_cut-3*dist_min:
            continue
         #make a loop over the last monomer of the tetramer
         #tetramers to be added
         list_tmp=[]
         trim_out=False
         #this is the array storing information about previous tested monomer
         dat3_tmp=[-100,-100,-100,-100]
         for id3 in range(id2+1,nmono_max):
            #get the information about current monomer
            dat3=mlist[id3,2:6]

            #check if we have moved to a new cell
            if (dat3[0]!=dat3_tmp[0]) or (dat3[1]!=dat3_tmp[1]) or (dat3[2]!=dat3_tmp[2]):
               #check if some of the previous tetramers is below cut-off
               if trim_out==True:
                  print(len(list_tmp),' tetramers added')
                  #add the tetramers to the global list
                  for elem in list_tmp:
                     list_tot.append(elem)
               #start making a new list for this cell
               list_tmp=[]
               trim_out=False
               dat3_tmp=dat3

            #dist03=get_distance_from_idx(id0,id3,mlist,coord,cell)
            #if the xshift is 0 or larger, use the shift directly
            if dat3[0]>=0:
               dist03=dist_ar[0,dat3[3],dat3[0],dat3[1]+2*shells,dat3[2]+2*shells]
            else: 
               #if the xshift is negative, swap the monomers and shift
               dist03=dist_ar[dat3[3],0,-dat3[0],-dat3[1]+2*shells,-dat3[2]+2*shells]
            #dist13=get_distance_from_idx(id1,id3,mlist,coord,cell)
            if dat3[0]>=dat1[0]:
               dist13=dist_ar[dat1[3],dat3[3],dat3[0]-dat1[0],dat3[1]-dat1[1]+2*shells,dat3[2]-dat1[2]+2*shells]
            else: 
               #if the xshift is negative, swap the monomers and shift
               dist13=dist_ar[dat3[3],dat1[3],dat1[0]-dat3[0],dat1[1]-dat3[1]+2*shells,dat1[2]-dat3[2]+2*shells]

            #dist23=get_distance_from_idx(id2,id3,mlist,coord,cell)
            if dat3[0]>=dat2[0]:
               dist23=dist_ar[dat2[3],dat3[3],dat3[0]-dat2[0],dat3[1]-dat2[1]+2*shells,dat3[2]-dat2[2]+2*shells]
            else: 
               #if the xshift is negative, swap the monomers and shift
               dist23=dist_ar[dat3[3],dat2[3],dat2[0]-dat3[0],dat2[1]-dat3[1]+2*shells,dat2[2]-dat3[2]+2*shells]

            tot_dist=dist01+dist02+dist12+dist03+dist13+dist23
            list_tmp.append([id0,id1,id2,id3,dist01,dist02,dist03,dist12,dist13,dist23,tot_dist])
            if tot_dist<=dist_cut:
               trim_out=True
               #print(tot_dist,get_vec(id1,mlist),get_vec(id2,mlist))
               #print (str(id0)+' '+str(id1)+' '+str(id2)+' '+str(id3))
               #f.write(get_vec_dir(id0,mlist)+get_vec_dir(id1,mlist)+get_vec_dir(id2,mlist)+get_vec_dir(id3,mlist)+' '+str(dist01)+' '+str(dist02)+' '+str(dist03)+' '+str(tot_dist)+'\n')
               #list_tot.append([id0,id1,id2,id3,dist01,dist02,dist03,tot_dist])

   with open(file_out,'w+') as f:
      for tetr in list_tot:
         sid0=get_vec_dir(int(tetr[0]),mlist)
         sid1=get_vec_dir(int(tetr[1]),mlist)
         sid2=get_vec_dir(int(tetr[2]),mlist)
         sid3=get_vec_dir(int(tetr[3]),mlist)
         f.write(sid0+sid1+sid2+sid3+' '+' '.join(map(str,tetr[4:11]))+'\n')


   list_tot.sort(key=lambda x: x[10])
   #print(list_tot) 
   with open('sort_'+file_out,'w+') as f:
      for tetr in list_tot:
         sid0=get_vec_dir(int(tetr[0]),mlist)
         sid1=get_vec_dir(int(tetr[1]),mlist)
         sid2=get_vec_dir(int(tetr[2]),mlist)
         sid3=get_vec_dir(int(tetr[3]),mlist)
         f.write(sid0+sid1+sid2+sid3+' '+' '.join(map(str,tetr[4:11]))+'\n')

#################################################
#   get_chematrix(mer,mlist,coord,elems,cell)   #
#################################################
#
# calculates the eigenvalues of coulomb or distance matrix
# for a cluster defined by a list of monomers (stored as list of indices
# into monomer list mlist)
#
#################################################
#
# mer list of integers which correspond to individual monomers in the crystal
#     with position according to mlist 
# mlist list of monomers (shell, XYZ, monomerid)
# coord array storing the coordinates of monomers
# elems list with monomer indices
# cell cell coordinates
#################################################
def get_chematrix(mer,mlist,coord,elems,cell):
   size=0
   for idx in range(len(mer)):
      #the monomer number is 
      #print(mlist[idx,5])
      #but now all our monomers have identical length,
      #it'd be better to either have a list of arrays with their coords
      #or use a vector storing the number of atoms
      size+=coord[mlist[mer[idx],5],:,:].shape[0]
      #print(mlist[idx,5])
      #print(size)
      #print(get_vec_dir(mer[idx],mlist))
      #print(get_mono_struc_from_idx(mer[idx],mlist,coord,cell))
   #create a matrix which will store the distances/coulomb interactions/...
   chmat=np.zeros((size,size),dtype=np.float64)
   #the x and y indices are combined indices of molecules*atoms_in_molecules + atom index
   #this is the molecules*atoms_in_molecules index
   basex=0
   #imerx is the molecules index
   for imerx in range(len(mer)):
     #get coordinates of m1
     m1=get_mono_struc_from_idx(mer[imerx],mlist,coord,cell)
     e1=elems[mlist[mer[imerx],5]]
     #print (imerx, e1)
     basey=0
     for imery in range(len(mer)):
       #get coordinates of m2
       m2=get_mono_struc_from_idx(mer[imery],mlist,coord,cell)
       e2=elems[mlist[mer[imery],5]]
       #loop over atoms of molecule 1
       for iat1 in range(m1.shape[0]):
         z1=get_Z(e1[iat1])
         #loop over atoms of molecule 2
         for iat2 in range(m2.shape[0]):
           z2=get_Z(e2[iat2])
           #print(m2[idy,:])
           #the coordinates are stored  in m1[idx,:] and m2[idy,:]
           #so we only need to use the formula to fill the chmat
           #print (iat1,iat2,np.sqrt((m1[iat1,0]-m2[iat2,0])**2+(m1[iat1,1]-m2[iat2,1])**2+(m1[iat1,2]-m2[iat2,2])**2))
           #using coulomb or distance matrix does not seem to make a difference, 
           if basex+iat1==basey+iat2:
             #chmat[basex+iat1,basey+iat2]=0. #5*z1**2.4
             chmat[basex+iat1,basey+iat2]=0.5*z1**2.4
           else:
             chmat[basex+iat1,basey+iat2]=np.sqrt((m1[iat1,0]-m2[iat2,0])**2+(m1[iat1,1]-m2[iat2,1])**2+(m1[iat1,2]-m2[iat2,2])**2) 
             #chmat[basex+iat1,basey+iat2]=z1*z2/np.sqrt((m1[iat1,0]-m2[iat2,0])**2+(m1[iat1,1]-m2[iat2,1])**2+(m1[iat1,2]-m2[iat2,2])**2)
             #chmat[basex+iat1,basey+iat2]=1./np.sqrt((m1[iat1,0]-m2[iat2,0])**2+(m1[iat1,1]-m2[iat2,1])**2+(m1[iat1,2]-m2[iat2,2])**2)
             #print(np.sqrt((m1[iat1,0]-m2[iat2,0])**2+(m1[iat1,1]-m2[iat2,1])**2+(m1[iat1,2]-m2[iat2,2])**2))
       basey+=m2.shape[0]
     basex+=m1.shape[0]
   #exit()
   #print(chmat)
   w,v = np.linalg.eig(chmat)
   w=w.real
   w.sort()
   #print(w)
   return w

#####################################################################
#   sift_nmer_list(listfile,mlist,coord,elems,cell,epsilon=1.e-6)   #
#####################################################################
#
# find symmetry equivalent nmers in a distance ordered list of nmers 
# the nmers are read from file created by gen_mers
#
# To do box summation in the future we need to produce the original list
# with symmetry equivalent mers groupped together. As this was a bug.
# Imagine we have several mers with the same distance, say there are 
# two equivalent groups, but they are not ordered, the third can belong
# to the first or second group, same for each of the next one.
# We print each inequivalent mer with the number of symmetry equivalent
# mers and when expanding the symmetry we take the mer directories 
# from the original list of mers (without symmetry applied). 
# This means that we could assign wrong mers to the groups and, moreover,
# we could use wrong name for some mer.
# Example:
# mer1 inequiv, group1
# mer2 inequiv, group2
# mer3 equiv, group1
# sift_nmer_list() =>
# mer1 multip 2  energy1
# mer2 multip 1  energy2
# expand sym =>
# mer1 energy1
# mer2 energy1
# mer3 energy2
# so that the energy assigned to mer3 is not energy1 but energy2
# also I see that the mer directory actually appears twice in these cases, 
# not clear why yet but also wrong
#
#######################################################
#
# listfile file listing the nmers
# mlist monomer list
# coord array with monomer coordinates
# elems names of elements of monomers
# cell unit cell vectors 
#######################################################
def sift_nmer_list(listfile,mlist,coord,elems,cell,epsilon=1.e-6):
   #open the list file and read all the mers in it
   with open(listfile,'r') as f:
      flin=f.readlines()
  
   #the list is ordered according to distance,
   #so for all nmers with the same distance we will 
   #build a list of inequivalent nmers. First, put the first into inequiv list
   #compare the next ones to those in the list, if equivalent, then say so,
   #if inequiv, add them to list
   #clear everything if the distance changes by more than epsilon
   #epsilon=1.e-6 
   #print (epsilon)
   #list of inequivalent nmers for each distance
   inequiv=[]
   #list of eigenfunctions of distance matrix 
   eig_list=[]
   #number of equivalent
   neq_list=[]
   #keep the complete line for the nonequivalent so that they are easy to print
   line_list=[]
   #complete line of equivalent to be able to print the nosym list of mers 
   #in correct order (with groups of equivalent together)
   all_list=[]
   #we start with meaningless distance
   dist_cur=-1
   
   shuf=open('shuffle_'+listfile,'w+')
   with open('sym_'+listfile,'w+') as outf:
      for line in flin:  
         fspl=line.split()
         #get directories
         dirs=fspl[0].strip('/').split('/')
         idxs=[]
         #make a list with monomer indices for the current nmer
         #TODO obviously replace by for dir in dirs
         for i in range(len(dirs)):
            idxs.append(get_idx_str(mlist,dirs[i]))
         if len(dirs) == 2:
            #for dimer the distance is on second place
            dist=float(fspl[1])
         elif len(dirs) == 3:
            #for trimer distance is printed on the fifth place
            #get distance of the nmer
            dist=float(fspl[4])
         elif len(dirs) == 4:
            #for tetramer distance is printed on the eigth place
            #get distance of the nmer
            dist=float(fspl[7])
         else:
            print('mbelib confused in sift_nmer_list ')
         
         #check if the distance has changed from previous significantly
         if abs(dist-dist_cur) > epsilon:
            #print('new distance'+str(dist))
            #print('Register new inequivalent mer ', str(idxs))
            #we moved to a new set of mers
            #clear the inequiv list of mers and their chematrix eigenvalues
            #store them in a file
            for inmer in range(len(inequiv)):
               print ('Found '+str(neq_list[inmer])+' equivalent mers for '+str(inequiv[inmer]))
               #write the inequivalent nmer to a file with the multiplicity included
               outf.write(line_list[inmer]+' '+str(neq_list[inmer])+'\n')
               #write the equivalent nmers to new sort_* file
               #neq_list[inmer] gives the number of equivalent nmers for inequivalent mer inmer
               for ieq in range(neq_list[inmer]):
                  shuf.write(all_list[inmer][ieq]+'\n')
            inequiv=[]
            eig_list=[]
            neq_list=[]
            line_list=[]
            all_list=[]
            #add the current mer as the first element
            inequiv.append(idxs)
            #we have now one equivalent mer for the new mer
            neq_list.append(1)
            #get the chematrix eigenvalues and add them to list 
            eigvals=get_chematrix(idxs,mlist,coord,elems,cell)
            eig_list.append(eigvals)
            line_list.append(line.strip())
            #we want to keep collecting all the lines to be able to print them in the end
            #we create a list for first monomer and use it as a first element of the list
            #if the subsequent line will be equivalent it will be appended to the first 
            #list element, if it will be inequivalent a new list in the main list will be
            #created
            #that might be valid?
            tmp_list=[line.strip()]
            all_list.append(tmp_list)
            #print('New mer ', str(idxs))
            #print('mer')
            #print(inequiv)
            #print(eig_list)
            dist_cur=dist
         else:
            #the distance has not changed
            #we have possibly an inequivalent or equivalent mer
            #get the current eigenvals and compare them to those in the inequiv list
            eigvals=get_chematrix(idxs,mlist,coord,elems,cell)
            #go through the list of inequivalent mers 
            #assume that the mer is inequivalent
            new=1
            for inmer in range(len(inequiv)):
               #print(inmer) 
               #print(eigvals)
               #print(eig_list[inmer])
               #print(eigvals-eig_list[inmer])
               tot=np.sum(np.absolute(eigvals-eig_list[inmer]))
               #print(tot)
               if tot < epsilon:
                  #we have a mer equivalent to inmer
                  #print(str(idxs)+' found equivalent to '+str(inequiv[inmer]))
                  new=0
                  equiv=inmer
                  neq_list[inmer]+=1
                  #add the current line to the appropriate all_list element
                  all_list[inmer].append(line.strip())
            #we have really a new inequivalent mer
            #store its eigenvalues and indices
            if new==1:
               #print('Found new inequivalent mer ', str(idxs))
               inequiv.append(idxs)
               eig_list.append(eigvals)
               neq_list.append(1)
               line_list.append(line.strip())
               #add the line to the all_list as a new list element
               tmp_list=[line.strip()]
               all_list.append(tmp_list)
      #we are at the end of the file, we need to print the stuff from the memory
      for inmer in range(len(inequiv)):
         #print ('Found '+str(neq_list[inmer])+' equivalent mers for '+str(inequiv[inmer]))
         #write the inequivalent nmer to a file with the multiplicity included
         outf.write(line_list[inmer]+' '+str(neq_list[inmer])+'\n')
         #write the equivalent nmers to new sort_* file
         for ieq in range(neq_list[inmer]):
            shuf.write(all_list[inmer][ieq]+'\n')
         #exit()
      
 
   
####################################################################################
#   expand_sym_nmer_list(list_file, enlist_file, out_enlist_file)   #
####################################################################################
#
# expand file with energies obtained for symmetrized lists
# needs to read the original list, symmetrized list, and file with energies 
#
# it goes through the energy file and for each line reads the symmetry factor
# and uses as many lines from the non-symmetry list to reprint the energies
#
#######################################################
#
# list_file           file listing the nmers
# sym_list_file       file listing the sym nmers (not needed)
# enlist_file         file with energies
# out_enlist_file     file to put energies to original list
#######################################################
def expand_sym_nmer_list(list_file, enlist_file, out_enlist_file):
   # open the list file and read all the mers in it
   with open(list_file,'r') as f:
      flin=f.readlines()
   # open the list file with energies of SYMs including the symmetry factors and read its contents
   with open(enlist_file,'r') as f:
      sflin=f.readlines()
   #TODO check, taken from other mbe.py file
   print(sflin[-1])

   # guess the number of mers
   test=sflin[0].split('/')
   nmer=0
   for part in test:
      #we need to check for finite length of the string as two // create an empty string 
      #  and part[0] is not valid
      if len(part) >0:
         if part[0]=='s':
            nmer+=1
         else:
            break
   #print('We have ',nmer,'-mer')
    
   # position of symmetry factor in the energy file
   if nmer==2:
     facpos=3
   elif nmer==3:
     facpos=6
   elif nmer==4:
     facpos=9
   else:
     print('need to fill in the position of the symmetry factor in expand_sym_nmer_list ')
     exit()

   # index into the non-sym list
   idx=0

   fout=open(out_enlist_file,'w') 

   #loop over the SYMs list
   for line in sflin:
      #get the symmetry factor
      symac=line.split()[facpos]
      print(line.split()[0]+' '+symac)
      #no symmetry companion for this -mer, be quick
      if symac=='1':
         #copy the line to output
         fout.write(line)
         #increment index into the non-Sym file
         idx+=1
      else:
         for i in range(int(symac)):
            #need to build the line
            lineadd=''
            #get the IDs of the mers 
            mers=flin[idx].split()[0]
            #starting bit with mers
            lineadd+=mers
            #basis set info
            #remove empty list elements and take the one after nmers
            basis=list(filter(None,line.split('/')))[nmer]
            lineadd+='/'+basis+'/ '
            # rest of the line with energies
            rest= line.split()[1:]
            #replace the symmetry factor
            rest[facpos-1]='1'
            #make a string
            lineadd+=' '.join(rest)
            #add to file
            fout.write(lineadd+'\n')
            idx+=1
 

###########################################################
#   nmer_box_reorder(enin_file, enout_file, dist_pos)   #
###########################################################
#
# read a list of mers with energies and reorder it so that
# all last monomers from one cell are together
#
# goes through the list, for each mer takes the string defining the mer
# and searches the rest of the list for the same string
#
#######################################################
#
# enin_file  file with unSYM mers and their energies
# enout_file output file with unSYM mers and their energies
# dist_pos   position of distance in the enin_file
#######################################################
def nmer_box_reorder(enin_file, enout_file, dist_pos):
   # open the file with energies
   with open(enin_file,'r') as f:
      flin=f.readlines()

   #this is the reordered list
   new_list=[]
   f=open(enout_file,'w')
   while len(flin)!= 0:
      print(len(flin))
      #We need to find the length of the substring to search for 
      #This is one way to do this
      part=flin[0].split('/')[0:3]
      len1=len(part[0])+len(part[1])+len(part[2])-4
      #print(part)
      #print(flin[0][0:len1])
      query=flin[0][0:len1]
      #list of indices to the list of energies
      idx_list=[]
      #list of rows to remove/add to the new list
      to_add_list=[]
      #search for all the occurences of the cell
      #do not delete stuff now to avoid some errors due to shifting
      #(might not be necessary, sure)
      for idx in range(len(flin)):
         if query in flin[idx]:
            idx_list.append(idx)
            to_add_list.append(flin[idx])
 
      new_dist=to_add_list[0].split()[dist_pos]
      #add all the mers from one cell to the new list
      for elem in to_add_list:
         #fake the distance 
         parts=elem.split()
         parts[dist_pos]=new_dist
         elem=' '.join(map(str,parts))
         new_list.append(elem)
         f.write(elem+'\n')
      
      for ind in reversed(idx_list):
         del flin[ind]
      
    
#################################################################
#   nmer_box_reorder_general(enin_file, enout_file, dist_pos)   #
#################################################################
#
# read a list of mers with energies and reorder it so that
# all last monomers from one cell are together
#
# goes through the list, for each mer takes the string defining the mer
# and searches the rest of the list for the same string
#
# compared to previous routine this one guesses the number of mers 
#
#######################################################
#
# enin_file  file with unSYM mers and their energies
# enout_file output file with unSYM mers and their energies
# dist_pos   position of distance in the enin_file
#######################################################
def nmer_box_reorder_general(enin_file, enout_file, dist_pos):
   # open the file with energies
   with open(enin_file,'r') as f:
      flin=f.readlines()

   #use the first line to find out how many mers we have
   part=flin[0].split('/')
   nmers=0
   #loop over the parts and check if they contain 'shell'
   #this should identify each molecule in a nmer
   for test in part:
      if "shell" in test:
         nmers+=1

   new_list=[]
   f=open(enout_file,'w')
   while len(flin)!= 0:
      print(len(flin))
      #We need to find the length of the substring to search for 
      #This is one way to do this
      part=flin[0].split('/')[0:nmers]
      #find the length of line to cut
      #this identifies the shell
      len1=0
      for i in range(nmers):
         len1+=len(part[i])
      #for dimer we need to subtract 5, not sure
      #TODO: probably needs to be changed when tetramer list changed
      len1-=5
      query=flin[0][0:len1]
      #list of indices to the list of energies
      idx_list=[]
      #list of rows to remove/add to the new list
      to_add_list=[]
      #search for all the occurences of the cell
      #do not delete stuff now to avoid some errors due to shifting
      #(might not be necessary, sure)
      for idx in range(len(flin)):
         if query in flin[idx]:
            idx_list.append(idx)
            to_add_list.append(flin[idx])
 
      new_dist=to_add_list[0].split()[dist_pos]
      #add all the mers from one cell to the new list
      for elem in to_add_list:
         #fake the distance 
         parts=elem.split()
         parts[dist_pos]=new_dist
         elem=' '.join(map(str,parts))
         new_list.append(elem)
         f.write(elem+'\n')
      
      for ind in reversed(idx_list):
         del flin[ind]


#####################################################################
#   filter_NN_nmer_list(listfile,dist_lim)   #
#####################################################################
#
# read the list of nmers and divide it into lists according to number
# of nearest neighbours (NN)
# NN is defined by distance which is compared to distances written 
# in the list file
#
#######################################################
#
# listfile file listing the nmers
# dist     distance to define nearest neighbours
#######################################################
def filter_NN_nmer_list(listfile,dist_lim):
   #open the list file and read all the mers in it
   with open(listfile,'r') as f:
      flin=f.readlines()
  
   #get the number of distances from the list file
   #the list file contains
   #fragment distances tot_dist symmetry_fac_(optional)
   # guess the fragment size "nfrag"
   # and set the number of groups "ngroup" based on contact-distance
   # the guess is made based on the index of column that contains the data
   nparts=len(flin[0].split(' '))
   if nparts==3 or nparts==2:
     nfrag=2
     ngroup=2
   elif nparts==5 or nparts==6:
     nfrag=3
     ngroup=4
   elif nparts==8 or nparts==9:
     nfrag=4
     ngroup=7
   else: 
     print ('filter_NN_nmer_list: can not figure out number of mers')
     print(flin[0])
     print(nparts)
   #check that distance is correct for the first line
   #tot_dist=flin[0].split(' ')[nfrag]
   #TODO get the distances to floats, sum them and compare to tot_dist

   #loop over listfile's lines stored in flin and write them to 
   #the correct file

   #open the files for each NN group
   fs=[]
   for i in range(ngroup):
     filename='NC'+str(i)+'_'+str(dist_lim)+'_'+listfile
     outf=open(filename,'w+')
     fs.append(outf)

   for line in flin:
     #get the distances in fragment
     dists=line.split(' ')[1:ngroup]
     #check how many distances are below cut-off
     nNN=0
     for dist in dists:
       if float(dist)<dist_lim:
         nNN+=1
     if nNN>ngroup:
       print('error in filter_NN_nmer_list, too many NNs',dists)
     #write the line to correct file
     fs[nNN].write(line)
    
   for i in range(ngroup):
     fs[i].close()

##################################
##################################
#                                #
#   getting many-body energies   #


# BIG TODO: to make these n-independent

############################
#   get_2body(data)   #
############################
# 
# function to return the 2-body interaction
# the energies on input are ordered as dimer,
# monomers
# with the order 
# ab a b 
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_2body(data):
   #for this dimer, the mers are ordered as
   # ab  a   b 
   # 0    1  2 
   E2b12=float(data[0])-float(data[1])-float(data[2])
   #print "2body ", E2b12,  E2b13, E2b23

   return E2b12

############################
#   get_3body(data)   #
############################
# 
# function to return the 3-body interactions
# the energies on input are ordered as trimer, 
# dimer, monomer
# with the exact order as
# abc bc ac ab a b c 
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_3body(data):
   #for this trimer, the mers are ordered as
   # abc  bc ac  ab  a   b c 
   # 0    1  2   3   4   5 6 
   E2b12=float(data[3])-float(data[4])-float(data[5])
   E2b13=float(data[2])-float(data[4])-float(data[6])
   E2b23=float(data[1])-float(data[5])-float(data[6])
   #print "2body ", E2b12,  E2b13, E2b23

   E3b=float(data[0])-float(data[4])-float(data[5])-float(data[6])
   E3bna=float(data[0])-float(data[4])-float(data[5])-float(data[6])-E2b12-E2b13-E2b23
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   return E3bna, E2b12, E2b13, E2b23, E3b

############################
#   get_3body_312(data)   #
############################
# 
# function to return the 3-body interactions
# the energies on input are ordered as trimer, 
# dimer, monomer
# with the exact order as
# abc a b c bc ac ab  
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_3body_312(data):
   #for this trimer, the mers are ordered as
   # abc  a  b   c   bc ac  ab
   # 0    1  2   3   4   5  6 
   E2b12=float(data[6])-float(data[1])-float(data[2])
   E2b13=float(data[5])-float(data[1])-float(data[3])
   E2b23=float(data[4])-float(data[2])-float(data[3])
   #print "2body ", E2b12,  E2b13, E2b23

   E3b=float(data[0])-float(data[1])-float(data[2])-float(data[3])
   E3bna=float(data[0])-float(data[1])-float(data[2])-float(data[3])-E2b12-E2b13-E2b23
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   return E3bna, E2b12, E2b13, E2b23, E3b

############################
#   get_3body_123(data)   #
############################
# 
# function to return the 3-body interactions
# the energies on input are ordered as monomer 
# dimer, trimer
# with the exact order as
# a b c bc ac ab abc
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_3body_123(data):
   #for this trimer, the mers are ordered as
   # a  b   c   bc ac  ab abc
   # 0   1  2   3   4   5  6 
   E2b12=float(data[5])-float(data[0])-float(data[1])
   E2b13=float(data[4])-float(data[0])-float(data[2])
   E2b23=float(data[3])-float(data[1])-float(data[2])
   #print "2body ", E2b12,  E2b13, E2b23

   E3b=float(data[6])-float(data[0])-float(data[1])-float(data[2])
   E3bna=float(data[6])-float(data[0])-float(data[1])-float(data[2])-E2b12-E2b13-E2b23
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   return E3bna, E2b12, E2b13, E2b23, E3b

############################
#   get_4body_4312(data)   #
############################
# 
# function to return the 4-body interactions
# the energies on input are ordered as tetramer, 
# trimer, monomer, dimer 
# with the exact order as
# abcd bcd acd abd abc a b c d ab ac ad bc bd cd
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_4body_4312(data):
   #for this tetramer, the mers are ordered as
   #abcd bcd acd abd abc a b c d ab ac ad bc bd cd
   # 0    1  2   3   4   5 6 7 8 9  10 11 12 13 14
   E2b12=float(data[9])-float(data[5])-float(data[6])
   E2b13=float(data[10])-float(data[5])-float(data[7])
   E2b14=float(data[11])-float(data[5])-float(data[8])
   E2b23=float(data[12])-float(data[6])-float(data[7])
   E2b24=float(data[13])-float(data[6])-float(data[8])
   E2b34=float(data[14])-float(data[7])-float(data[8])
   #print "2body ", E2b12,  E2b13, E2b14, E2b23, E2b24, E2b34

   E3b123=float(data[4])-float(data[5])-float(data[6])-float(data[7])
   E3b124=float(data[3])-float(data[5])-float(data[6])-float(data[8])
   E3b134=float(data[2])-float(data[5])-float(data[7])-float(data[8])
   E3b234=float(data[1])-float(data[6])-float(data[7])-float(data[8])
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   E3b123na=E3b123-E2b12-E2b13-E2b23
   E3b124na=E3b124-E2b12-E2b14-E2b24
   E3b134na=E3b134-E2b13-E2b14-E2b34
   E3b234na=E3b234-E2b23-E2b24-E2b34
   #print "3body ", E3b123na, E3b124na, E3b134na, E3b234na

   E4b=float(data[0])-float(data[5])-float(data[6])-float(data[7])-float(data[8])
   #print E4b
   E4bna=E4b-E2b12-E2b13-E2b14-E2b23-E2b24-E2b34-E3b123na-E3b124na-E3b134na-E3b234na
   #print "4body ", E4bna
   return E4bna, E2b12, E2b13, E2b14, E2b23, E2b24, E2b34, E3b123na, E3b124na, E3b134na, E3b234na, E4b

############################
#   get_4body_4123(data)   #
############################
# 
# function to return the 4-body interactions
# the energies on input are ordered as tetramer, 
# trimer, monomer, dimer 
# with the exact order as
# abcd a b c d ab ac ad bc bd cd bcd acd abd abc
# 0    1 2 3 4 5  6  7  8  9  10 11  12  13  14
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_4body_4123(data):
   #for this tetramer, the mers are ordered as
   # abcd a b c d ab ac ad bc bd cd bcd acd abd abc
   # 0    1 2 3 4 5  6  7  8  9  10 11  12  13  14
   E2b12=float(data[5])-float(data[1])-float(data[2])
   E2b13=float(data[6])-float(data[1])-float(data[3])
   E2b14=float(data[7])-float(data[1])-float(data[4])
   E2b23=float(data[8])-float(data[2])-float(data[3])
   E2b24=float(data[9])-float(data[2])-float(data[4])
   E2b34=float(data[10])-float(data[3])-float(data[4])
   #print "2body ", E2b12,  E2b13, E2b14, E2b23, E2b24, E2b34

   E3b123=float(data[14])-float(data[1])-float(data[2])-float(data[3])
   E3b124=float(data[13])-float(data[1])-float(data[2])-float(data[4])
   E3b134=float(data[12])-float(data[1])-float(data[3])-float(data[4])
   E3b234=float(data[11])-float(data[2])-float(data[3])-float(data[4])
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   E3b123na=E3b123-E2b12-E2b13-E2b23
   E3b124na=E3b124-E2b12-E2b14-E2b24
   E3b134na=E3b134-E2b13-E2b14-E2b34
   E3b234na=E3b234-E2b23-E2b24-E2b34
   #print "3body ", E3b123na, E3b124na, E3b134na, E3b234na

   E4b=float(data[0])-float(data[1])-float(data[2])-float(data[3])-float(data[4])
   #print E4b
   E4bna=E4b-E2b12-E2b13-E2b14-E2b23-E2b24-E2b34-E3b123na-E3b124na-E3b134na-E3b234na
   #print "4body ", E4bna
   return E4bna, E2b12, E2b13, E2b14, E2b23, E2b24, E2b34, E3b123na, E3b124na, E3b134na, E3b234na, E4b

############################
#   get_4body_MMM(data)   #
############################
# 
# function to return the 4-body interactions
# the energies on input are ordered as tetramer, 
# trimer, monomer, dimer 
# with the exact order as
# abcd a b c ab ac bc d  ad bd cd abc abd acd bcd
# 0    1 2 3 4  5  6  7  8  9  10 11  12  13  14
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_4body_MMM(data):
   #for this tetramer, the mers are ordered as
   E2b12=float(data[4])-float(data[1])-float(data[2])
   E2b13=float(data[5])-float(data[1])-float(data[3])
   E2b23=float(data[6])-float(data[2])-float(data[3])
   E2b14=float(data[8])-float(data[1])-float(data[7])
   E2b24=float(data[9])-float(data[2])-float(data[7])
   E2b34=float(data[10])-float(data[3])-float(data[7])
   #print "2body ", E2b12,  E2b13, E2b14, E2b23, E2b24, E2b34

   E3b123=float(data[11])-float(data[1])-float(data[2])-float(data[3])
   E3b124=float(data[12])-float(data[1])-float(data[2])-float(data[7])
   E3b134=float(data[13])-float(data[1])-float(data[3])-float(data[7])
   E3b234=float(data[14])-float(data[2])-float(data[3])-float(data[7])
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   E3b123na=E3b123-E2b12-E2b13-E2b23
   E3b124na=E3b124-E2b12-E2b14-E2b24
   E3b134na=E3b134-E2b13-E2b14-E2b34
   E3b234na=E3b234-E2b23-E2b24-E2b34
   #print "3body ", E3b123na, E3b124na, E3b134na, E3b234na

   E4b=float(data[0])-float(data[1])-float(data[2])-float(data[3])-float(data[7])
   #print E4b
   E4bna=E4b-E2b12-E2b13-E2b14-E2b23-E2b24-E2b34-E3b123na-E3b124na-E3b134na-E3b234na
   #print "4body ", E4bna
   return E4bna, E2b12, E2b13, E2b14, E2b23, E2b24, E2b34, E3b123na, E3b124na, E3b134na, E3b234na, E4b


############################
#   get_5body_54312(data)   #
############################
# 
# function to return the 5-body interactions
# the energies on input are ordered as pentamer, tetramer, 
# trimer, monomer, dimer 
# with the exact order as
# abcde bcde acde abde abce abcd abc abd abe acd ace ade bcd bce bde cde a b c d e ab ac ad ae bc bd be cd ce de
############################
# 
# data an array with the energies (stored as whatever, str or float)
# on output all the interactions are returned
############################
def get_5body_54312(data):
   #for this pentamer, the mers are ordered as

   # abcde bcde acde abde abce abcd abc abd abe acd ace ade bcd bce bde cde a  b  c  d  e  ab ac ad ae bc bd be cd ce de
   # 0     1    2    3     4    5    6   7   8   9  10  11  12  13  14  15  16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
   E2b12=float(data[21])-float(data[16])-float(data[17])
   E2b13=float(data[22])-float(data[16])-float(data[18])
   E2b14=float(data[23])-float(data[16])-float(data[19])
   E2b15=float(data[24])-float(data[16])-float(data[20])
   E2b23=float(data[25])-float(data[17])-float(data[18])
   E2b24=float(data[26])-float(data[17])-float(data[19])
   E2b25=float(data[27])-float(data[17])-float(data[20])
   E2b34=float(data[28])-float(data[18])-float(data[19])
   E2b35=float(data[29])-float(data[18])-float(data[20])
   E2b45=float(data[30])-float(data[19])-float(data[20])
   #print "2body ", E2b12,  E2b13, E2b14, E2b23, E2b24, E2b34

   E3b123=float(data[6])-float(data[16])-float(data[17])-float(data[18])
   E3b124=float(data[7])-float(data[16])-float(data[17])-float(data[19])
   E3b125=float(data[8])-float(data[16])-float(data[17])-float(data[20])
   E3b134=float(data[9])-float(data[16])-float(data[18])-float(data[19])
   E3b135=float(data[10])-float(data[16])-float(data[18])-float(data[20])
   E3b145=float(data[11])-float(data[16])-float(data[19])-float(data[20])
   E3b234=float(data[12])-float(data[17])-float(data[18])-float(data[19])
   E3b235=float(data[13])-float(data[17])-float(data[18])-float(data[20])
   E3b245=float(data[14])-float(data[17])-float(data[19])-float(data[20])
   E3b345=float(data[15])-float(data[18])-float(data[19])-float(data[20])
   #print "3body ", E3b123, E3b124, E3b134, E3b234

   E3b123na=E3b123-E2b12-E2b13-E2b23
   E3b124na=E3b124-E2b12-E2b14-E2b24
   E3b125na=E3b125-E2b12-E2b15-E2b25
   E3b134na=E3b134-E2b13-E2b14-E2b34
   E3b135na=E3b135-E2b13-E2b15-E2b35
   E3b145na=E3b145-E2b14-E2b15-E2b45
   E3b234na=E3b234-E2b23-E2b24-E2b34
   E3b235na=E3b235-E2b23-E2b25-E2b35
   E3b245na=E3b245-E2b24-E2b25-E2b45
   E3b345na=E3b345-E2b34-E2b35-E2b45
   #print "3body ", E3b123na, E3b124na, E3b134na, E3b234na

   E4b1234=float(data[5])-float(data[16])-float(data[17])-float(data[18])-float(data[19])
   E4b1235=float(data[4])-float(data[16])-float(data[17])-float(data[18])-float(data[20])
   E4b1245=float(data[3])-float(data[16])-float(data[17])-float(data[19])-float(data[20])
   E4b1345=float(data[2])-float(data[16])-float(data[18])-float(data[19])-float(data[20])
   E4b2345=float(data[1])-float(data[17])-float(data[18])-float(data[19])-float(data[20])

   #print E4b
   E4bna1234=E4b1234-E2b12-E2b13-E2b14-E2b23-E2b24-E2b34-E3b123na-E3b124na-E3b134na-E3b234na
   E4bna1235=E4b1235-E2b12-E2b13-E2b15-E2b23-E2b25-E2b35-E3b123na-E3b125na-E3b135na-E3b235na
   E4bna1245=E4b1245-E2b12-E2b14-E2b15-E2b24-E2b25-E2b45-E3b124na-E3b125na-E3b145na-E3b245na
   E4bna1345=E4b1345-E2b13-E2b14-E2b15-E2b34-E2b35-E2b45-E3b134na-E3b135na-E3b145na-E3b345na
   E4bna2345=E4b2345-E2b23-E2b24-E2b25-E2b34-E2b35-E2b45-E3b234na-E3b235na-E3b245na-E3b345na

   E5b=float(data[0])-float(data[16])-float(data[17])-float(data[18])-float(data[19])-float(data[20])
   E5bna=E5b-E4bna1234-E4bna1235-E4bna1245-E4bna1345-E4bna2345-E3b123na-E3b124na-E3b125na-E3b134na-E3b135na-E3b145na-E3b234na-E3b235na-E3b245na-E3b345na
   E5bna=E5bna-E2b12-E2b13-E2b14-E2b15-E2b23-E2b24-E2b25-E2b34-E2b35-E2b45

   #print "4body ", E4bna
   return E5bna, E2b12, E2b13, E2b14, E2b15, E2b23, E2b24, E2b25, E2b34, E2b35, E2b45, E3b123na, E3b124na, E3b125na, E3b134na, E3b135na, E3b145na, E3b234na, E3b235na, E3b245na, E3b345na, E4bna1234, E4bna1235, E4bna1245, E4bna1345, E4bna2345, E5b

def get_HF_MP2_Molpro_file(filen,n):
   with open(filen,'r') as f:
     filelist=f.readlines()
    
   return get_HF_MP2_Molpro(filelist,n)

#TODO: get that from kartoffel?
###########################################
#   read_stuff_Molpro(filelist,n)   #
###########################################
#
# function to grep the energies from the file (stored as a list or whatever)
#
###########################################
#
# filelist: a list with lines of the file
# method: TODO give the method as a parameter
# n: number of mers
###########################################
def read_stuff_Molpro(filelist,method,n):
   stuff=[]
   countMP2=-1
   countMPF12=-1
   for line in filelist:
      if "!RHF STATE  1.1 Ene" in line:
         #print line,
          tmp=line.split()
          HF.append(tmp[4])
      if "New reference energy" in line:
          #print line,
          tmp=line.split()
          CABS.append(tmp[3])
      #get MP2 energies
      if "DF-MP2-F12 correlation energies:" in line:
         countMP2=3
         countMPF12=6
      if countMP2==0:
         tmp=line.split()
         MP2.append(tmp[3])
      #get MP2-F12 energies
      if countMPF12==0:
         tmp=line.split()
         MPF12.append(tmp[3])
      countMP2-=1
      countMPF12-=1
   CABS=CABS[1::2]
   
   return HF, CABS, MP2, MPF12

def read_old_HF_Molpro(filelist):
   HF=[]
   for line in filelist:
      if "!RHF STATE 1.1 Ene" in line:
         #print line,
         tmp=line.split()
         HF.append(tmp[4])
   return HF

def read_MP2_Molpro(filelist):
   MP=[]
   for line in filelist:
      if " MP2 correlation energy" in line:
         #print line,
         tmp=line.split()
         MP.append(tmp[3])
   return MP

def read_MP3_Molpro(filelist):
   MP=[]
   for line in filelist:
      if " MP3 correlation energy" in line:
         tmp=line.split()
         MP.append(tmp[3])
   return MP

def read_MP4_Molpro(filelist):
   MP=[]
   for line in filelist:
      if " Total correlation energy" in line:
         #print line,
         tmp=line.split()
         MP.append(tmp[3])
   return MP



###########################################
#   read_HF_MP2_Molpro(filelist,n)   #
###########################################
#
# function to grep the energies from the file (stored as a list or whatever)
#
###########################################
#
# filelist: a list with lines of the file
# n: number of mers
###########################################
def read_HF_MP2_Molpro(filelist,n):
   HF=[]
   CABS=[]
   MP2=[]
   MPF12=[]
   countMP2=-1
   countMPF12=-1
   for line in filelist:
      if "!RHF STATE  1.1 Ene" in line:
         #print line,
          tmp=line.split()
          HF.append(tmp[4])
      if "New reference energy" in line:
          #print line,
          tmp=line.split()
          CABS.append(tmp[3])
      #get MP2 energies
      if "DF-MP2-F12 correlation energies:" in line:
         countMP2=3
         countMPF12=6
      if countMP2==0:
         tmp=line.split()
         MP2.append(tmp[3])
      #get MP2-F12 energies
      if countMPF12==0:
         tmp=line.split()
         MPF12.append(tmp[3])
      countMP2-=1
      countMPF12-=1
   CABS=CABS[1::2]
   
   return HF, CABS, MP2, MPF12

###########################################
#   read_RPA_MMM(filelist)   #
###########################################
#
# function to grep the energies from given file (stored as a list or whatever)
#
###########################################
#
# filelist: a list with lines of the file
###########################################
def read_RPA_MMM(filelist):
   DFT=[]
   EXX=[]
   RSE=[]
   RPA=[]
   for line in filelist:
      if "Converged energy" in line:
          tmp=line.split()
          DFT.append(tmp[2])
      if "HF contribution (EtotHF) " in line:
          tmp=line.split()
          EXX.append(tmp[3])
      if "Singles correction (EcSingles) " in line:
          tmp=line.split()
          RSE.append(tmp[3])
      if "Direct RPA correlation (EcRPA) " in line:
          tmp=line.split()
          RPA.append(tmp[4])
   EXX=EXX[1::2]
   RSE=RSE[1::2]
   
   return DFT, EXX, RSE, RPA

###########################################
#   read_HF_MP2_Molpro_F12new(filelist,n)   #
###########################################
#
# function to grep the energies from the file (stored as a list or whatever)
#
###########################################
#
# filelist: a list with lines of the file
# n: number of mers
###########################################
def read_HF_MP2_Molpro_F12new(filelist,df=False):
   HF=[]
   CABS=[]
   MP2=[]
   MPF12=[]
   countMP2=-1
   countMPF12=-1
   for line in filelist:
      if "!RHF STATE  1.1 Ene" in line:
         #print line,
          tmp=line.split()
          HF.append(tmp[4])
      if "New reference energy" in line:
          #print line,
          tmp=line.split()
          CABS.append(tmp[3])
      #get MP2 energies
      if df==True:
          if "DF-MP2-F12 correlation energies:" in line:
             countMP2=3
             countMPF12=6
      else:
          if " MP2-F12 correlation energies:" in line:
             countMP2=3
             countMPF12=6

      if countMP2==0:
         tmp=line.split()
         MP2.append(tmp[3])
      #get MP2-F12 energies
      if countMPF12==0:
         tmp=line.split()
         MPF12.append(tmp[3])
      countMP2-=1
      countMPF12-=1
   CABS=CABS[1::2]
   
   return HF, CABS, MP2, MPF12

###########################################
#   read_CC_F12_Molpro(filelist)   #
###########################################
#
# function to grep the energies from the file (stored as a list or whatever)
#
###########################################
#
# filelist: a list with lines of the file
# n: number of mers
###########################################
def read_CC_F12_Molpro(filelist, return_list=False, only_MP2=False):
   HF=[]
   CABS=[]
   MP2=[]
   MF12=[]
   CC=[]
   CCF12b=[]
   Tu=[]
   fac=[]
   Ts=[]
   F12a=[]
   F12b=[]
   count=-1
   count2=-1
   countMP2=-1
   countMF12=-1
   for line in filelist:
      if "!RHF STATE  1.1 Ene" in line or "!RHF STATE 1.1 Ene" in line:
         #print line,
         tmp=line.split()
         HF.append(tmp[4])

      if "New reference energy" in line:
         #print line,
         tmp=line.split()
         CABS.append(tmp[3])

      #get MP2 energies
      if "DF-MP2-F12 correlation energies:" in line:
         countMP2=3
         countMF12=6
      if countMP2==0:
         tmp=line.split()
         MP2.append(tmp[3])
      #get MP2-F12 energies
      if countMF12==0:
         tmp=line.split()
         MF12.append(tmp[3])
 
      if only_MP2!=True:
         #get CCSD correlation energy
         if "CCSD correlation energy" in line:
            #print line,
            tmp=line.split()
            CC.append(tmp[3])
    
         #get CCSD correlation energy
         if "CCSD-F12b correlation energy" in line or "CCSD-F12c correlation energy" in line:
            #print line,
            tmp=line.split()
            CCF12b.append(tmp[3])
    
         #get scaled triples energy
         #TODO check correctness that (scaled) is not necessary in string
         if "Triples (T) contribution" in line:
            if "scaled" in line:
               #print line,
               tmp=line.split()
               Ts.append(tmp[4])
            else:
               #print line,
               tmp=line.split()
               Ts.append(tmp[3])
   
    
         #get scaling factor for triples
         if "Scale factor for triples energy" in line:
            #print line,
            tmp=line.split()
            fac.append(tmp[5])
   
         #get CCSD(T)-F12a correlation energies
         if "F12a corrections for ansatz F12/3C(FIX) added to CCSD energy" in line:
            count=7
         if count==0:
            #print line,
            tmp=line.split()
            F12a.append(tmp[3])
   
         #get CCSD(T)-F12b correlation energies
         if "F12b corrections for ansatz F12/3C(FIX) added to CCSD(T)-F12a energy" in line:
            count2=7
         if count2==0:
            #print line,
            tmp=line.split()
            F12b.append(tmp[3])

      count-=1
      count2-=1
      countMP2-=1
      countMF12-=1
 
   #take only every other line for CABS 
   CABS=CABS[1::2]
   if only_MP2 != True:
      if (len(Ts)==2*len(fac)):
         Ts=Ts[1::2]
      for i in range(len(fac)):
          Tu.append(str(float(Ts[i])/float(fac[i])))
      if return_list==True:
         return [HF, CABS, MP2, MF12, CC, CCF12b, Ts, Tu, F12a, F12b]
      else:
         return HF, CABS, MP2, MF12, CC, CCF12b, Ts, Tu, F12a, F12b
   else:
      if return_list==True:
         return [HF, CABS, MP2, MF12]
      else:
         return HF, CABS, MP2, MF12

###########################################
#   read_RPA_CC_Molpro(filelist)   #
###########################################
#
# function to grep the energies from the file (stored as a list or whatever)
#
###########################################
#
# filelist: a list with lines of the file
# n: number of mers
###########################################
def read_RPA_CC_Molpro(filelist,return_list=False):
   HF=[]
   MP2=[] 
   RPA=[]
   count=-1
   for line in filelist:
      if "!RHF STATE  1.1 Ene" in line or "!RHF STATE 1.1 Ene" in line:
         #print line,
         tmp=line.split()
         HF.append(tmp[4])

      #get MP2 energies
      if "E_c (MP2-RCCD)" in line:
         tmp=line.split()
         MP2.append(tmp[3])

      #get RPA energies
      if "E_c (dRPA-I-RCCD)" in line:
         tmp=line.split()
         RPA.append(tmp[3])

   if return_list==True:
      return [HF, MP2, RPA]
   else:
      return HF, MP2, RPA


###########################################
#   read_several_Molpro(filelist)   #
###########################################
#
# function to grep energies of DFT calcs from file
# the calculations are peformed with all the theoretical approximations 
# for the largest cluster, then for smaller again one structure, several energies
#
# the returned structure will be a list with [[E1_AB,E1_A,E1_B],[E2_AB,E2_A,E2_B]]
# where 1,2,... is index of method
#
###########################################
#
# filelist: a list with lines of the file
# 
###########################################
def read_several_Molpro(filelist):
   
   #we'll loop over the output and search for lines with energy (in DFT/HF format)
   #and make them to list, 
   #if we find "Recomputing" then we'll start from zero and append the newly 
   #found energies to the existing ones 
   #this should lead to a list of lists

   #TODO this routine expects that everything worked just fine
   #print(len(filelist))
  
   #list for lists of energies
   Es=[]
   #let's keep track of how many structures we have in total
   nstruc=0
   idx=0
   for line in filelist:
      if line.find("Recomputing") != -1:
         #print(line)
         nstruc+=1
         #print(nstruc)
         idx=0
      if nstruc==1:
         if (line.find("STATE 1.1 Energy") !=-1) or (line.find("STATE  1.1 Energy")!=-1):
            #print('found', nstruc, idx)
            ees=[]
            #print(line)
            en=line.split()[4]
            ees.append(en)
            Es.append(ees)
      else:
         if (line.find("STATE 1.1 Energy") !=-1) or (line.find("STATE  1.1 Energy")!=-1):
            #print('found', nstruc, idx)
            en=line.split()[4]
            Es[idx].append(en)
            idx+=1
   #print(Es)      
   return Es


###########################################
#   read_LCCSDT_MRCC(filelist)   #
###########################################
#
# function to grep the energies from MRCC output files
# the files are cated together
#
###########################################
#
# filelist: a list with lines of the file
# n: number of mers
###########################################
def read_LCCSDT_MRCC(filelist,return_list=False):
   HF=[]
   MP2c=[] 
   CCc=[]
   CCTc=[]
   CCT=[]
   count=-1
   for line in filelist:
      if "Reference energy [au]:" in line:
         #print line,
         tmp=line.split()
         HF.append(tmp[3])

      #get LMP2 energies
      if "LMP2 correlation energy [au]:" in line:
         tmp=line.split()
         MP2c.append(tmp[4])

      #get LNO-CCSD correlation energy
      if "CCSD correlation energy + 0.5 MP2 corrections [au]" in line:
         tmp=line.split()
         CCc.append(tmp[8])

      #get LNO-CCSD(T) correlation energy
      if "CCSD(T) correlation energy + MP2 corrections [au]" in line:
         tmp=line.split()
         CCTc.append(tmp[7])

      #get LNO-CCSD(T) energies
      if "Total LNO-CCSD(T) energy" in line:
         tmp=line.split()
         CCT.append(tmp[7])

   if return_list==True:
      return [HF, MP2c, CCc, CCTc, CCT]
   else:
      return HF, MP2c, CCc, CCTc, CCT


###########################################
#   read_several_Psi4(filelist)   #
###########################################
#
# function to grep energies of DFT calcs from file
# the calculations are peformed with all the theoretical approximations 
# for the largest cluster, then for smaller again one structure, several energies
#
# the returned structure will be a list with [[E1_AB,E1_A,E1_B],[E2_AB,E2_A,E2_B]]
# where 1,2,... is index of method
#
###########################################
#
# filelist: a list with lines of the file
# 
###########################################
def read_several_Psi4(filelist):
   
   #we'll loop over the output and search for lines with energy (in DFT/HF format)
   #and make them to list, 
   #if we find "Recomputing" then we'll start from zero and append the newly 
   #found energies to the existing ones 
   #this should lead to a list of lists

   #TODO this routine expects that everything worked just fine
   #print(len(filelist))
  
   #list for lists of energies
   Es=[]
   #let's keep track of how many structures we have in total
   nstruc=0
   idx=0
   idx_max=0
   #track the number of electrons to find out when we start a new structure
   nel=0
   for line in filelist:
      #if we find that the structure changed, start
      if line.find("Electrons    =") != -1:
         nel_new=line.split()[2]
         #we either started a new structure or started completely
         if nel_new != nel:
            #store the number of electron for this set of calculations
            nel=nel_new
            #store the number of energies for one structure
            #when we go from dimer to monomer (or similar for larger clusters)
            if nstruc==1:
               idx_max=idx
            #in either case we increase the index of mer
            nstruc+=1
            #print(nstruc)
            #reset the monomer count
            idx=0
         #if we find a new structure (line with Electrons) 
         #and the idx of calculation gets to the largest one
         #we also set it to zero 
         if idx==idx_max:
            idx=0
      #first structure
      if nstruc==1:
         if (line.find("Total Energy =") !=-1):
            #print('found', nstruc, idx)
            ees=[]
            #print(line)
            en=line.split()[3]
            ees.append(en)
            Es.append(ees)
            idx+=1
      #subsequent structures
      else:
         if (line.find("Total Energy =") !=-1):
            #print('found', nstruc, idx)
            en=line.split()[3]
            #add the energy to the correct monomer
            Es[idx].append(en)
            idx+=1
   #print(Es)      
   return Es



##########################################
#   get_nbody_all(HF,CABS,MP2,MPF12,n)   #
##########################################
#
# function to calculate the n-body interactions
###########################################
#
# HF, CABS, MP2, MPF12: lists with the energies of mers from a single run
# n: number of mers
#
###########################################
def get_nbody_all(HF, CABS, MP2, MPF12, n):

   if n==5:
      #check that the length is 15 for 4 body interactions
      if len(HF)==31 and len(CABS)==31 and len(MP2)==31 and len(MPF12)==31:
         HFval=get_5body_54312(HF)
         CABSval=get_5body_54312(CABS)
         MP2val=get_5body_54312(MP2)
         MPF12val=get_5body_54312(MPF12)
         return True, CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   #we have tetramers
   if n==4:
      #check that the length is 15 for 4 body interactions
      if len(HF)==15 and len(CABS)==15 and len(MP2)==15 and len(MPF12)==15:
         HFval=get_4body_4312(HF)
         CABSval=get_4body_4312(CABS)
         MP2val=get_4body_4312(MP2)
         MPF12val=get_4body_4312(MPF12)
         return True, CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   #we have trimers
   elif n==3:
      #check that the length is 7 for 3 body interactions
      if len(HF)==7 and len(CABS)==7 and len(MP2)==7 and len(MPF12)==7:
         HFval=get_3body(HF)
         CABSval=get_3body(CABS)
         MP2val=get_3body(MP2)
         MPF12val=get_3body(MPF12)
         return True, CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   #we have dimers
   elif n==2:
      #check that the length is 3 for 2 body interactions
      if len(HF)==3 and len(CABS)==3 and len(MP2)==3 and len(MPF12)==3:
         HFval=get_2body(HF)
         CABSval=get_2body(CABS)
         MP2val=get_2body(MP2)
         MPF12val=get_2body(MPF12)
         return True, CABSval, MPF12val,  HFval, MP2val
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   else:
      #uff, that looks illegal, 
      print ('Number of mers not implemented')
      return False, -1, -1, -1, -1

##########################################
#   get_nbody_all_new(energies,n)   #
##########################################
#
# function to calculate the n-body interactions
###########################################
#
# energies: list of lists with the energies of mers from a single run
# n: number of mers 
# incl_raw: do we want to return also the raw data?
#
#####
# returns list of energies and true if correctly calculated
###########################################
def get_nbody_all_new(energies, n, incl_raw=False):

   if n==5:
      #check that the length is 15 for 4 body interactions
      #if len(HF)==31 and len(CABS)==31 and len(MP2)==31 and len(MPF12)==31:
      correct=True
      for Emet in energies:
         if len(Emet)!=31:
            correct=False
      if correct==True:
         Eall=[]
         for Emet in energies:
            Emetval=get_5body_54312(Emet)
            Eall.append(Emetval[0])
            #HFval=get_5body_54312(HF)
            #CABSval=get_5body_54312(CABS)
            #MP2val=get_5body_54312(MP2)
            #MPF12val=get_5body_54312(MPF12)
         return True, Eall # , CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, []
   if n==4:
      correct=True
      #check first that the number of energies is correct for all the terms
      for Emet in energies:
         if len(Emet)!=15:
            correct=False
      #if we have all the data, get the non-additive energy
      if correct==True:
         Eall=[]
         for Emet in energies:
            Emetval=get_4body_4312(Emet)
            Eall.append(Emetval[0])
         return True, Eall 
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, []
   if n==3:
      correct=True
      #check first that the number of energies is correct for all the terms
      for Emet in energies:
         if len(Emet)!=7:
            correct=False
      #if we have all the data, get the non-additive energy
      if correct==True:
         Eall=[]
         for Emet in energies:
            Emetval=get_3body(Emet)
            Eall.append(Emetval[0])
         return True, Eall 
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, []
   if n==2:
      correct=True
      #check first that the number of energies is correct for all the terms
      #print('energies')
      #print(energies)
      for Emet in energies:
         if len(Emet)!=3:
            correct=False
      #if we have all the data, get the non-additive energy
      if correct==True:
         Eall=[]
         for Emet in energies:
            Emetval=get_2body(Emet)
            Eall.append(Emetval)
         return True, Eall 
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, []

   else:
      #uff, that looks illegal, 
      print ('Number of mers not implemented')
      return False, []

###########################################
#   get_tetra_HF_MP2_Molpro(filelist,n)   #
###########################################
#
# function to grep the energies from the file (stored as a list or whatever)
# and calculate the n-body interactions
###########################################
#
# filelist: a list with lines of the file
# n: number of mers
###########################################
def get_HF_MP2_Molpro(filelist,n):
   HF=[]
   CABS=[]
   MP2=[]
   MPF12=[]
   countMP2=-1
   countMPF12=-1
   for line in filelist:
      if "!RHF STATE  1.1 Ene" in line:
         #print line,
          tmp=line.split()
          HF.append(tmp[4])
      if "New reference energy" in line:
          #print line,
          tmp=line.split()
          CABS.append(tmp[3])
      #get MP2 energies
      if "DF-MP2-F12 correlation energies:" in line:
         countMP2=3
         countMPF12=6
      if countMP2==0:
         tmp=line.split()
         MP2.append(tmp[3])
      #get MP2-F12 energies
      if countMPF12==0:
         tmp=line.split()
         MPF12.append(tmp[3])
      countMP2-=1
      countMPF12-=1
   CABS=CABS[1::2]
   #we have pentamers

   if n==5:
      #check that the length is 15 for 4 body interactions
      if len(HF)==31 and len(CABS)==31 and len(MP2)==31 and len(MPF12)==31:
         HFval=get_5body_54312(HF)
         CABSval=get_5body_54312(CABS)
         MP2val=get_5body_54312(MP2)
         MPF12val=get_5body_54312(MPF12)
         return True, CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   #we have tetramers
   if n==4:
      #check that the length is 15 for 4 body interactions
      if len(HF)==15 and len(CABS)==15 and len(MP2)==15 and len(MPF12)==15:
         HFval=get_4body_4312(HF)
         CABSval=get_4body_4312(CABS)
         MP2val=get_4body_4312(MP2)
         MPF12val=get_4body_4312(MPF12)
         return True, CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   #we have trimers
   elif n==3:
      #check that the length is 7 for 3 body interactions
      if len(HF)==7 and len(CABS)==7 and len(MP2)==7 and len(MPF12)==7:
         HFval=get_3body(HF)
         CABSval=get_3body(CABS)
         MP2val=get_3body(MP2)
         MPF12val=get_3body(MPF12)
         return True, CABSval[0], MPF12val[0],  HFval[0], MP2val[0]
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   #we have dimers
   elif n==2:
      #check that the length is 3 for 2 body interactions
      #print(HF, CABS)
      #print(MP2,MPF12)
      if len(HF)==3 and len(CABS)==3 and len(MP2)==3 and len(MPF12)==3:
         HFval=get_2body(HF)
         CABSval=get_2body(CABS)
         MP2val=get_2body(MP2)
         MPF12val=get_2body(MPF12)
         return True, CABSval, MPF12val,  HFval, MP2val
      else:
         print ('Not all required energies found, check also molpro version!')
         return False, -1, -1, -1, -1
   else:
      #uff, that looks illegal, 
      print ('Number of mers not implemented')
      return False, -1, -1, -1, -1


def basis_extrap(val0,val1,bas0,bas1):
    return (bas1**3*val1-bas0**3*val0)/(bas1**3-bas0**3)

def basis_extrap5(val0,val1,bas0,bas1):
    return (bas1**5*val1-bas0**5*val0)/(bas1**5-bas0**5)

def basis_extrap7(val0,val1,bas0,bas1):
    return (bas1**7*val1-bas0**7*val0)/(bas1**7-bas0**7)

###############################################
####    ROUTINES FOR PROCESSING RESULTS    ####
###############################################

##
#
# dirfile	directory with the file
# enfile	the file with the results
# idxdist	column index of distances
# idxdata	column index of first data
# idxsym	column index of symmetry factor
# reverse	should we do the summation beckwards?
# hassym	do the data contain symmetry factor?
#

def sum_energies_file(dirfile,enfile,idxdist,idxdata,idxsym,reverse=False,hassym=True):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()
   
   hf=0.
   mp=0.
   ca=0.
   f12=0.
   #open the output file
   if reverse==False:
      outfile=dirfile+'/part_sum_'+enfile
   else:
      outfile=dirfile+'/part_revsum_'+enfile
   with open(outfile,'w+') as out:
      #with open('box15_sum.dat','w+') as out:
      if reverse==False:
         for line in data:
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               if hassym:
                  ca+=float(lsp[idxdata])*float(lsp[idxsym])
                  f12+=float(lsp[idxdata+1])*float(lsp[idxsym])
                  hf+=float(lsp[idxdata+2])*float(lsp[idxsym])
                  mp+=float(lsp[idxdata+3])*float(lsp[idxsym])
               else:
                  ca+=float(lsp[idxdata])
                  f12+=float(lsp[idxdata+1])
                  hf+=float(lsp[idxdata+2])
                  mp+=float(lsp[idxdata+3])
                  
            out.write(lsp[idxdist]+' '+str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp)+'\n')
      else:
         for line in reversed(data):
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               if hassym:
                  ca+=float(lsp[idxdata])*float(lsp[idxsym])
                  f12+=float(lsp[idxdata+1])*float(lsp[idxsym])
                  hf+=float(lsp[idxdata+2])*float(lsp[idxsym])
                  mp+=float(lsp[idxdata+3])*float(lsp[idxsym])
               else:
                  ca+=float(lsp[idxdata])
                  f12+=float(lsp[idxdata+1])
                  hf+=float(lsp[idxdata+2])
                  mp+=float(lsp[idxdata+3])
             
            out.write(lsp[idxdist]+' '+str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp)+'\n')
   
      #t(str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp))
   return ca, f12, hf, mp

#####################
# sum_energies_file1
#####################
#
# dirfile	directory with the file
# enfile	the file with the results
# idxdist	column index of distances
# idxdata	column index of first data
# idxsym	column index of symmetry factor
# nelem		how many data columns to processes
# reverse	should we do the summation beckwards?
# hassym	do the data contain symmetry factor?
#

def sum_energies_file1(dirfile,enfile,idxdist,idxdata,idxsym,nelem,reverse=False,hassym=True):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()
   
   #the total sums will be stored here
   results=np.zeros((nelem))

   #open the output file
   if reverse==False:
      outfile=dirfile+'/part_sum_'+enfile
   else:
      outfile=dirfile+'/part_revsum_'+enfile
   with open(outfile,'w+') as out:
      #with open('box15_sum.dat','w+') as out:
      if reverse==False:
         for line in data:
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               if hassym:
                  for j in range(nelem):
                     results[j]+=float(lsp[idxdata+j])*float(lsp[idxsym])
                     #f12+=float(lsp[idxdata+1])*float(lsp[idxsym])
               else:
                  for j in range(nelem):
                     results[j]+=float(lsp[idxdata+j])
                     #ca+=float(lsp[idxdata])
                  
            #out.write(lsp[idxdist]+' '+str(results[0])+' '+str(results[1])+' '+str(results[2])+' '+ str(results[3])+'\n')
            out.write(lsp[idxdist]+' '+' '.join(map(str,results[:]))+'\n')
      else:
         for line in reversed(data):
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               if hassym:
                  for j in range(nelem):
                     results[j]+=float(lsp[idxdata+j])*float(lsp[idxsym])
               else:
                  for j in range(nelem):
                     results[j]+=float(lsp[idxdata+j])
                     #ca+=float(lsp[idxdata])
             
            #out.write(lsp[idxdist]+' '+str(results[0])+' '+str(results[1])+' '+str(results[2])+' '+ str(results[3])+'\n')
            out.write(lsp[idxdist]+' '+' '.join(map(str,results[:]))+'\n' )
   
      #t(str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp))
   return results




###############################################
####      RPA PROCESSING RESULTS           ####
###############################################
#### by PNK

def sum_energies_file_rpa(dirfile,enfile,idxdist,idxdata,idxsym,reverse=False):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()

   rpa=0.
   dft=0.
   hf=0.
   c=0.
   crpas=0.
   crpax=0.
   crpaxx=0.
   #open the output file
   if reverse==False:
      outfile=dirfile+'/part_sum_'+enfile
   else:
      outfile=dirfile+'/part_revsum_'+enfile
   with open(outfile,'w+') as out:
      #with open('box15_sum.dat','w+') as out:
      if reverse==False:
         for line in data:
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               rpa+=float(lsp[idxdata])*float(lsp[idxsym])
               dft+=float(lsp[idxdata+1])*float(lsp[idxsym])
               hf+=float(lsp[idxdata+2])*float(lsp[idxsym])
               c+=float(lsp[idxdata+3])*float(lsp[idxsym])
               crpas+=float(lsp[idxdata+4])*float(lsp[idxsym])
               crpax+=float(lsp[idxdata+5])*float(lsp[idxsym])
               crpaxx+=float(lsp[idxdata+6])*float(lsp[idxsym])
            out.write(lsp[idxdist]+' '+str(rpa)+' '+str(dft)+' '+str(hf)+' '+str(c)+' '+str(crpas)+' '+str(crpax)+' '+str(crpaxx)+'\n')
      else:
         for line in reversed(data):
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               rpa+=float(lsp[idxdata])*float(lsp[idxsym])
               dft+=float(lsp[idxdata+1])*float(lsp[idxsym])
               hf+=float(lsp[idxdata+2])*float(lsp[idxsym])
               c+=float(lsp[idxdata+3])*float(lsp[idxsym])
               crpas+=float(lsp[idxdata+4])*float(lsp[idxsym])
               crpax+=float(lsp[idxdata+5])*float(lsp[idxsym])
               crpaxx+=float(lsp[idxdata+6])*float(lsp[idxsym])
            out.write(lsp[idxdist]+' '+str(rpa)+' '+str(dft)+' '+str(hf)+' '+str(c)+' '+str(crpas)+' '+str(crpax)+' '+str(crpaxx)+'\n')

      #t(str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp))
   return rpa, dft, hf, c, crpas, crpax, crpaxx

###############################################
####      DIFFERENT BASIS-SETS RPA PROCESSING RESULTS           ####
###############################################
#### by PNK

def sum_different_basis_energies_file_rpa(dirfile,enfile,idxdist,idxdata,idxsym,reverse=False):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()

   rpa=0.
   rpaxx=0.
   rpax=0.
   dftxx=0.
   dftx=0.
   hfxx=0.
   hfx=0.
   rsexx=0.
   rsex=0.
   crpa=0.
   crpax=0.
   crpaxx=0.
   #open the output file
   if reverse==False:
      outfile=dirfile+'/part_sum_'+enfile
   else:
      outfile=dirfile+'/part_revsum_'+enfile
   with open(outfile,'w+') as out:
      #with open('box15_sum.dat','w+') as out:
      if reverse==False:
         for line in data:
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               rpa+=float(lsp[idxdata])*float(lsp[idxsym])
               rpaxx+=float(lsp[idxdata+1])*float(lsp[idxsym])
               rpax+=float(lsp[idxdata+2])*float(lsp[idxsym])
               dftxx+=float(lsp[idxdata+3])*float(lsp[idxsym])
               dftx+=float(lsp[idxdata+4])*float(lsp[idxsym])
               hfxx+=float(lsp[idxdata+5])*float(lsp[idxsym])
               hfx+=float(lsp[idxdata+6])*float(lsp[idxsym])
               rsexx+=float(lsp[idxdata+7])*float(lsp[idxsym])
               rsex+=float(lsp[idxdata+8])*float(lsp[idxsym])
               crpa+=float(lsp[idxdata+9])*float(lsp[idxsym])
               crpaxx+=float(lsp[idxdata+10])*float(lsp[idxsym])
               crpax+=float(lsp[idxdata+11])*float(lsp[idxsym])
            out.write(lsp[idxdist]+' '+str(rpa)+' '+str(rpaxx)+' '+str(rpax)+' '+str(dftxx)+' '+str(dftx)+' '+str(hfxx)+' '+str(hfx)+' '+str(rsexx)+' '+str(rsex)+' '+str(crpa)+' '+str(crpaxx)+' '+str(crpax)+'\n')
      else:
         for line in reversed(data):
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               rpa+=float(lsp[idxdata])*float(lsp[idxsym])
               rpaxx+=float(lsp[idxdata+1])*float(lsp[idxsym])
               rpax+=float(lsp[idxdata+2])*float(lsp[idxsym])
               dftxx+=float(lsp[idxdata+3])*float(lsp[idxsym])
               dftx+=float(lsp[idxdata+4])*float(lsp[idxsym])
               hfxx+=float(lsp[idxdata+5])*float(lsp[idxsym])
               hfx+=float(lsp[idxdata+6])*float(lsp[idxsym])
               rsexx+=float(lsp[idxdata+7])*float(lsp[idxsym])
               rsex+=float(lsp[idxdata+8])*float(lsp[idxsym])
               crpa+=float(lsp[idxdata+9])*float(lsp[idxsym])
               crpaxx+=float(lsp[idxdata+10])*float(lsp[idxsym])
               crpax+=float(lsp[idxdata+11])*float(lsp[idxsym])
            out.write(lsp[idxdist]+' '+str(rpa)+' '+str(rpaxx)+' '+str(rpax)+' '+str(dftxx)+' '+str(dftx)+' '+str(hfxx)+' '+str(hfx)+' '+str(rsexx)+' '+str(rsex)+' '+str(crpa)+' '+str(crpaxx)+' '+str(crpax)+'\n')

      #t(str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp))
   return rpa, rpaxx, rpax, dftxx, dftx, hfxx, hfx, rsexx, rsex, crpa, crpaxx, crpax

###############################################
####    CCSD(T) PROCESSING RESULTS         ####
###############################################
#### by PNK

def sum_energies_file_ccsd(dirfile,enfile,idxdist,idxdata,idxsym,reverse=False):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()

   hf=0.
   ccsd=0.
   ca=0.
   ccsdf12=0.
   tr=0.
   #open the output file
   if reverse==False:
      outfile=dirfile+'/part_sum_'+enfile
   else:
      outfile=dirfile+'/part_revsum_'+enfile
   with open(outfile,'w+') as out:
      #with open('box15_sum.dat','w+') as out:
      if reverse==False:
         for line in data:
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               ca+=float(lsp[idxdata])*float(lsp[idxsym])
               ccsdf12+=float(lsp[idxdata+1])*float(lsp[idxsym])
               hf+=float(lsp[idxdata+2])*float(lsp[idxsym])
               ccsd+=float(lsp[idxdata+3])*float(lsp[idxsym])
               tr+=float(lsp[idxdata+4])*float(lsp[idxsym])
            out.write(lsp[idxdist]+' '+str(ca)+' '+str(ccsdf12)+' '+str(hf)+' '+str(ccsd)+' '+str(tr)+'\n')
      else:
         for line in reversed(data):
            lsp=line.split()
            #add the calculated data
            if lsp[1]=='DONE':
               ca+=float(lsp[idxdata])*float(lsp[idxsym])
               ccsdf12+=float(lsp[idxdata+1])*float(lsp[idxsym])
               hf+=float(lsp[idxdata+2])*float(lsp[idxsym])
               ccsd+=float(lsp[idxdata+3])*float(lsp[idxsym])
               tr+=float(lsp[idxdata+4])*float(lsp[idxsym])
            out.write(lsp[idxdist]+' '+str(ca)+' '+str(ccsdf12)+' '+str(hf)+' '+str(ccsd)+' '+str(tr)+'\n')

      #t(str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp))
   return ca, ccsdf12, hf, ccsd, tr

###################################
####    sum_energies_file_NN   ####
###################################
#
# read the results and produce separate sums according to number of neighbours
# (distance to define neighbour is given and the distances need to be in the results file)
#
#########################
#
# dirfile	directory with the file
# enfile	the file with the results
# idxdist	column index of distances
# idxdata	column index of first data
# idxsym	column index of symmetry factor
# dist_lim      distance limit for distinguishing compact and spread fragments
# nelem         number of energies to process (MP2: 4, CCSD(T)-F12: 5, RPA:7)
# reverse	should we do the summation beckwards?
# hassym	do the data contain symmetry factor?
#

def sum_energies_file_NN(dirfile,enfile,idxdist,idxdata,idxsym,dist_lim,nelem,reverse=False,hassym=True):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()

   # guess the fragment size "nfrag"
   # and set the number of groups "ngroup" based on contact-distance
   # the guess is made based on the index of column that contains the data
   if idxdist==2:
      nfrag=2
      ngroup=2
   elif idxdist==5:
      nfrag=3
      ngroup=4
   elif idxdist==8:
      nfrag=4
      ngroup=7
   else:
      print('sum_energies_file_NN: can not guess the fragment size')
      print('will exit now')
      exit()

   # this array stores the energies for each group
   results=np.zeros((ngroup,nelem))
 
   #loop over the data first time and find group for each mer
   #first create an array that will store the number of close contacts
   igroup=np.zeros(len(data),dtype=np.int)

   #loop over the data and assign monomers "i" into groups acc. to number of neighbours
   for i in range(len(data)):
       #get the distances
       if nfrag==2:
          #for dimers, the dimer distance is just the total distance
          lens=data[i].split()[idxdist]
          if float(lens)<dist_lim: 
             igroup[i]=1
          else:
             igroup[i]=0
       else:
          #for larger mers, take the part in front of total distance
          #careful, the retard upper index means we take one less (which we want)
          #the distances should start at index 2
          lens=data[i].split()[idxdist-ngroup+1:idxdist]
          #loop over the distances and get the number of contacts
          for x in lens:
             if float(x)<dist_lim:
                igroup[i]+=1

   #not sure how to make an array of files so let's open one group file at one time...
   #we start by searching for mers with no contact
   #loop over number of contacts
   for i in range(ngroup):
      #name of output file
      if reverse==False:
         outfile=dirfile+'/part_sum_NC'+str(i)+'_'+enfile
         outfile2=dirfile+'/NC'+str(i)+'_'+enfile
      else:
         #if the sum is done from the largest cut-off, add "rev" to the name
         outfile=dirfile+'/part_revsum_NC'+str(i)+'_'+enfile
         outfile2=dirfile+'/NC'+str(i)+'_'+enfile
      out2=open(outfile2,'w+')
      idx=0
      #open output file
      with open(outfile,'w+') as out:
         if reverse==False:
            #loop over data
            for line in data:
               #does this line belong to the current NN group?
               if igroup[idx]==i:
                  lsp=line.split()
                  if lsp[1]=='DONE':
                     #symmetry factor present?
                     if hassym:
                        #loop over all values
                        for j in range(nelem):
                           results[i,j]+=float(lsp[idxdata+j])*float(lsp[idxsym])
                     else:
                        #loop over all values
                        for j in range(nelem):
                           results[i,j]+=float(lsp[idxdata+j])
                  #out.write(lsp[idxdist]+' '+str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp)+'\n')   
                  out.write(lsp[idxdist]+' '+' '.join(map(str,results[i,:]))+'\n')
                  out2.write(line)
               idx+=1
         else:
            #loop over data
            for line in reversed(data):
               #does this line belong to the current NN group?
               if igroup[idx]==i:
                  lsp=line.split()
                  if lsp[1]=='DONE':
                     #symmetry factor present?
                     if hassym:
                        #loop over all values
                        for j in range(nelem):
                           results[i,j]+=float(lsp[idxdata+j])*float(lsp[idxsym])
                     else:
                        #loop over all values
                        for j in range(nelem):
                           results[i,j]+=float(lsp[idxdata+j])
                  #out.write(lsp[idxdist]+' '+str(ca)+' '+str(f12)+' '+str(hf)+' '+ str(mp)+'\n')   
                  out.write(lsp[idxdist]+' '+' '.join(map(str,results[i,:]))+'\n')
                  out2.write(line)
               idx+=1
   return results

###############################
####    get_3b_ATM_fac     ####
###############################
##
# read a trimer results file and calculate the ATM factor assuming average distance 
# is a reasonable approximation and C_9 is isotropic
#
# dirfile	directory with the file
# enfile	the file with the results
# idxdist	column index of distances
# idxdata	column index of first data
# idxsym	column index of symmetry factor
# hassym	do the data contain symmetry factor?
#


#def get_3b_ATM_fac(dirfile,enfile,idxdist,idxdata,idxsym,hassym=True):
def get_3b_ATM_fac(dirfile,enfile,idxdist):

   with open(dirfile+'/'+enfile,'r') as inp:
      data=inp.readlines()

   with open(dirfile+'/atm_'+enfile,'w+') as out:
      #loop over the data
      for lin in data:
          #careful! python indices
          lens=lin.split()[idxdist-3:idxdist]
          r0=float(lens[0])
          r1=float(lens[1])
          r2=float(lens[2])
          a0=0.5*(r1**2+r2**2-r0**2)/(r1*r2)
          a1=0.5*(r2**2+r0**2-r1**2)/(r2*r0)
          a2=0.5*(r0**2+r1**2-r2**2)/(r0*r1)
          fac=(1+3*a0*a1*a2)/(r0*r1*r2)**3
          out.write(lin.strip()+' '+str(fac)+'\n')
  

######################
####   elem_order ####
######################
# Orders atoms in an Atoms structure according to element
# number 
#
# molec: ASE Atoms structure to order elementwise
# 
def elem_order(molec):
    import ase
    mol_order=molec.copy()
    del mol_order[:]
    #loop over all possible elements
    for iel in range(100):
        #loop over all the atoms
        for iat in range(len(molec)):
            if molec.get_atomic_numbers()[iat] == iel:
                mol_order.append(molec[iat])
    return mol_order

#######################
####   monomerise  ####
#######################
# Turns a solid or molecular cluster into monomers
# uses ASE to get information about cell and so on
#
# 
# a       input ASE Atoms structure 
# ignore  a list of contacts to ignore, 
#         could be useful for H transfer
# bond_single if there is a lone atom, find a molecule for it
#
# TODO: check PBC
# TODO: get atoms together for PBC (move atoms between unit cells)
# TODO: implement ignore
# TODO: define smallest cluster instead of bond_single
#######################
def monomerise(a, ignore=[], bond_single=False):
    # to use Atoms
    import ase
    # covalent radii needed to check distance
    from ase.data import covalent_radii
    # math for sqrt
    import math   
    
    #make a nparray from positions
    pos=a.get_positions()
    elem=a.get_atomic_numbers()
    #pbc is a [T/F T/F T/] structure
    #if it is false, the get cell will return zero, 
    #so we should be ok with these statements
    pbc=a.get_pbc()
    cell=a.get_cell()
   
    #create connectivity matrix that stores if atoms are 
    #connected or not
    connect=np.zeros((len(elem),len(elem)), dtype=bool)
    connected=np.zeros(len(elem), dtype=int)
 
    #loop over pairs of atoms and set-up the connectivity matrix based on distance
    #and covalent radii
    for i in range(len(elem)):
       for j in range(i+1,len(elem)):
          #get distance considering pbc if necessary 
          #(hopefully works ootb when pbc=False)
          dist=a.get_distance(i,j,mic=True)

          #check if distance is small
          if dist<1.1*(covalent_radii[elem[i]]+covalent_radii[elem[j]]):
             connect[i,j]=True
             connect[j,i]=True
             connected[i]=1
             connected[j]=1
       #we did not find a neighbour for this atom
       #check if bond_single is set and we must find a molecule for this atom
       if connected[i]==0 and bond_single==True: 
          print('will search for a neighbour for atom', i)
          distmin=1000.
          idxmin=-1
          #find an index of atom with the shortest distance to the lone atom i
          for j in range(len(elem)):
             if i!=j: 
                dist=a.get_distance(i,j,mic=True)
                if dist<distmin:
                   distmin=dist
                   idxmin=j
          if idxmin!=-1:
             print('found a neighbour ', idxmin, ' with a distance ', distmin)
             connect[i,idxmin]=True
             connect[idxmin,i]=True
             connected[i]=1
             connected[idxmin]=1
          else:
             print('did not find a neighbour atom for atom ', i)

    print(connect)

    #now find which atoms are connected with each other and remember it
    #create an array which stores to which group/molecule atom belongs
    imol=np.arange(len(elem))

    #loop over atom pairs and sort them to the groups
    #not sure if necessary, remove the print after testing
    for ii in range(len(elem)):
      nchange=0
      for i in range(len(elem)):
        for j in range(i+1,len(elem)):
           if connect[i,j]==True:
              #if they are connected and belong to different group 
              #set them to the same group, use the lower value
              if (imol[j]<imol[i]):
                 imol[i]=imol[j]
                 nchange+=1
              elif (imol[j]>imol[i]):
                 imol[j]=imol[i]
                 nchange+=1
      print(ii, nchange)
      print(imol)
 
    #find how many monomers we have 
    #we use a little trick, or we are lazy, one of that
    #create an array length number of atoms filled with zeros
    nmono_ar=np.zeros(len(elem))  
    #loop over all the atoms and set nmono_ar element of their group to one
    for idx in range(len(elem)):
       nmono_ar[imol[idx]]=1
    #sum the elements of the array to find the number of monomers
    nmono=int(np.sum(nmono_ar))
    print('found ', nmono, ' monomers')

    #the division is now basically done
    #the issue is that we have a discontinuous array of molecular indices
    #so we now get an array with continuous indices
    
    #transformation ?
    trafo=np.full(len(elem),-1)
    lfix=np.full(len(elem),False)
    #copy ?
    imol2=imol.copy()
 
    for idx in range(nmono):
       #we currently don't know on which monomer we are working
       #if icur is nonnegative the value is the molecular index
       icur=-1
       for iat in range(len(elem)):
          #loop over atoms and look for elements with imol2 value 
          #different than -1, easy at the start
          if imol2[iat] != -1:
             #we have an atom which belongs to an unprocessed group
             
             if (icur==-1):
                #declare that we are working on the molecule imol2[iat]
                icur=imol2[iat]
                #what used to be icur (discontinuous) will be now idx (continuous)
                trafo[icur]=idx
                #don't consider this atom anymore
                imol2[iat]=-1
                #don't move this atom across cells
                lfix[iat]=True
             else:
                #icur has some value which corresponds to some molecule index
                #if we have atom with the same index we remove it as it will
                #not form a new molecule
                #(yes, this is a bit convoluted, we could have set it to some group now)
                if imol2[iat]==icur:
                   imol2[iat]=-1

    print('old sorting ', imol)
    #divide the atoms into the continous group now
    for iat in range(len(elem)):
       #this looks safe to do as we work only on one element of imol
       #trafo can be [0,4,5], it gives the lowest atom indices for each monomer
       imol[iat]=trafo[imol[iat]]

    print('new sorting ', imol)

    #create list for storing monomer structures
    co=[]
    for idx in range(nmono):
       co.append(a.copy())

    #loop over monomers and delete the atoms which do not belong to it
    #a bit silly, could have used append I guess
    #or some Atoms create stuff
    for imono in range(nmono):
       one=co[imono]
       del one[[atom.index for atom in one if imol[atom.index] != imono] ]
  
    #return shifted cluster, monomers, transformation and imol
    #last two are needed when dividing other arrays into monomer parts
    return a, co, trafo, imol
 


