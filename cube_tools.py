#!/usr/bin/env python
import numpy as np 
from os.path import isfile 
from sys import exit
from scipy import ndimage

'''
-----------
Cube Tools
-----------

A python library and tool to read in and manipulate Gaussian cube files. This code allows you to:
    Read and write Gaussian cube files
    Translate and rotate cube data
    Integrate around a particular atom
    Integrate around a sphere
    Integrate around the whole cube file
    Take the planar average



Version      Date              Coder          Changes
=======   ==========        ===========       =======
0.1       07/04/2017        A. El-Sayed       Initial Version
'''


class cube():
    ''' Cube Class'''

    def __init__(self,fname=None):
        if fname != None:
            try:
                self.read_cube(fname)
            except IOError as e:
                print "File used as input: %s" % fname
                print "File error ({0}): {1}".format(e.errno, e.strerror)
                self.terminate_code()
        else:
            self.default_values()

    def terminate_code(self):
        print "Code terminating now"
        exit()

    def default_values(self):
            self.natoms=0
            self.comment1=0
            self.comment2=0
            self.origin=np.array([0,0,0])
            self.nX=0
            self.nY=0
            self.nZ=0
            self.X=0
            self.Y=0
            self.Z=0
            self.atoms=['0']
            self.atomsXYZ=[0,0,0]
            self.data=[0]


    def read_cube(self,fname):
        """Class docstring."""
        fin = open(fname, 'r')
        self.filename = fname
        self.comment1 = fin.readline() #Save 1st comment
        self.comment2 = fin.readline() #Save 2nd comment
        nOrigin = fin.readline().split() # Number of Atoms and Origin 
        self.natoms = int(nOrigin[0]) #Number of Atoms
        self.origin = np.array([float(nOrigin[1]),float(nOrigin[2]),float(nOrigin[3])]) #Position of Origin
        nVoxel = fin.readline().split() #Number of Voxels
        self.NX = int(nVoxel[0])
        self.X = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
        nVoxel = fin.readline().split() #
        self.NY = int(nVoxel[0])
        self.Y = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
        nVoxel = fin.readline().split() #
        self.NZ = int(nVoxel[0])
        self.Z = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
        self.atoms = []
        self.atomsXYZ = []
        for atom in range(self.natoms):
            line= fin.readline().split()
            self.atoms.append(line[0])
            self.atomsXYZ.append(map(float,[line[2], line[3], line[4]]))
        self.data = np.zeros((self.NX,self.NY,self.NZ))
        i=0
        for s in fin:
            for v in s.split():
               self.data[i/(self.NY*self.NZ), (i/self.NZ)%self.NY, i%self.NZ] = float(v)
               i+=1
        if i != self.NX*self.NY*self.NZ: raise NameError, "FSCK!"

    def write_cube(self,fname,comment='Cube file written by CubeToolz\nCubeToolz 0.1'):
        '''
        Output a Gaussian Cube file
        '''
        try:
            fout=open(fname,'w')
        except IOError as e:
            print "File used as output does not work: %s" % fname
            print "File error ({0}): {1}".format(e.errno, e.strerror)
            self.terminate_code()
        if len(comment.split('\n')) != 2:
            print 'Comment line NEEDS to be two lines!'
            self.terminate_code()
        fout.write('%s\n' % comment)
        fout.write("%4d %.6f %.6f %.6f\n" % (self.natoms, self.origin[0], self.origin[1], self.origin[2]))
        fout.write("%4d %.6f %.6f %.6f\n" % (self.NX, self.X[0], self.X[1], self.X[2]))
        fout.write("%4d %.6f %.6f %.6f\n" % (self.NY, self.Y[0], self.Y[1], self.Y[2]))
        fout.write("%4d %.6f %.6f %.6f\n" % (self.NZ, self.Z[0], self.Z[1], self.Z[2]))
        for atom,xyz in zip(self.atoms,self.atomsXYZ):
            fout.write("%s %d %6.3f %6.3f %6.3f\n" % (atom, 0, xyz[0], xyz[1], xyz[2]))
        for ix in xrange(self.NX):
           for iy in xrange(self.NY):
               for iz in xrange(self.NZ):
                   fout.write("%.5e " % self.data[ix,iy,iz]),
                   if (iz % 6 == 5): fout.write('\n')
               fout.write("\n")
        fout.close()

    def square_cube(self):
        '''
        Some bullshit
        '''
        self.data=self.data**2
        return None

    def rotate_cube(self,angle,axes=(0,1)):
        '''
        Some bullshit
        '''
        self.data=ndimage.rotate(self.data,angle,axes=axes,mode='wrap')
        return None

    def translate_cube(self,


Pooo



wweee


Whyyy what is happeninng
