#!/usr/bin/env python
import time
import numpy as np 
from os.path import isfile 
from sys import exit
from scipy import ndimage
from scipy.ndimage.filters import gaussian_filter
from scipy.constants import physical_constants
import argparse,copy
import pdb
from argparse import RawTextHelpFormatter

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
    Take the planar average (with an option to Gaussian broaden)



Version      Date              Coder          Changes
=======   ==========        ===========       =======
0.1       07/04/2017        A. El-Sayed       Initial Version
0.2       13/07/2017        A. El-Sayed       Option to run as main program or use as library. All functions written in
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
        return None

    def terminate_code(self):
        print "Code terminating now"
        exit()
        return None

    def default_values(self):
        self.natoms=0
        self.comment1=0
        self.comment2=0
        self.origin=np.array([0,0,0])
        self.NX=0
        self.NY=0
        self.NZ=0
        self.X=0
        self.Y=0
        self.Z=0
        self.atoms=['0']
        self.atomsXYZ=[0,0,0]
        self.data=[0]
        return None


    def read_cube(self,fname):
        """Class docstring."""
        with open(fname, 'r') as fin:
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
        return None

    def write_cube(self,fname,comment='Cube file written by CubeToolz\nCubeToolz 0.1'):
        '''
        Output a Gaussian Cube file
        '''
        try:
            with open(fname,'w') as fout:
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
        except IOError as e:
            print "File used as output does not work: %s" % fname
            print "File error ({0}): {1}".format(e.errno, e.strerror)
            self.terminate_code()
        return None

    def square_cube(self):
        '''
        Some bullshit
        '''
        self.data=self.data**2
        return None



    def rotate_cube(self,angle,axes=None):
        '''
        Some bullshit
        '''
        self.data = ndimage.rotate(self.data,angle,axes=axes,mode='wrap')
        self.NX,self.NY,self.NZ = np.shape(self.data)
        return None

    def rotate_atoms(self, angle, v=np.array([0,0,1]), rotate_cell=False):
        angle *= np.pi / 180
        v /= np.linalg.norm(v)
        if isinstance(v,list):
            v = np.array(v)
        c = np.cos(angle)
        s = np.sin(angle)
        center = np.array([0,0,0])
        center = np.array([np.max(np.array(self.atomsXYZ)[:,0])/2.,np.max(np.array(self.atomsXYZ)[:,1])/2.,np.max(np.array(self.atomsXYZ)[:,2])/.2])

        p = self.atomsXYZ - center
        self.atomsXYZ[:] = (c * p - np.cross(p, s * v) +
                                       np.outer(np.dot(p, v), (1.0 - c) * v) +
                                       center)

        if rotate_cell:
            rotcell = np.array([self.X,self.Y,self.Z]) 
            rotcell[:] = (c * rotcell -
                          np.cross(rotcell, s * v) +
                          np.outer(np.dot(rotcell, v), (1.0 - c) * v))
            #self.set_cell(rotcell)
            self.X,self.Y,self.Z = rotcell[0],rotcell[1],rotcell[2]
        return None

    def rotate(self,angle,axis):
        '''
        Rotates atoms and cube file together.
        Need to rotate atoms then put them in center of cell! Remember that
        '''

    def translate_cube(self,tVector):
        self.data = ndimage.shift(self.data,tVector,mode='wrap')
        return None

    def planar_average(self,axis):
        bohrM=physical_constants['Bohr radius'][0]
        bohrA=physical_constants['Bohr radius'][0]*1e10
        if axis == 'x':
            yz_area=np.linalg.norm(np.cross((self.NY*self.Y*bohrM),(self.NZ*bohrM*self.Z)))
            PlanAv=np.array([[nx*self.X[0]*bohrA,np.sum(self.data[nx,::])/(self.NY*self.NZ)] for nx in xrange(self.NX)])
        elif axis == 'y':
            xy_area=np.linalg.norm(np.cross((self.NX*self.X*bohrM),(self.NZ*bohrM*self.Z)))
            PlanAv=np.array([[ny*self.Y[1]*bohrA,np.sum(self.data[:,ny,:])/(self.NX*self.NZ)] for ny in xrange(self.NY)])
        elif axis == 'z':
            xy_area=np.linalg.norm(np.cross((self.NY*self.Y*bohrM),(self.NX*bohrM*self.X)))
            PlanAv=np.array([[nz*self.Z[2]*bohrA,np.sum(self.data[nz])/(self.NY*self.NZ)] for nz in xrange(self.NZ)])
        else:
            print '%s' % 'No axis specified! Planar average will return zero and fail.'
            PlanAv = 0.0
        return PlanAv

    def planar_averageG(self,axis,sigma):
        PlanAvG = self.planar_average(axis)
        PlanAvG[:,1] = gaussian_filter(PlanAvG[:,1],sigma)
        return PlanAvG

    def cube_int(self):
        angstrom=physical_constants['Bohr radius'][0]*1e10
        vol=np.linalg.det(np.array([self.X,self.Y,self.Z]))
        edensity=np.sum(self.data)
        nelectron=vol*edensity
        
        print 'Number of electrons: %.7g' % (nelectron)
        return nelectron

    def cube_int_atom(self,atomID,radius):
        nelectron = 0.0
        voxelMatrix = [self.X,self.Y,self.Z]
        vol = np.linalg.det(voxelMatrix)
        atomXYZ = self.atomsXYZ[atomID]

        initial = time.time()
        for x in xrange(self.NX):
            for y in xrange(self.NY):
                for z in xrange(self.NZ):
                   xPos,yPos,zPos = [x * self.X[0] / 1.88,y * self.Y[1] / 1.88,z * self.Z[2] / 1.88]
                   distance = np.sqrt((xPos - atomXYZ[0])**2 + (yPos - atomXYZ[1])**2 +(zPos - atomXYZ[2])**2)
                   if distance <= radius:    nelectron += self.data[x][y][z] * vol
        final=time.time()
        forTime = final - initial
        print 'for loop: %.2f s' % forTime
   
        return nelectron

    def cube_int_ref(self,ref,radius):
        nelectron = 0.0
        ind = 0
        voxelMatrix = [self.X,self.Y,self.Z]
        vol = np.linalg.det(voxelMatrix)

        initial = time.time()
        for x in xrange(self.NX):
            for y in xrange(self.NY):
                for z in xrange(self.NZ):
                   xPos,yPos,zPos = [x * self.X[0] / 1.88,y * self.Y[1] / 1.88,z * self.Z[2] / 1.88]
                   distance = np.sqrt((xPos - self.atomsXYZ[atomID][0])**2 + (yPos - self.atomsXYZ[atomID][1])**2 +(zPos - self.atomsXYZ[atomID][2])**2)
                   if distance <= radius:
                       nelectron += self.data[x][y][z] * vol
                   ind += 1
        final=time.time()
        forTime = final - initial
        print 'for loop: %.2f s' % forTime
   
        return nelectron

    def time_loop(self):
        initial=time.time()
        for x in xrange(self.NX):
            for y in xrange(self.NY):
                for z in xrange(self.NZ):
                    a=self.data[x][y][z]
        final=time.time()
        rtime = final - initial
        print 'for loop: %.2f s' % rtime
       
        initial=time.time()
        nelec = 0.0
        for ind,point in enumerate(np.nditer(self.data)):
            nelec += point
        final=time.time()
        rtime = final - initial
        print 'NP iter: %.2f s' % rtime
        print nelec

        initial=time.time()
        vol=np.linalg.det(np.array([self.X,self.Y,self.Z]))
        it = np.nditer(self.data, flags=['multi_index'])
        nelectron = 0.0
        radius = 100
        atomID = self.atomsXYZ[100] 
        while not it.finished:
            nelec += it[0]
            xPos = it.multi_index[0] * self.X[0] / 1.88
            yPos = it.multi_index[1] * self.Y[1] / 1.88
            zPos = it.multi_index[2] * self.Z[2] / 1.88
            distance = 8 #np.sqrt((xPos - atomID[0])**2 + (yPos - atomID[1])**2 +(zPos - atomID[2])**2)
            if distance <= radius:
                nelectron += float(it[0])#points[ind] 
            it.iternext()
        final=time.time()
        rtime = final - initial
        print nelectron*vol
        print 'NP iter: %.2f s' % rtime

        return None

def add_cubes(files):
    cubes = [cube(fin) for fin in files]
    print "====== Adding cube files ======"
    cube_out = copy.deepcopy(cubes[0])
    
    for ctmp in cubes[1:]:
        cube_out.data += ctmp.data 
    print "====== Writing output cube as diff.cube ======"
    cube_out.write_cube('diff.cube') 
    return cube_out

def square_cubes(files):
    cubes = [cube(fin) for fin in files]
    print "====== Squaring cube files ======"
    [ctmp.square_cube() for ctmp in cubes] 
    
    print "====== Writing output cubes as squareN.cube ======"
    if len(cubes) == 1:
        cubes[0].write_cube('square.cube') 
    else:
        for ind,cout in enumerate(cubes):
            cout.write_cube('square%d.cube' % ind) 
    return None 

def translate_cubes(files,tVector):
    cubes = [cube(fin) for fin in files]
    print "====== Squaring cube files ======"
    [ctmp.translate_cube(tVector) for ctmp in cubes] 
    
    print "====== Writing output cubes as translateN.cube ======"
    if len(cubes) == 1:
        cubes[0].write_cube('translate.cube') 
    else:
        for ind,cout in enumerate(cubes):
            cout.write_cube('translate%d.cube' % ind) 
    return None 



def main():
    parser = argparse.ArgumentParser(description="A python library and tool to read in and manipulate Gaussian cube files. This code allows you to:\n    Read and write Gaussian cube files\n    Translate and rotate cube data\n    Integrate around a particular atom\n    Integrate around a sphere\n    Integrate around the whole cube file",formatter_class=RawTextHelpFormatter)
#    Take the planar average")

    parser.add_argument("Files",help="Cube files used in program",nargs = '+')
    parser.add_argument("-a","--add",help="Add two or more cube files together",action = "store_true")
    parser.add_argument("-s","--square",help="Square a cube file",action = "store_true")
    parser.add_argument("-t","--translate",help="Translate a cube file. Requires a translation vector as an argument.", nargs = 3,type=float)
    args = parser.parse_args()

    if args.add:
        if len(args.Files) >= 2: 
            add_cubes(args.Files)
        else:
            print "Error: To use the add function, two or more cube files need to be specified."
    if args.square:
        if args.Files:
            square_cubes(args.Files)
        else:
           print "Error: At least one cube file needed to calculate its square."
    if args.translate:
        if args.Files:
            translate_cubes(args.Files,args.translate) 
        
    return None

if __name__ == '__main__':
    None#main()
