#!/usr/bin/env python
import time
import numpy as np
from os.path import isfile
from sys import exit,argv
from scipy import ndimage
from scipy.ndimage.filters import gaussian_filter
from scipy.constants import physical_constants
import argparse,copy
import pdb
from argparse import RawTextHelpFormatter
from skimage import transform
from scipy.spatial.transform import Rotation as R

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
0.3       26/11/2017        A. El-Sayed       Added supercube function. More command line functions.
0.4       19/12/2018        A. El-Sayed       Fixed cube integration around an atom and reference point. Can be called from command line
0.5       12/03/2019        A. El-Sayed       Added cube rotation from the command line
'''
__version__ = 0.3

class cube():
    '''
    Cube Class:
    Includes a bunch of methods to manipulate cube data
    '''

    def __init__(self,fname=None):
        if fname != None:
            try:
                self.read_cube(fname)
            except IOError as e:
                print( "File used as input: %s" % fname )
                print( "File error ({0}): {1}".format(e.errno, e.strerror))
                self.terminate_code()
        else:
            self.default_values()
        return None

    def terminate_code(self):
        print( "Code terminating now")
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
        """
        Method to read cube file. Just needs the filename
        """

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
                self.atomsXYZ.append(list(map(float,[line[2], line[3], line[4]])))
            self.data = np.zeros((self.NX,self.NY,self.NZ))
            i= int(0)
            for s in fin:
                for v in s.split():
                    self.data[int(i/(self.NY*self.NZ)), int((i/self.NZ)%self.NY), int(i%self.NZ)] = float(v)
                    i+=1
            # if i != self.NX*self.NY*self.NZ: raise NameError, "FSCK!"
        return None

    def write_cube(self,fname,comment='Cube file written by CubeToolz\nCubeToolz %3.1f' % __version__):
        '''
        Write out a Gaussian Cube file
        '''
        try:
            with open(fname,'w') as fout:
                if len(comment.split('\n')) != 2:
                    print( 'Comment line NEEDS to be two lines!')
                    self.terminate_code()
                fout.write('%s\n' % comment)
                fout.write("%4d %.6f %.6f %.6f\n" % (self.natoms, self.origin[0], self.origin[1], self.origin[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NX, self.X[0], self.X[1], self.X[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NY, self.Y[0], self.Y[1], self.Y[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NZ, self.Z[0], self.Z[1], self.Z[2]))
                for atom,xyz in zip(self.atoms,self.atomsXYZ):
                    fout.write("%s %d %6.3f %6.3f %6.3f\n" % (atom, 0, xyz[0], xyz[1], xyz[2]))
                for ix in range(self.NX):
                   for iy in range(self.NY):
                       for iz in range(self.NZ):
                           fout.write("%.5e " % self.data[ix,iy,iz]),
                           if (iz % 6 == 5): fout.write('\n')
                       fout.write("\n")
        except IOError as e:
            print( "File used as output does not work: %s" % fname)
            print( "File error ({0}): {1}".format(e.errno, e.strerror))
            self.terminate_code()
        return None

    def square_cube(self,power=2):

        '''
        Function to raise cube data to a power. Squares cube data by default.
        '''
        self.data=self.data**power
        print( power )
        return None

    def rotate_cube(self, angle, axes=None):
        '''
        Rotate cube data around a plane. The plane is defined in the axes variable with origin set as point. For example, to rotate along the xy plane, axes would be defined as (0,1).
        '''
        if 0 not in axes:
            rotAxis = 'x'
        elif 1 not in axes:
            rotAxis = 'y'
        elif 2 not in axes:
            rotAxis = 'z'

        # Rotate the atoms
        r = R.from_euler(rotAxis, angle, degrees=True)
        self.atomsXYZ = r.apply(self.atomsXYZ)

        #Centre of new cell
        centreNewCell = np.sum(r.apply([self.X * self.NX, self.Y * self.NY, self.Z * self.NZ]) / 2., axis = 0)

        # Rotate the cube data
        self.data = ndimage.rotate(self.data, angle, axes=axes, mode='wrap')
        self.NX,self.NY,self.NZ = np.shape(self.data)

        # Move atoms' centre of mass to centre of cell
        newCentre = ((self.X * self.NX + self.Y * self.NY + self.Z * self.NZ)) / 2.
        centreDiff = centreNewCell - newCentre

        self.atomsXYZ = self.atomsXYZ - centreDiff

        # Make sure atoms are in cell
        self.atomsXYZ = self.atomsXYZ % (self.X * self.NX + self.Y * self.NY  + self.Z * self.NZ)

        return None


    def translate_cube(self,tVector):
        '''
        Translate cube data by some vector. The vector is given as a list to the tVector function.
        '''
        self.data = ndimage.shift(self.data,tVector,mode='wrap')
        return None

    def planar_average(self,axis):
        '''
        Calculate the planar average along an axis. The axis is given as a string of either x,y or z.
        '''
        bohrM=physical_constants['Bohr radius'][0]
        bohrA=physical_constants['Bohr radius'][0]*1e10
        if axis == 'x':
            yz_area=np.linalg.norm(np.cross((self.NY*self.Y*bohrM),(self.NZ*bohrM*self.Z)))
            PlanAv=np.array([[nx*self.X[0]*bohrA,np.sum(self.data[nx,::])/(self.NY*self.NZ)] for nx in range(self.NX)])
        elif axis == 'y':
            xy_area=np.linalg.norm(np.cross((self.NX*self.X*bohrM),(self.NZ*bohrM*self.Z)))
            PlanAv=np.array([[ny*self.Y[1]*bohrA,np.sum(self.data[:,ny,:])/(self.NX*self.NZ)] for ny in range(self.NY)])
        elif axis == 'z':
            xy_area=np.linalg.norm(np.cross((self.NY*self.Y*bohrM),(self.NX*bohrM*self.X)))
            PlanAv=np.array([[nz*self.Z[2]*bohrA,np.sum(self.data[nz])/(self.NY*self.NZ)] for nz in range(self.NZ)])
        else:
            print( '%s' % 'No axis specified! Planar average will return zero and fail.')
            PlanAv = 0.0
        return PlanAv

    def planar_averageG(self,axis,sigma):
        '''
        Broaden the planar average along an axis. The axis is given as a string of either x,y or z. A broadening value is also needed.
        '''
        PlanAvG = self.planar_average(axis)
        PlanAvG[:,1] = gaussian_filter(PlanAvG[:,1],sigma)
        return PlanAvG

    def cube_int(self):
        '''
        Integrate the entire cube data.
        '''
        angstrom=physical_constants['Bohr radius'][0]*1e10
        vol=np.linalg.det(np.array([self.X,self.Y,self.Z]))
        edensity=np.sum(self.data)
        nelectron=vol*edensity

        #print 'Number of electrons: %.7g' % (nelectron)
        return nelectron

    def cube_int_atom(self,atomID,radius):
        '''
        Integrate the cube data in a sphere around a particular atom. Needs the atom number (note that atom 0 is the first atom). Also needs a radius of the sphere.
        '''
        nelectron = 0.0
        voxelMatrix = np.array([self.X,self.Y,self.Z])
        radius *= 1 / (physical_constants['Bohr radius'][0] * 1e10)
        vol = np.linalg.det(voxelMatrix)
        atomXYZ = np.array(self.atomsXYZ[atomID]) - self.origin

        initial = time.time()
        for x in range(self.NX):
            for y in range(self.NY):
                for z in range(self.NZ):
                   pos = np.array([x * self.X[0],y * self.Y[1],z * self.Z[2]])
                   distance = np.linalg.norm(pos - atomXYZ)
                   if distance <= radius:
                       nelectron += self.data[x][y][z] * vol
        final=time.time()
        forTime = final - initial
        #print 'for loop: %.2f s' % forTime

        return nelectron


    def cube_int_ref(self,ref,radius):
        '''
        Integrate the cube data in a sphere around a point. Needs the atom number (note that atom 0 is the first atom). Also needs the point as a list.
        '''
        nelectron = 0.0
        voxelMatrix = [self.X,self.Y,self.Z]
        vol = np.linalg.det(voxelMatrix)
        ref = np.array(ref)
        ref = ref * (1 / (physical_constants['Bohr radius'][0] * 1e10))

        initial = time.time()
        for x in range(self.NX):
            for y in range(self.NY):
                for z in range(self.NZ):
                   pos = np.array([x * self.X[0],y * self.Y[1],z * self.Z[2]])
                   distance = np.linalg.norm(pos - ref)
                   if distance <= radius:
                       nelectron += self.data[x][y][z] * vol
        final = time.time()
        forTime = final - initial
        #print 'for loop: %.2f s' % forTime

        return nelectron

    def super_cube(self,new_size):
        '''
        Function to make a new cube supercell. Takes in 3D list of how big the supercell should be.
        '''
        cell = np.array((self.NX * self.X, self.NY * self.Y, self.Z * self.NZ))
        new_data = np.zeros([new_size[0]*self.NX,new_size[1]*self.NY,new_size[2]*self.NZ])
        new_xyz = self.atomsXYZ
        n_newcells = np.prod(new_size)
        new_xyz = np.array(tuple(self.atomsXYZ)*n_newcells)
        counter = 0
        for x in range(new_size[0]):
            for y in range(new_size[1]):
                for z in range(new_size[2]):
                        new_data[x*self.NX:(x+1)*self.NX,y*self.NY:(y+1)*self.NY,z*self.NZ:(z+1)*self.NZ] += self.data
                        new_xyz[counter*self.natoms:(counter+1)*self.natoms] = self.atomsXYZ + ( (np.array(cell[0]) * x ) + (np.array(cell[1]) * y) + (np.array(cell[2]) * z))
                        counter += 1
        new_data = transform.rescale(new_data,1/np.mean(new_size),order=3)
        self.data = new_data
        self.X = ((self.NX * self.X) / float(np.shape(new_data)[0])) * new_size[0]
        self.Y = ((self.NY * self.Y) / float(np.shape(new_data)[1])) * new_size[1]
        self.Z = ((self.NZ * self.Z) / float(np.shape(new_data)[2])) * new_size[2]
        self.NX, self.NY, self.NZ = np.shape(new_data)
        self.atomsXYZ = new_xyz
        self.atoms *= len(self.atomsXYZ)
        self.natoms = len(new_xyz)
        return None


def add_cubes(files):
    cubes = [cube(fin) for fin in files]
    print( "====== Adding cube files ======")
    cube_out = copy.deepcopy(cubes[0])

    for ctmp in cubes[1:]:
        cube_out.data += ctmp.data
    print( "====== Writing output cube as diff.cube ======")
    cube_out.write_cube('diff.cube')
    return cube_out

def diff_cubes(files):
    cubes = [cube(fin) for fin in files]
    print( "====== Subtracting cube files ======")
    cube_out = copy.deepcopy(cubes[0])

    for ctmp in cubes[1:]:
        cube_out.data -= ctmp.data
    print( "====== Writing output cube as diff.cube ======")
    cube_out.write_cube('diff.cube')
    return cube_out

def mult_cubes(files):
    cubes = [cube(fin) for fin in files]
    print( "====== Multiplying cube files ======")
    cube_out = copy.deepcopy(cubes[0])

    for ctmp in cubes[1:]:
        cube_out.data *= ctmp.data
    print( "====== Writing output cube as mult.cube ======")
    cube_out.write_cube('mult.cube')
    return cube_out

def square_cubes(files,power):
    cubes = [cube(fin) for fin in files]
    print( "====== Squaring cube files ======")
    [ctmp.square_cube(power) for ctmp in cubes]

    print( "====== Writing output cubes as squareN.cube ======")
    if len(cubes) == 1:
        cubes[0].write_cube('square.cube')
    else:
        for ind,cout in enumerate(cubes):
            cout.write_cube('square%d.cube' % ind)
    return None

def translate_cubes(files,tVector):
    cubes = [cube(fin) for fin in files]
    print( "====== Squaring cube files ======")
    [ctmp.translate_cube(tVector) for ctmp in cubes]

    print( "====== Writing output cubes as translateN.cube ======")
    if len(cubes) == 1:
        cubes[0].write_cube('translate.cube')
    else:
        for ind,cout in enumerate(cubes):
            cout.write_cube('translate%d.cube' % ind)
    return None

def rotate_cube(files, angle, axes):
    cubes = cube(files)
    print("====== Rotating Cube ======")
    cubes.rotate_cube(angle, axes = axes)
    print("====== Writing output cubes as rotated.cube ======")
    cubes.write_cube('rotated.cube')
    return None

def expand_cell(files,new_size):
    cube_in = cube(files[0])
    cube_in.super_cube(new_size)
    cube_in.write_cube('expand_%dx%dx%d.cube' % (new_size[0],new_size[1],new_size[2]))
    return None

def cube_integrate(files):
    cube_in = cube(files[0])
    cube_int_total = cube_in.cube_int()
    print( 'Integral of cube file is: %9.3f' % cube_int_total)
    return None

def cube_integrate_atom(files, atomId, radius):
    cube_in = cube(files[0])
    atomId = int(atomId)
    radius = float(radius)
    cube_int_total = cube_in.cube_int_atom(atomId, radius)
    print( 'Integral around chosen atom is: %12.9f' % cube_int_total)
    return None

def cube_integrate_ref(files, x, y, z, radius):
    cube_in = cube(files[0])
    xyz = np.array([float(x), float(y), float(z)])
    radius = float(radius)
    cube_int_total = cube_in.cube_int_ref(xyz, radius)
    print( 'Integral around chosen reference point is: %12.9f' % cube_int_total)
    return None

def planar_average_cube(files,vector):
    cube_in = cube(files[0])
    planav = cube_in.planar_average(vector[0])
    return planav

def main():
    parser = argparse.ArgumentParser(description="A python library and tool to read in and manipulate Gaussian cube files. This code allows you to:\n    Read and write Gaussian cube files\n    Translate and rotate cube data\n    Integrate around a particular atom\n    Integrate around a sphere\n    Integrate around the whole cube file",formatter_class=RawTextHelpFormatter)
#    Take the planar average")

    parser.add_argument("Files",help="Cube files used in program",nargs = '+')
    parser.add_argument("-a","--add",help="Add two or more cube files together",action = "store_true")
    parser.add_argument("-s","--subtract",help="Subtract two or more cube files together",action = "store_true")
    parser.add_argument("-M","--multiply",help="Multiply two or more cube files together",action = "store_true")
    parser.add_argument("-p","--power",help="Raise the cube file to a certain power. Any number of cube files can be specified and they will all be raised to the power defined. Default is to square the cube file(s).",nargs='?',const=2,type=int)
    parser.add_argument("-i","--integrate",help="Integrate over the entire cube file.")
    parser.add_argument("-ia","--integrateatom",help="Integrate a sphere around a particular atom. Needs atom id and a radius to integrate within.", nargs=2 )
    parser.add_argument("-ir","--integrateref",help="Integrate a sphere around a reference xyz. Supply x y z r, where r is the radius.", nargs=4 )
    parser.add_argument("-e","--expand",help="Make a supercell of the specified cube file",nargs=3,type=float)
    parser.add_argument("-m","--mean",help="Calculate planar average of a cube file along a particular axis. Arguments are x,y or z.",nargs=1,type=str)
    parser.add_argument("-t","--translate",help="Translate a cube file. Requires a translation vector as an argument.", nargs = 3,type=float)
    parser.add_argument("-r","--rotate",help="Rotate a cube file. Requires an angle and an axis around which to rotate as an argument. The axis is provided with the -ax flag.", nargs = 1,type=float)
    parser.add_argument("-ax","--axis",help="Axis parameter", nargs = 2,type=int)
    if len(argv) <= 2:
        parser.print_help()

    args = parser.parse_args()

    if args.add:
        if len(args.Files) >= 2:
            add_cubes(args.Files)
        else:
            print( "Error: To use the add function, two or more cube files need to be specified.")
    if args.subtract:
        if len(args.Files) >= 2:
            diff_cubes(args.Files)
        else:
            print( "Error: To use the subtract function, two or more cube files need to be specified.")
    if args.multiply:
        if len(args.Files) >= 2:
            mult_cubes(args.Files)
        else:
            print( "Error: To use the multiply function, two or more cube files need to be specified.")
    if args.power:
        if args.Files:
            square_cubes(args.Files,args.power)
        else:
           print( "Error: At least one cube file needed to calculate its square.")
    if args.translate:
        if args.Files:
            translate_cubes(args.Files,args.translate)
    if args.expand:
        if args.Files:
            print( type(args.expand[0]))
            expand_cell(args.Files,list(map(int,args.expand)))
    if args.mean:
        if args.Files:
            PlanAv = planar_average_cube(args.Files,args.mean)
            fout = open('planav.dat','w')
            for line in PlanAv:
                fout.write('%9.4f  %9.4f\n' % (line[0],line[1]))
    if args.integrate:
        if args.Files:
             cube_integrate(args.Files)
    if args.integrateatom:
        if args.Files:
            cube_integrate_atom(args.Files, args.integrateatom[0], args.integrateatom[1])
    if args.integrateref:
        if args.Files:
            cube_integrate_ref(args.Files, args.integrateref[0], args.integrateref[1], args.integrateref[2], args.integrateref[3])
    if args.rotate:
        if args.Files:
            if args.axis:
                axes = (args.axis[0], args.axis[1])
                rotate_cube(args.Files[0], args.rotate[0], axes)
            else:
                print('No axis provided. Please provide an axis using the -ax flag.')

    return None

if __name__ == '__main__':
    main()
