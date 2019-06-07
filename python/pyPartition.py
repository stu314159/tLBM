#!/usr/bin/env python3
#pyPartition.py
"""
provide class libraries for partitioning
goal of this code is to provide a reasonable geometric partition
to the pre-processing libraries.  The output will be a list of length Nx*Ny*Nz
containing the integer of which partition each lattice point lives.

"""

import sys
sys.path.insert(1,'.')

import partition_suggestion as ps
import partition_compare as pc
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK

NO_PYMETIS=0
try:
  from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected
except ImportError:
  NO_PYMETIS=1


#import PartitionHelper as PH
import numpy as np
import sys

class Lattice(object):
    """
       define the layout and adjacency of the LBM lattice
    """
    def __init__(self,Nx,Ny,Nz):
        """
            basic constructor
            Nx - number of lattice points in the x-direction
            Ny - number of lattice points in the y-direction
            Nz - number of lattice points in the z-direction
            
        """
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz
        self.ex = []; self.ey = []; self.ez = []
        self.bbSpd = []; self.w = []; 
        #self.initialize_adjDict(); # just do it.
        self.adjDict = None
        self.cutSize = None  # must ask for cut size after partitioning
        self.partition = None # must give extra inputs for the partition

    def get_dims(self):
        return [self.Nx, self.Ny, self.Nz]

    def get_ex(self):
        return self.ex[:]

    def get_ey(self):
        return self.ey[:]

    def get_ez(self):
        return self.ez[:]

    def get_numSpd(self):
        return len(self.ex)

    def get_bbSpd(self):
        return self.bbSpd[:]

    def get_w(self):
        return self.w[:]

    def get_nnodes(self):
        return self.Nx*self.Ny*self.Nz

    def initialize_adjDict(self):  
        
        self.adjDict = pc.set_adjacency(self.Nx,self.Ny,self.Nz,self.ex,self.ey,self.ez)
        #self.partHelper = PH.PartitionHelper(self.Nx,self.Ny,self.Nz,self.get_numSpd());
        #self.adjDict = {}
        #self.partHelper.setAdjacency(self.adjDict); # now self.adjDict is populated
        # hopefully this took less time than before.
        


    def compute_cutSize(self):
        if self.adjDict == None: # verify that the adjacency list has been initialized
            raise ValueError('adjacency list must be initialized before getting cut size')
        else:
            self.cutSize = pc.count_cuts(self.adjDict,self.partition.get_partition()) #yeah -- looks overly complicated

        return self.cutSize

    def get_cutSize(self):
        return self.cutSize

    def set_Partition(self, numParts = 1, numTrials = 2000, style = '1D'):
        """
          numParts = number of partitions
          numTrials = number of random 3D partition permutations should be tested
          style = ['1D', '3D','metis']

        """
        self.partition = Partitioner(self.Nx, self.Ny, self.Nz, 
                                     numParts = numParts, adjList = self.adjDict,
                                     numTrials = numTrials,
                                     style = style)

    

class D3Q15Lattice(Lattice):
    """
      D3Q15 Lattice
    """
    def __init__(self,Nx,Ny,Nz):
        """
          D3Q15 Lattice
        """
        super(D3Q15Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1]
        self.ey = [0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1]
        self.ez = [0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1]
        self.bbSpd = [0,2,1,4,3,6,5,14,13,12,11,10,9,8,7]
        self.w = [2./9.,1./9.,1./9,1./9.,1./9.,1./9.,1./9.,
	       1./72.,1./72.,1./72.,1./72.,
	       1./72.,1./72.,1./72.,1./72.]

class D3Q19Lattice(Lattice):
    """
    """
    def __init__(self,Nx,Ny,Nz):
        super(D3Q19Lattice,self).__init__(Nx,Ny,Nz)
        self.ex =  [0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0]
        self.ey = [0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1]
        self.ez = [0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1]
        self.bbSpd = [0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15]
        self.w = [2./9.,1./9.,1./9,1./9.,1./9.,1./9.,1./9.,
	       1./72.,1./72.,1./72.,1./72.,
	       1./72.,1./72.,1./72.,1./72.]

   
class D3Q27Lattice(Lattice):
    """
    """
    def __init__(self,Nx,Ny,Nz):
        super(D3Q27Lattice,self).__init__(Nx,Ny,Nz)
        self.ex = [0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1]; self.ex=np.array(self.ex,dtype=np.float32)
        self.ey = [0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1]; self.ey=np.array(self.ey,dtype=np.float32)
        self.ez = [0,0,0,0,0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]; self.ez=np.array(self.ez,dtype=np.float32)
        self.bbSpd = [0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15,26,25,24,23,22,21,20,19]
        self.w = [8./27.,2./27.,2./27.,2./27.,2./27.,2./27.,2./27.,
                  1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
                  1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
                  1./216.,1./216.,1./216.,1./216.,
                  1./216.,1./216.,1./216.,1./216.]

class Partitioner:
    """
     the class that will do the work to select and obtain a partition
    """

    def __init__(self,Nx,Ny,Nz,numParts,adjList,numTrials = 2000, style = '1D'):
        """
          Nx - number of lattice points in the x-direction (int)
          Ny - number of lattice points in the y-direction (int)
          Nz - number of lattice points in the z-direction (int)
          numParts - number of partitions to form (int)
          adjList - adjacency list (dictionary)
          numTrials - number of attempts that the randomized
                      partition advisor should use to find 
                      a good partitioning
          style - '1D', '3D','metis' partition style
                      
        """
        self.Nx = Nx; self.Ny = Ny; self.Nz = Nz
        self.numParts = numParts
        self.numTrials = numTrials
        self.style = style
        self.adjList = adjList

        if style=='1D':
            self.px = 1; self.py = 1; self.pz = self.numParts;
        elif style=='3D':
            [self.px,self.py,self.pz] = ps.part_advisor(self.Nx,self.Ny,self.Nz,
                                                        self.numParts,
                                                        self.numTrials)
                    
        if (style == '1D' or style == '3D'):  
            self.part_vert = pc.set_geometric_partition(self.Nx, self.Ny, self.Nz,
                                                self.px, self.py, self.pz)
        else:
          if (NO_PYMETIS==1):
            print("pymetis partitioning selected but not available")
            sys.exit()
          [cuts,  self.part_vert] = part_graph(self.numParts,self.adjList) 



    def get_partition(self):
        """
          give access to partition
        """

        return self.part_vert[:]

    
    def get_partition_sizes(self):
        """
          give access to partition sizes
        """
        return [self.px, self.py, self.pz]

    def write_vtk(self,file_name = 'partition_pa.vtk'):
        """
          write out a vtk file to allow visualization of the partitioning
        """
        dims = [self.Nx, self.Ny, self.Nz]
        origin = [0., 0., 0.]
        spacing = [0.1, 0.1, 0.1] #<-- no need for this to correspond to actual physical spacing
        writeVTK(self.part_vert,'partitions',file_name,dims,origin,spacing)


    def write_partition(self):
        """
         write the partition information to parts.lbm
        """
        
#        np_pv = np.array(self.part_vert,dtype=np.int32);
#        fn = 'parts.lbm';
#        np_pv.astype('int32').tofile(fn)
        parts = open('parts.lbm','w')
        for p in self.part_vert:
            parts.write('%d \n'% p)

        parts.close()
   
        


if __name__=="__main__":
    """
      put testing code here  
    """
    #Nx = 447; Ny = 447; Nz = 838;from pymetis import part_graph #<-- requires that the PrgEnv-intel module be selected
    Nx = 10; Ny = 10; Nz = 10;
    print("nnodes = %g" % (Nx*Ny*Nz))
    lat15 = D3Q15Lattice(Nx, Ny, Nz);
    print("initializing the adjacency list")
    lat15.initialize_adjDict();
    numProcs = 4;
    print("setting the partition")
    lat15.set_Partition(numParts = numProcs,numTrials = 8000,style = '3D')

    print("compute cut size")
    lat15.compute_cutSize()
    print("cut size = %g" % lat15.get_cutSize())

    print("writing vtk file for partition")
    lat15.partition.write_vtk('partition_pa.vtk')

    print("re-set partition for 1D geometric")
    lat15.set_Partition(numParts = numProcs,style='1D')
    print("getting new cut size")
    lat15.compute_cutSize()
    print("cut size for 1D geometric partition = %g" % lat15.get_cutSize())

    print("writing vtk file for 1D partition")
    lat15.partition.write_vtk('partition_1D.vtk')

    print("re-set partition for metis")
    lat15.set_Partition(numParts= numProcs, style = 'metis')
    lat15.compute_cutSize()
    print("cut size for metis partition = %g" % lat15.get_cutSize())
    print("writing vtk file for metis partition")
    lat15.partition.write_vtk('partition_metis.vtk')

    print("writing metis partition to disk")
    lat15.partition.write_partition()
    
