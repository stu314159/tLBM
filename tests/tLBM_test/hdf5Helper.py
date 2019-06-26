import numpy as np
import h5py

def writeH5(pressure,u,v,w,velmag,filename):
  """
  Write the h5 file that will save the information needed in proper structure.
  pressure = numpy array with pressure values
  u,v,w = numpy array with velocity data
  filename = string with desired filename
  dims = 3-tuple with the number of rank of each dimension
  """

  f = h5py.File(filename,'w')

  # Store velocity data into the velo_group of h5 file
  velo_group = f.create_group("velo_group")
  x_velo = velo_group.create_dataset("x_velo",data=u)
  y_velo = velo_group.create_dataset("y_velo",data=v)
  z_velo = velo_group.create_dataset("z_velo",data=w)
  velmag = velo_group.create_dataset("velmag",data=velmag)

  # Store velocity data into the velo_group of h5 file
  pres_group = f.create_group("pres_group")
  presmag = pres_group.create_dataset("presmag",data=pressure)

  f.close()

def writeXdmf(dims,dx,filename,h5_file):
  """
  Write the xmf file, that describes the hdf5 data, to be read by Paraview.
  filename = string with the desired filename
  dims = 3-tuple with the number of rank in each dimension (z,y,x)
  """

  f = open(filename,'w')
  f.write('<?xml version="1.0" ?>\n')
  f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
  f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.1">\n')
  f.write('<Domain>\n')

  f.write('<Grid Name="my_Grid" GridType="Uniform">\n')
  f.write('<Topology TopologyType="3DCoRectMesh" Dimensions="%d %d %d">\n'%(dims[0],dims[1],dims[2]))
  f.write('</Topology>\n')

  f.write('<Geometry GeometryType="Origin_DxDyDz">\n')
  f.write('<DataItem Dimensions="3" NumberType="Integer" Format="XML">\n')
  f.write('0 0 0\n') 
  f.write('</DataItem>\n')
  f.write('<DataItem Dimensions="3" NumberType="Integer" Format="XML">\n')
  f.write('%g %g %g\n'%(dx,dx,dx))
  f.write('</DataItem>\n')
  f.write('</Geometry>\n')

  f.write('<Attribute Name="velocity" AttributeType="Vector" Center="Node">\n')
  f.write('<DataItem ItemType="Function" Function="JOIN($0, $1, $2)" Dimensions="%d %d %d 3">\n'%(dims[0],dims[1],dims[2]))
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Float" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/x_velo\n')
  f.write('%s:/velo_group/x_velo\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Float" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/y_velo\n')
  f.write('%s:/velo_group/y_velo\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Float" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/z_velo\n')
  f.write('%s:/velo_group/z_velo\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('</DataItem>\n')
  f.write('</Attribute>\n')

  f.write('<Attribute Name="pressure" AttributeType="Scalar" Center="Node">\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Float" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/pres_group/presmag\n')
  f.write('%s:/pres_group/presmag\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('</Attribute>\n')

  f.write('<Attribute Name="velocityMagnitude" AttributeType="Scalar" Center="Node">\n')
  f.write('<DataItem Dimensions="%d %d %d" NumberType="Float" Format="HDF">\n'%(dims[0],dims[1],dims[2]))
  #f.write('out'+str(i)+'.h5:/velo_group/velmag\n')
  f.write('%s:/velo_group/velmag\n'%h5_file)
  f.write('</DataItem>\n')
  f.write('</Attribute>\n')

  f.write('</Grid>\n')
  f.write('</Domain>\n')
  f.write('</Xdmf>\n')

  f.close()
