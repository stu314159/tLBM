#FluidChannel.py
"""
Class implementation file for the Python class FluidChannel
Depends on vtkHelper module for geometry visualization functionality
"""
import math
import numpy as np
from vtkHelper import saveStructuredPointsVTK_ascii as writeVTK
import scipy.io


class LatticeSubset(object):
    """
    object that will represent a subset of lattice points.
    
    A fluid channel may have lattice subsets.  Each lattice
    subset will be responsible for identifying its members
    and providing facilities for them.  These facilities
    can be specialized in derived classes.
    """
    def __init__(self):
        """
        
        """
        
    def get_members(self,X,Y,Z):
        """
        given X,Y, and Z coordinates of a lattice
        return the indices of lattice points
        that are a member of the subset.
        """
        
        return []


class YZ_Slice(LatticeSubset):
    """
    subset of a lattice comprising a 2D slice in the YZ plane
    
    """
    
    def __init__(self,Xo,Ymin,Ymax,Zmin,Zmax):
        """
        pass in X-position and YZ bounds for slice
        """
        super(YZ_Slice,self).__init__();
        self.Xo = Xo
        self.Ymin = Ymin;
        self.Ymax = Ymax;
        self.Zmin = Zmin;
        self.Zmax = Zmax;
        
    def get_members(self,X,Y,Z):
        """
        return indices of lattice points that are within dx/2 of specified YZ plane
        """
        x = np.array(X);
        y = np.array(Y);
        z = np.array(Z);
        
        # get lattice spacing.
        xVals = np.unique(x); #sorted unique values
        dx = xVals[1] - xVals[0]; #difference between two sorted unique values
        
        obst = np.where(y >= self.Ymin);
        obst = np.intersect1d(obst[:],np.where(y<=self.Ymax));
        obst = np.intersect1d(obst[:],np.where(z>=self.Zmin));
        obst = np.intersect1d(obst[:],np.where(z<=self.Zmax));
        obst = np.intersect1d(obst[:],np.where(x>=(self.Xo - dx/2.)));
        obst = np.intersect1d(obst[:],np.where(x<=(self.Xo+dx/2.)));
        
        obst = obst.astype(np.int)
        return obst[:]
        
    



class EmptyChannel:  
    """
     a channel with nothing in it
    """
    def __init__(self,Lo):
        """
         constructor
        """
        self.Lo = Lo


    def get_Lo(self):
        """
         set Lo if need be ?
        """
        return self.Lo

    def get_obstList(self,X,Y,Z):
        """
          for an empty channel - no obstacles 
        """
        return []

class ChannelCavity(EmptyChannel):
   """
     a channel with a cavity part way down the length
      - implicitly the channel begins at z = 0
      - and the depth of the channel indicates the maximum value
      - of Y that will be defined as "floor"  
      - the bottom of the cavity will be defined as Y = 0
      - this version of the cavity will span the width (X-direction)
      - of the channel.
   """

   def __init__(self,depth,z_start,z_end):
       self.depth = depth
       self.z_start = z_start
       self.z_end = z_end

   def get_Lo(self):
       return self.depth

   def get_obstList(self,X,Y,Z):
      """ return a list of indices within the boundary of the channel floor
      """
      #x = np.array(X); 
      y = np.array(Y); z = np.array(Z);
      cav1 = np.where(z >= self.z_start)
      cav2 = np.where(z <= self.z_end)
      ol = np.setxor1d(cav1[:],cav2[:])
      cav3 = np.where(y <= self.depth)
      ol = np.intersect1d(ol[:],cav3[:])
      return ol[:]

class GridObst(EmptyChannel):
    """
    a channel with a turbulence-inducing grid at the entrance
    of the channel.  The grid has a circular hole to accomodate
    passage of a drive-shaft for a propeller
    
    """
    def __init__(self,gridZ,xT,yT,zT,xPitch,yPitch,hX,hY,hD):
        """
        gridZ - z-coordinate of center of grid
        xT - thickness of grid webs in x-direction (vertical webs)
        yT - thickness of grid webs in y-direction (horizontal webs)
        zT - thickness of grid in z-direction
        yPitch - pitch of horizontal grids in Y-direction
        xPitch - pitch of vertical grids in the X-direction
        hX - x-coordinate of grid hole
        hY - y-coordinate of grid hole
        hD - diameter of grid hole
        """
        
        self.gridZ = gridZ;
        self.xT = xT;
        self.yT = yT;
        self.zT = zT;
        self.yPitch = yPitch;
        self.xPitch = xPitch;
        self.hX = hX;
        self.hY = hY;
        self.hD = hD;
        
    def get_Lo(self,):
        """
        to do: talk to Luksa about non-dimensionalization
        and most appropriate choice for this
        non-dimensionalization.
        """
        return self.yPitch
    
    def get_obstList(self,X,Y,Z):
        """
        
        """
        x = np.array(X);
        y = np.array(Y);
        z = np.array(Z);
        
        xMax = np.max(x);
        xMin = np.min(x); # expect this to be zero
        yMax = np.max(y);
        yMin = np.min(y); # expect this to be zero
        xPitch = self.xPitch
        yPitch = self.yPitch
        xT = self.xT
        yT = self.yT
        zT = self.zT
        gridZ = self.gridZ
        hX = self.hX
        hY = self.hY
        hD = self.hD
        
        
        obst_list = [];
        # get x-center of vertical grids
        xC_vGrids = np.linspace(xMin+xT/2.,xMax-xT/2.,((xMax-xMin)/xPitch)+1);
        for i in range(len(xC_vGrids)):
            distX = np.abs(x - xC_vGrids[i]);
            distZ = np.abs(z - gridZ);
            gridObstA = np.where((distX < xT/2.))
            gridObstB = np.where((distZ < zT/2.))
            gridObst = np.intersect1d(gridObstA,gridObstB);
            obst_list = np.union1d(obst_list[:],gridObst)
        
        
        # get y-center of horizontal grids
        yC_hGrids = np.linspace(yMin+yT/2.,yMax-yT/2.,((yMax - yMin)/yPitch)+1);
        for i in range(len(yC_hGrids)):
            distY = np.abs(y - yC_hGrids[i]);
            distZ = np.abs(z - gridZ);
            gridObstA = np.where((distY < yT/2.))
            gridObstB = np.where((distZ < zT/2.))
            gridObst = np.intersect1d(gridObstA,gridObstB);
            obst_list = np.union1d(obst_list[:],gridObst)
        
        # remove grids within the hole region
        distH = np.sqrt((y - hY)**2. + (x - hX)**2.)
        obstH = np.where(distH < hD/2.)
        obst_list = np.setdiff1d(obst_list[:],obstH)
            
        obst_list = obst_list.astype(np.int)
        return obst_list[:]
        
        

class StraightPipe(EmptyChannel):
    """
    a square channel where the non-solid nodes
    constitute a circular pipe.
    
    """
    def __init__(self,x_c,y_c,diameter):
        self.x_c = x_c;
        self.y_c = y_c;
        self.diameter = diameter;
        
    def get_Lo(self):
        return self.diameter;
        
    def get_obstList(self,X,Y,Z):
        x = np.array(X); y = np.array(Y); 
        dist = (x - self.x_c)**2 + (y - self.y_c)**2
        
        return list(np.where(dist >= (self.diameter/2.0)**2))
    
class SphereObstruction(EmptyChannel):
    """
     a channel with a sphere obstruction
    """

    def __init__(self,r,x_c,y_c,z_c):
        """
          just need to define the radius and position of the center of the obstacle.
          it is up to the caller to verify that the object will fit within the intended
          channel.  If it does not fit, the obstacle will effectively be
          truncated at the channel boundaries
         
        """
        self.r = r
        self.x_c = x_c
        self.y_c = y_c
        self.z_c = z_c
 
    def get_Lo(self):
        return self.r*2.

    def get_obstList(self,X,Y,Z):
        """
            return a list of all indices all indices within boundary of sphere 
        """

        x = np.array(X); y = np.array(Y); z = np.array(Z);
        dist = (x - self.x_c)**2 + (y - self.y_c)**2 + (z - self.z_c)**2
       
        return list(np.where(dist < self.r**2))
    

class WallMountedBrick(EmptyChannel):
    """
    a channel with a brick mounted to the wall (y = min)
    """
    
    def __init__(self,x_c,z_c,L,W,H):
        """
        x_c = x-coordinate of the brick centroid
        z_c = z-coordinate of the brick centroid
        L = length of the brick (in the Z-direction)
        W = width of the brick (in the X-direction)
        H = height of the brick (in the Y-direction)
        """
        self.x_c = x_c
        self.z_c = z_c
        self.L = L
        self.W = W
        self.H = H
        
    def get_Lo(self):
        
        return self.H
    
    def get_obstList(self,X,Y,Z):
        """
        return a list of all indices within the boundary of the brick
        """
        x = np.array(X); y = np.array(Y); z = np.array(Z);
        
        inH = np.where(y<=self.H);
        inZa = np.where(z>=(self.z_c - self.L/2.));
        inZb = np.where(z<=(self.z_c + self.L/2.));
        inZ = np.intersect1d(inZa,inZb);
        inXa = np.where(x>=(self.x_c - self.W/2.));
        inXb = np.where(x<=(self.x_c + self.W/2.));
        inX = np.intersect1d(inXa,inXb);
        obst = np.intersect1d(inH,inZ);
        obst = np.intersect1d(obst[:],inX);
        
        return obst[:]
        


class TwinJet(EmptyChannel):
    """
    a channel with block obstructions allowing for simulation of 
    the twin-jet problem
  
    """
  
    def __init__(self,y1,x1,W,a,S,L):
        """
        y1 = height (m) of first channel in twin-jet
        x1 = value of smallest magnitude x-coordinate (assume width
        in x-direction
        W = width (m) width of each jet channel
        a = width (m) of the openings of both jets
        S = pitch (m) distance between centroids of twin jets (should be
        greater than a)
        L = length (m) length of the channel for twin jets prior to openings
        
        assumes the overall domain starts at z = 0.
        assumes overall domain is wider (in x-direction) than x1+W
        assumes overall domain is higher (in y-direction) than y1+S+a/2.
        """
    
        self.y1 = y1
        self.x1 = x1
        self.W = W
        self.a = a
        self.S = S
        self.L = L
    
    def get_Lo(self):
        """ 
           characteristic length is the hydraulic diameter. (flow area/wetted
           perimeter)
        """
        flow_area = self.a * self.W
        wetted_perimeter = 2.*self.a+2.*self.W
        return flow_area/wetted_perimeter
    
    def get_obstList(self,X,Y,Z):
        
       #x = np.array(X); 
        y = np.array(Y); 
        z = np.array(Z);
        obst_l = np.where(z < self.L)
        obst_h = np.where(z > 0.2)
        obst = np.intersect1d(obst_l[:],obst_h[:])
        y_dist1 = np.abs(y - (self.y1+self.a/2.))
        ch1 = np.where(y_dist1<self.a/2.)
        ch1 = np.intersect1d(obst[:],ch1[:])
        obst = np.setxor1d(obst[:],ch1[:])
        y_dist2 = np.abs(y - (self.y1+self.a/2.+self.S))
        ch2 = np.where(y_dist2<self.a/2.)
        ch2 = np.intersect1d(obst[:],ch2[:])
        obst = np.setxor1d(obst[:],ch2[:])
        return obst[:]

class GolfBall(EmptyChannel):
    """
     a channel with a golf ball obstacle
    """

    def __init__(self,SO,d_dimp,rd_dimp,N_e,N_a):
        """
           SO - pass in a sphericle obstacle as one of the arguments
           d_dimp = diameter of the dimples on the golf ball
           rd_dimp = radial distance of the center of the dimple from the center
                     of the golf ball
           N_e = number of dimples along all [0,pi] elevation angles 
           N_e = number of dimples along all [0,2pi] azimuthal angles
        """
        self.sphere = SO;
        self.d_dimp = d_dimp;
        self.rd_dimp = rd_dimp;
        self.N_e = N_e;
        self.N_a = N_a;

    def get_Lo(self):
        return self.sphere.get_Lo()


    def get_obstList(self,X,Y,Z):
        """
           return the obst list for the golf ball
        """
        obst_list1 = self.sphere.get_obstList(X,Y,Z)
        el_angles = np.linspace(0.,np.pi,self.N_e)
        
        x = np.array(X); y = np.array(Y); z = np.array(Z);
        print("removing the dimples")
        # start removing dimples
        iel = 0;
        for el in el_angles:
            iel+=1
        # for each elevation, we will get a different number of dimples
            N_az_el = np.floor(self.N_a*np.sin(el))+1;
            if N_az_el == 1:
                N_az_el+=1
            
            az_angles = np.linspace(0.,2.*np.pi, N_az_el, endpoint = False)
            print("removing dimples in elevation %g of %g" % (iel, len(el_angles)))
            iaz = 0;
            for az in az_angles:
              iaz+=1
              print("removing dimple %g of %g on this elevation" % (iaz,len(az_angles)))
              # get coordinates of the center of the spherical dimple
              y_c_d = self.sphere.y_c + self.rd_dimp*np.cos(el);
              z_c_d = self.sphere.z_c + self.rd_dimp*np.sin(az)*np.sin(el);
              x_c_d = self.sphere.x_c + self.rd_dimp*np.cos(az)*np.sin(el);
 
              dist = (x - x_c_d)**2 + (y - y_c_d)**2 + (z - z_c_d)**2
              dimples = np.where(dist <= ((self.d_dimp/2.))**2)
              obst_list1 = np.setxor1d(obst_list1[:],
                  np.intersect1d(obst_list1[:],dimples[:]))
             

        return obst_list1[:] 
        


class EllipticalScourPit(EmptyChannel):
    """
     a channel with an elliptical scour pit with prescribed properties
     corresponds to case 3 of Bryan's geometry_desc.m
    """

    def __init__(self,x_c,z_c,cyl_rad):
        """
          constructor giving the x and z coordinates of the scour pit along with
          the radius of the cylindrical piling
        """
        self.x_c = x_c
        self.z_c = z_c
        self.cyl_rad = cyl_rad

    def get_Lo(self):
        return self.cyl_rad*2.

    def get_obstList(self,X,Y,Z):
        """
         return a list of all indices of lattice points within the boundaries of the
         scour pit obstacle
        """
       
        ellip_a = 2.*2.*self.cyl_rad
        ellip_b = 2.*self.cyl_rad
        ellip_c = 8.*self.cyl_rad
        ellip_x = self.x_c
        ellip_z = self.z_c + self.cyl_rad
        ellip_y = ellip_b 

        floor_part = np.array(np.where(Y < ellip_b)).flatten()

        dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
        cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())

        scour_pit = np.array(np.where( (X - ellip_x)**2/(ellip_a**2) + 
                        (Y - ellip_y)**2/(ellip_b**2) +
                        (Z - ellip_z)**2/(ellip_c**2) <= 1.)).flatten()

        # remove the scour pit from the floor
        obst_list = np.setxor1d(floor_part[:], 
                        np.intersect1d(floor_part[:],scour_pit[:]))


        # then add the cylinder
        obst_list = np.union1d(obst_list[:],cyl_part[:])
        
        return list(obst_list[:])

class ConeScourPit(EmptyChannel):
    """
    a channel with a conical scour pit determined by the angle of repose of the soil particles (assumed to be river sand, phi=30 deg).  
    """

    def __init__(self,x_c,z_c,cyl_rad):
        """
          constructor giving the x and z coordinates of the scour pit along with the radius of the cylindrical piling
        """
        self.x_c = x_c
        self.z_c = z_c
        self.cyl_rad = cyl_rad

    def get_Lo(self):
        return self.cyl_rad*2.

    def get_obstList(self,X,Y,Z):
        """
         return a list of all indices of lattice points within the boundaries of the conical scour pit obstacle.  x_s is defined in 'Scour at marine structures' by Richard Whitehouse, 1998.  Assumes river sand with phi (angle of repose) equal to 30 degrees.  h_cone is equal to rad_cone*tan(30) = rad_cone*0.57735
        """
       
        x_c_cone = self.x_c
        z_c_cone = self.z_c
        y_c_cone = 0
        x_s = 2.25*2*self.cyl_rad
        rad_cone = x_s + self.cyl_rad
        h_cone = rad_cone*0.57735

        floor_part = np.array(np.where(Y < h_cone)).flatten()

        dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
        cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())

        scour_pit = np.array(np.where( (X - x_c_cone)**2 + (Z - z_c_cone)**2 <= ((self.cyl_rad/cone)/(h_cone))**2*(Y - y_c_cone)**2))

        # remove the scour pit from the floor
        obst_list = np.setxor1d(floor_part[:], 
                        np.intersect1d(floor_part[:],scour_pit[:]))


        # then add the cylinder
        obst_list = np.union1d(obst_list[:],cyl_part[:])
        
        return list(obst_list[:])

class SinglePile(EmptyChannel):
    """
    a channel with a single pile, no scour.  Used for comparison to both elliptical and conical scour pits.
    """

    def __init__(self,x_c,z_c,cyl_rad):
        """
          constructor giving the x and z coordinates of the piling center along with the radius of the cylindrical piling
        """
        self.x_c = x_c
        self.z_c = z_c
        self.cyl_rad = cyl_rad

    def get_Lo(self):
        return self.cyl_rad*2.

    def get_obstList(self,X,Y,Z):
        """
         return a list of all indices of lattice points within the boundaries of the bed Bed thickness is equal to the diameter of the piling (2x radius)
        """
       
        #Bed
        floor_part = np.array(np.where(Y < 2*self.cyl_rad)).flatten()

        #Piling
        dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
        cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())


        # then add the cylinder
        obst_list = np.union1d(floor_part[:],cyl_part[:])
        
        return list(obst_list[:])

class WavyBed(EmptyChannel):
    """
    a channel with a single pile, Sin-wave bottom.  
    """

    def __init__(self,x_c,z_c,cyl_rad):
        """
          constructor giving the x and z coordinates of the piling center along with the radius of the cylindrical piling
        """
        self.x_c = x_c
        self.z_c = z_c
        self.cyl_rad = cyl_rad

    def get_Lo(self):
        return self.cyl_rad*2.

    def get_obstList(self,X,Y,Z):
        """
        waveh and wavel are used to characterize the sine wave for the bed.  shallower sin waves do better in remaining stable throughout the simulation at low Reynolds numbers.

        """
        waveh = 0.125
        wavel = 5        
        floor_part = np.array(np.where(Y < (waveh*np.sin(wavel*Z) + 2*self.cyl_rad))).flatten()
        
        #Piling
        dist = (X - self.x_c)**2 + (Z - self.z_c)**2;
        cyl_part = list(np.array(np.where( dist < self.cyl_rad**2)).flatten())


        # then add the cylinder
        obst_list = np.union1d(floor_part[:],cyl_part[:])
        
        return list(obst_list[:])


class PipeContract(EmptyChannel):
    """
    a single smooth pipe with diameter in, diam_in, through a contraction and leaving at diameter out, diam_out.  Contraction assumed to be 45 degrees.  Channel assumed to be 2 x 2 x 8.  Lo = diam_out (smaller diameter).  Contraction begins at z = 4.  For a clean pipe, diam_in = 1.8 and diam_out = 0.8.
    """

    def __init__(self,diam_in,diam_out):
        """
 constructor identifying diameters into and out of contraction.  Recommend diam_in = 1.8 and diam_out = 0.8
        """
        self.diam_in = diam_in
        self.diam_out = diam_out

    def get_Lo(self):
        return self.diam_out

    def get_obstList(self,X,Y,Z):
        """
   Define areas external to pipe.
        """
       #Pipe in - find all points exterior of large pipe
        pipe_in = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_in/2)**2)).flatten()
        pipe_in_stop = np.array(np.where(Z <= 4)).flatten()
        pipe_in = np.intersect1d(pipe_in[:],pipe_in_stop[:])
    
        #Contraction - find all points exterior of contraction
        r_cone = self.diam_out
        h_cone = self.diam_out	
        contraction = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (r_cone/h_cone)**2*(Z - (4 + h_cone))**2)).flatten()
        contraction_start = np.array(np.where(Z >= 4)).flatten()
        contraction_stop = np.array(np.where(Z <= 4 + .5*self.diam_out)).flatten()
        contraction = np.intersect1d(contraction[:],contraction_start[:])
        contraction = np.intersect1d(contraction[:],contraction_stop[:])
    
        #Pipe out - final all points exterior of smaller pipe
        pipe_out = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_out/2)**2)).flatten()
        pipe_out_start = np.array(np.where(Z >= 4 + .5*self.diam_out)).flatten()
        pipe_out = np.intersect1d(pipe_out[:],pipe_out_start[:])
    
    
        #Put the pieces together
    
        #pipe = pipe_in[:]
        pipe = np.union1d(contraction[:],pipe_in[:])
        pipe = np.union1d(pipe[:],pipe_out[:])
    
        obst_list = pipe[:]

       
        return list(obst_list[:])

class PipeExpand(EmptyChannel):
    """
    opposite of pipe contraction.  a single smooth pipe with diameter in, diam_in, through an expansion and leaving at diameter out, diam_out.  Expansion assumed to be 45 degrees.  Channel assumed to be 2 x 2 x 8.  Lo = diam_in (smaller diameter).  Expansion begins at z = 4.  Best works when diam_in = 0.8 and diam_out = 1.8
    """

    def __init__(self,diam_in,diam_out):
        """
          constructor identifying pipe diameters into and out of expansion.  Recommend diam_in = 0.8 and diam_out = 1.8
        """
        self.diam_in = diam_in
        self.diam_out = diam_out

    def get_Lo(self):
        return self.diam_in

    def get_obstList(self,X,Y,Z):
        """
   Define areas external to pipe.
        """
       #Pipe in - find all points exterior of small
        pipe_in = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_in/2)**2)).flatten()
        pipe_in_stop = np.array(np.where(Z <= 1.5 + 0.5*(self.diam_out - self.diam_in))).flatten()
        pipe_in = np.intersect1d(pipe_in[:],pipe_in_stop[:])
    
        #Expansion - find all points exterior of expansion
        r_cone = self.diam_in
        h_cone = self.diam_in	
        expansion = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (r_cone/h_cone)**2*(Z - 1.5)**2)).flatten()
        expansion_start = np.array(np.where(Z >= 1.5 + 0.5*(self.diam_out - self.diam_in)))
        #expansion_stop = np.array(np.where(Z <= 4)).flatten()
        expansion = np.intersect1d(expansion[:],expansion_start[:])
        #expansion = np.intersect1d(expansion[:],expansion_stop[:])
    
        #Pipe out - final all points exterior of smaller pipe
        pipe_out = np.array(np.where((X - 1)**2 + (Y - 1)**2 > (self.diam_out/2)**2)).flatten()
        pipe_out_start = np.array(np.where(Z >= 1.5 + 0.5*(self.diam_in - self.diam_out))).flatten()
        pipe_out = np.intersect1d(pipe_out[:],pipe_out_start[:])
    
    
        #Put the pieces together
    
        pipe = expansion[:]
        pipe = np.union1d(expansion[:],pipe_in[:])
        pipe = np.union1d(pipe[:],pipe_out[:])

        obst_list = pipe[:]

       
        return list(obst_list[:])

class PipeTurn(EmptyChannel):
    """
  Provides an s-shaped pipe of constant radius with two 180-degree turns constructed out of constant-radius tori.  Diameter needs to be 0.5 for 
    """

    def __init__(self,diam_in):
        """
          constructor providing pipe diameter for use in Lo.  Use 0.5.
        """
        self.diam_in = diam_in

    def get_Lo(self):
        return self.diam_in

    def get_obstList(self,X,Y,Z):
        """
   Define areas external to pipe.
        """
       #Pipe_1
        pipe_1 = np.array(np.where((X - 1)**2 + (Y - 4)**2 >= 0.5**2)).flatten()
        pipe_1_stop_z = np.array(np.where(Z <= 3.0)).flatten()
        pipe_1_stop_y = np.array(np.where(Y >= 3.25)).flatten()
        pipe_1_stop = np.intersect1d(pipe_1_stop_z[:],pipe_1_stop_y[:])
        pipe_1 = np.intersect1d(pipe_1[:],pipe_1_stop[:])
    
        #Turn_1
        turn_1 = np.array(np.where((0.75 - np.sqrt((Y - 3.25)**2 + (Z -3)**2))**2 + (X - 1)**2 >= 0.5**2)).flatten()
        turn_1_stop_z = np.array(np.where(Z >= 3.0)).flatten()
        turn_1_stop_y = np.array(np.where(Y>= 1.75)).flatten()
        turn_1_stop = np.intersect1d(turn_1_stop_z[:],turn_1_stop_y[:])
        turn_1 = np.intersect1d(turn_1[:],turn_1_stop[:])
    
        #Pipe_2
        pipe_2 = np.array(np.where((X - 1)**2 + (Y - 2.5)**2 >= 0.5**2)).flatten()
        pipe_2_start_z = np.array(np.where(Z >= 1.5)).flatten()
        pipe_2_start_y_up = np.array(np.where(Y <= 3.25)).flatten()
        pipe_2_start_y_down = np.array(np.where(Y >= 1.75)).flatten()
        pipe_2_start_y = np.intersect1d(pipe_2_start_y_up[:],pipe_2_start_y_down[:])	
        pipe_2_start = np.intersect1d(pipe_2_start_z[:],pipe_2_start_y[:])
        pipe_2 = np.intersect1d(pipe_2[:],pipe_2_start[:])
        pipe_2_stop_z = np.array(np.where(Z <= 3.0)).flatten()
        pipe_2_stop_y = np.array(np.where(Y <= 3.25)).flatten()
        pipe_2_stop = np.intersect1d(pipe_2_stop_z[:],pipe_2_stop_y[:])
        pipe_2 = np.intersect1d(pipe_2[:],pipe_2_stop[:])
    
        #Turn_2
        turn_2 = np.array(np.where((0.75 - np.sqrt((Y - 1.75)**2 + (Z -1.5)**2))**2 + (X - 1)**2 >= 0.5**2)).flatten()
        turn_2_stop_z = np.array(np.where(Z <= 1.5)).flatten()
        turn_2_stop_y = np.array(np.where(Y <= 3.25)).flatten()
        turn_2_stop = np.intersect1d(turn_2_stop_z[:],turn_2_stop_y[:])
        turn_2 = np.intersect1d(turn_2[:],turn_2_stop[:])
        
        #Pipe_3
        pipe_3 = np.array(np.where((X - 1)**2 + (Y - 1.0)**2 >= 0.5**2)).flatten()
        pipe_3_start_z = np.array(np.where(Z >= 1.5)).flatten()
        pipe_3_start_y = np.array(np.where(Y <= 1.75)).flatten()
        pipe_3_start = np.intersect1d(pipe_3_start_z[:],pipe_3_start_y[:])
        pipe_3 = np.intersect1d(pipe_3[:],pipe_3_start[:])	
    
        #Put the pieces together
    
        pipe = np.union1d(pipe_1[:],turn_1[:])
        pipe = np.union1d(pipe[:],pipe_2[:])
        pipe = np.union1d(pipe[:],turn_2[:])	
        pipe = np.union1d(pipe[:],pipe_3[:])
    
        obst_list = pipe[:]

        return list(obst_list[:])

class PipeOut(EmptyChannel):
    """
  Class consisting of a single pipe of diam_in and length length_in exiting a wall into an open space.  
    """

    def __init__(self,diam_in,length_in):
        """
        defines the diameter and length (z axis) of pipe leading to open area
        """
        self.diam_in = diam_in
        self.length_in = length_in

    def get_Lo(self):
        return self.diam_in

    def get_obstList(self,X,Y,Z):
        """
            Define solid areas around pipe.  Everything else will be open.  Ensure coordinates for center of circle match center of Lx-Ly.
        """
        #Pipe In
        pipe_in = np.array(np.where((X - 0.5*(4))**2 + (Y - 0.5*(4))**2 >= (0.5*self.diam_in)**2)).flatten()
        pipe_in_stop = np.array(np.where(Z <= self.length_in)).flatten()
        pipe_in = np.intersect1d(pipe_in[:],pipe_in_stop[:])


        obst_list = pipe_in[:]

        return list(obst_list[:])

class Butterfly(EmptyChannel):
    """
  A geometry class that defines a fully open butterfly valve within a pipe of diam=1.0.  
    """

    def __init__(self,diam):
        """
          constructor identifying pipe diameter.  Must be a 1 diam pipe inside a 1.2 x 1.2 x 8 channel.  Valve center at z = 3.
        """
        self.diam = diam


    def get_Lo(self):
        return self.diam

    def get_obstList(self,X,Y,Z):
        """
   Define solid areas

        """

	#Pipe
        pipe = np.array(np.where((X - 0.6)**2 + (Y - 0.6)**2 >= 0.5**2)).flatten()

	#Seat
        seat = np.array(np.where((X - 0.6)**2 + (Y - 0.6)**2 >= 0.42**2)).flatten()
        seat_start = np.array(np.where(Z >= 2.975)).flatten()
        seat_stop = np.array(np.where(Z <= 3.025)).flatten()
        seat = np.intersect1d(seat[:],seat_start[:])
        seat = np.intersect1d(seat[:],seat_stop[:])

	#Pivot
        pivot = np.array(np.where((X - 0.6)**2 + (Z - 3)**2 <= 0.075**2)).flatten()

	#Front Disc
        front_disc = np.array(np.where((Y - 0.6)**2 + (Z - 3)**2 <= 0.5**2)).flatten()
        front_disc_stop = np.array(np.where(Z <= 3.0)).flatten()
        front_disc_x_min = np.array(np.where(X >= 0.525)).flatten()
        front_disc_x_max = np.array(np.where(X <= 0.575)).flatten()

        front_disc = np.intersect1d(front_disc[:],front_disc_stop[:])
        front_disc = np.intersect1d(front_disc[:],front_disc_x_min[:])
        front_disc = np.intersect1d(front_disc[:],front_disc_x_max[:])

	#Back Disc
        back_disc = np.array(np.where((Y - 0.6)**2 + (Z - 3)**2 <= 0.5**2)).flatten()
        back_disc_start = np.array(np.where(Z >= 3.0)).flatten()
        back_disc_x_min = np.array(np.where(X >= 0.625)).flatten()
        back_disc_x_max = np.array(np.where(X <= 0.675)).flatten()

        back_disc = np.intersect1d(back_disc[:],back_disc_start[:])
        back_disc = np.intersect1d(back_disc[:],back_disc_x_min[:])
        back_disc = np.intersect1d(back_disc[:],back_disc_x_max[:])

	#Put the pieces together

        valve = np.union1d(pipe[:],seat[:])
        valve = np.union1d(valve[:],pivot[:])
        valve = np.union1d(valve[:],front_disc[:])
        valve = np.union1d(valve[:],back_disc[:])
	
        obst_list = valve[:]

        return list(obst_list[:])


class Tee(EmptyChannel):
    """
  establishes a single large pipe with a "tee" into a smaller pipe that loops up and around before rejoining the main line.  The Main line undergoes a contraction after the tee but before the rejoining secondary line.  diam_2 should be smaller than diam_1.
    """

    def __init__(self,diam_1,diam_2):
        """
        Constructor identifying the diameters of the two pipes.  Pipe 1 runs straight through from Z_min to Z_max.  Pipe 2 tees off and runs parallel to Pipe 1.  Pipe 1 enters/exits z planes at y = 1.  Pipe 2 runs at y = 3.  Assumes dimensions of space (X,Y,Z) is (2,4,8).
        """
        self.diam_1 = diam_1
        self.diam_2 = diam_2

    def get_Lo(self):
        return self.diam_1

    def get_obstList(self,X,Y,Z):
        """
   Define solid areas

        """

	#Pipe 1
        pipe_1a = np.array(np.where((X - 1)**2 + (Y - 1)**2 <= (self.diam_1/2)**2)).flatten()
        pipe_1a_stop = np.array(np.where(Z<=4.)).flatten()
        pipe_1a = np.intersect1d(pipe_1a[:],pipe_1a_stop[:])	
        pipe_1b = np.array(np.where((X - 1)**2 + (Y - 1)**2 <= (self.diam_1/4)**2)).flatten()
        pipe_1 = np.union1d(pipe_1a[:],pipe_1b[:])

	#Pipe 2 Tee Off
        tee_1 = np.array(np.where((X - 1)**2 + (Z - 1.5)**2 <= (self.diam_2/2)**2)).flatten()
        tee_1_start = np.array(np.where(Y >= 1)).flatten()
        tee_1_end = np.array(np.where(Y <= 3 - 0.5*self.diam_2)).flatten()
        tee_1 = np.intersect1d(tee_1[:],tee_1_start[:])
        tee_1 = np.intersect1d(tee_1[:],tee_1_end[:])

	#Pipe 2 Elbow 1
        elbow_1 = np.array(np.where((self.diam_2/2 - np.sqrt((Y - (3 - self.diam_2/2))**2 + (Z -(1.5 + self.diam_2/2))**2))**2 + (X - 1)**2 <= (self.diam_2/2)**2)).flatten()
        elbow_1_start = np.array(np.where(Y >= 3- 0.5*self.diam_2)).flatten()
        elbow_1_stop = np.array(np.where(Z <= 1.5 + self.diam_2/2)).flatten()
        elbow_1 = np.intersect1d(elbow_1[:],elbow_1_start[:])
        elbow_1 = np.intersect1d(elbow_1[:],elbow_1_stop[:])


	#Pipe 2
        pipe_2 = np.array(np.where((X - 1)**2 + (Y - 3)**2 <= (self.diam_2/2)**2)).flatten()
        pipe_2_start = np.array(np.where(Z >= 1.5 + self.diam_2/2)).flatten()
        pipe_2_stop = np.array(np.where(Z <= 5 - self.diam_2/2)).flatten()
        pipe_2 = np.intersect1d(pipe_2[:],pipe_2_start[:])
        pipe_2 = np.intersect1d(pipe_2[:],pipe_2_stop[:])

	#Pipe 2 Elbow 2
        elbow_2 = np.array(np.where((self.diam_2/2 - np.sqrt((Y - (3 - self.diam_2/2))**2 + (Z -(5- self.diam_2/2))**2))**2 + (X - 1)**2 <= (self.diam_2/2)**2)).flatten()
        elbow_2_start = np.array(np.where(Y >= 3- 0.5*self.diam_2)).flatten()
        elbow_2_stop = np.array(np.where(Z >= 5- self.diam_2/2)).flatten()
        elbow_2 = np.intersect1d(elbow_2[:],elbow_2_start[:])
        elbow_2 = np.intersect1d(elbow_2[:],elbow_2_stop[:])

	#Pipe 2 Tee In
        tee_2 = np.array(np.where((X - 1)**2 + (Z - 5)**2 <= (self.diam_2/2)**2)).flatten()
        tee_2_start = np.array(np.where(Y >= 1)).flatten()
        tee_2_end = np.array(np.where(Y <= 3 - 0.5*self.diam_2)).flatten()
        tee_2 = np.intersect1d(tee_2[:],tee_2_start[:])
        tee_2 = np.intersect1d(tee_2[:],tee_2_end[:])

        empty = np.array(np.where(Y>=0.)).flatten() 

	#Put the pieces together
        pipe = np.union1d(pipe_1[:],tee_1[:])
        pipe = np.union1d(pipe[:],elbow_1[:])
        pipe = np.union1d(pipe[:],pipe_2[:])
        pipe = np.union1d(pipe[:],elbow_2[:])
        pipe = np.union1d(pipe[:],tee_2[:])
        pipe = np.setxor1d(pipe[:], empty[:])	
        obst_list = pipe[:]

        return list(obst_list[:])

def fluid_properties(fluid_str):  
   """
   Return the physical density and kinematic viscosity for the prescribed
   fluid.
   
   """
   fluid_lib = {'water':(1000., 1.0e-6), 
                'glycol':(965.3,0.06/965.3),
                'glycerin':(1260.0,1.49/1260.0)}
   if fluid_str in list(fluid_lib.keys()):
     return fluid_lib[fluid_str]
   else:
     print('valid fluids are:')
     for keys in fluid_lib:
       print(" '%s' " % keys)
     raise KeyError('invalid fluid specified')

class LidDrivenCavity:
    def __init__(self,Lx_p = 1.,
                      Ly_p = 1., Lz_p = 1.,
                      fluid='water',
                      N_divs=11) :
      """
      class constructor - lid driven cavity.
      All surfaces xm,xp,ym,zm,zp are solid except
      the "lid" (yp) which is a moving surface
      (node type 5)
      """
      self.Lx_p = Lx_p
      self.Ly_p = Ly_p
      self.Lz_p = Lz_p
      self.N_divs = N_divs
      self.fluid = fluid
      
      # generate the geometry
      Lo = Ly_p # by convention (for now)
      self.Lo = Lo
      self.Ny = math.ceil((N_divs-1)*(Ly_p/Lo))+1
      self.Nx = math.ceil((N_divs-1)*(Lx_p/Lo))+1
      self.Nz = math.ceil((N_divs-1)*(Lz_p/Lo))+1
      self.nnodes = self.Nx*self.Ny*self.Nz
      print("Creating channel with %g lattice points." % self.nnodes)
      x = np.linspace(0.,Lx_p,self.Nx).astype(np.float32);
      y = np.linspace(0.,Ly_p,self.Ny).astype(np.float32);
      z = np.linspace(0.,Lz_p,self.Nz).astype(np.float32);
   
      Y,Z,X = np.meshgrid(y,z,x);
    
      self.x = np.reshape(X,int(self.nnodes))
      self.y = np.reshape(Y,int(self.nnodes))
      self.z = np.reshape(Z,int(self.nnodes))
      self.dx = x[1]-x[0]
      # get fluid properties from the included fluid library
      self.rho_p, self.nu_p = fluid_properties(fluid)
      self.set_cavity_walls()
      self.ndType = np.zeros((self.nnodes,),dtype=np.int32)
      self.ndType[self.solid_list]=1
      self.ndType[self.lid_list]=5
      
       # must have geometry set first
    def set_cavity_walls(self,walls=['left','right','bottom','west','east']): 
        """
         set up to 5 walls as solid walls for the simulation
        """
        solid_list_a = np.empty(0).flatten()
        solid_list_b = np.empty(0).flatten()
        solid_list_c = np.empty(0).flatten()
        solid_list_d = np.empty(0).flatten()
        solid_list_e = np.empty(0).flatten()

        for w in walls:
            if w=='right':
                solid_list_a = np.array(np.where((self.x==0.))).flatten()
            elif w=='left':
                solid_list_b = np.array(np.where((self.x > (self.Lx_p-self.dx/2.)))).flatten()
            elif w=='west':
                solid_list_d = np.array(np.where((self.z == 0.))).flatten()
            elif w=='bottom':
                solid_list_c = np.array(np.where((self.y == 0.))).flatten()
            elif w=='east':
                solid_list_e = np.array(np.where((self.z > (self.Lz_p - self.dx/2.)))).flatten()

        solid_list = np.array(np.union1d(solid_list_a,solid_list_b)); 
        solid_list = np.array(np.union1d(solid_list,solid_list_c));
        solid_list = np.array(np.union1d(solid_list,solid_list_e));
        self.solid_list = np.array(np.union1d(solid_list,solid_list_d))
        
        self.lid_list = np.array(np.where((self.y > (self.Ly_p-self.dx/2.)))).flatten()
#        self.lid_list = np.setxor1d(self.lid_list[:],
#            np.intersect1d(self.lid_list[:],self.solid_list[:]))
        
        
    def write_mat_file(self, geom_filename):
        """
          generate the mat file to interface with genInput.py.  Needs to save
          Lx_p, Ly_p, Lz_p, Lo, Ny_divs, rho_p, nu_p, and ndType.
          note that the snl and obst_list need to be combined into one list 
        """
        mat_dict = {}
        mat_dict['Lx_p'] = self.Lx_p
        mat_dict['Ly_p'] = self.Ly_p
        mat_dict['Lz_p'] = self.Lz_p
        mat_dict['Lo'] = self.Lo
        mat_dict['Ny_divs'] = self.N_divs
        mat_dict['rho_p'] = self.rho_p
        mat_dict['nu_p'] = self.nu_p

        mat_dict['ndType'] = list(self.ndType[:])

        scipy.io.savemat(geom_filename,mat_dict)
        
    def write_bc_vtk(self):
        """
         write node lists to properly formatted VTK files
        """
        print("Creating boundary condition arrays")
        # perhaps at some other time this class should be extended
        # to allow obstacles within the lid-driven cavity
#        obst_array = np.zeros(self.nnodes)
#        obst_array[list(self.obst_list)] = 100.
#
#        #print type(self.inlet_list)
#        inlet_array = np.zeros(self.nnodes)
#        inlet_array[list(self.inlet_list)] = 200.
#
#        outlet_array = np.zeros(self.nnodes)
#        outlet_array[list(self.outlet_list)] = 300.
        print("length of lid_list = %d \n"%(len(list(self.lid_list))))
        lid_array = np.zeros(int(self.nnodes))
        lid_array[list(self.lid_list)] = 700.

        solid_array = np.zeros(int(self.nnodes))
        solid_array[list(self.solid_list)] = 500.
        
        dims = [int(self.Nx), int(self.Ny), int(self.Nz)]
        origin = [0., 0., 0.]
        dx = self.x[1] - self.x[0]
        spacing = [dx, dx, dx] #uniform lattice
        
        print("Writing boundary conditions to VTK files")
        #writeVTK(inlet_array,'inlet','inlet.vtk',dims,origin,spacing)
        #writeVTK(outlet_array,'outlet','outlet.vtk',dims,origin,spacing)
        #writeVTK(obst_array,'obst','obst.vtk',dims,origin,spacing)
        writeVTK(lid_array,'lid','lid.vtk',dims,origin,spacing)
        writeVTK(solid_array,'solid','solid.vtk',dims,origin,spacing)
                         

class FluidChannel:
    def __init__(self,Lx_p=1.,
        Ly_p=1.,
        Lz_p=6.,
        fluid='water', 
        obst=EmptyChannel(1.),
        N_divs = 5,
        wallList=['left','right','top','bottom']):
        """
         class constructor
        """
        self.Lx_p = Lx_p
        self.Ly_p = Ly_p
        self.Lz_p = Lz_p
        self.N_divs = N_divs
        self.fluid = fluid
        self.obst = obst
        
        self.subsetList = []; # make a list of subsets visible in the class interface.
        

        # generate the geometry

        Lo = obst.get_Lo()

        self.Ny = math.ceil((N_divs-1)*(Ly_p/Lo))+1
        self.Nx = math.ceil((N_divs-1)*(Lx_p/Lo))+1
        self.Nz = math.ceil((N_divs-1)*(Lz_p/Lo))+1
        self.nnodes = self.Nx*self.Ny*self.Nz
        print("Creating channel with %g lattice points." % self.nnodes)
        x = np.linspace(0.,Lx_p,self.Nx).astype(np.float32);
        y = np.linspace(0.,Ly_p,self.Ny).astype(np.float32);
        z = np.linspace(0.,Lz_p,self.Nz).astype(np.float32);
   
        Y,Z,X = np.meshgrid(y,z,x);
    
        self.x = np.reshape(X,int(self.nnodes))
        self.y = np.reshape(Y,int(self.nnodes))
        self.z = np.reshape(Z,int(self.nnodes))
        
        self.dx = self.x[1]-self.x[0]
        
        self.set_pRef_indx(Lx_p/2.,Ly_p/2.,0.95*Lz_p);

        # get fluid properties from the included fluid library
        self.rho_p, self.nu_p = fluid_properties(fluid)

        # identify inlet and outlet nodes - 
        # require the user to set solid boundaries separately
        self.inlet_list = np.where(self.z==0)
        self.outlet_list = np.where(self.z==Lz_p)
        
        print("Getting obstacle list")
        # get obstacle list
        self.obst_list = self.obst.get_obstList(self.x[:],self.y[:],self.z[:])
        

        print("Generating channel solid boundaries")
        # set channel walls
        self.set_channel_walls(wallList)

        # now eliminate overlap between node lists

        self.inlet_list = np.setxor1d(self.inlet_list[:],
            np.intersect1d(self.inlet_list[:],self.solid_list[:]))
        self.inlet_list = np.setxor1d(self.inlet_list[:],
            np.intersect1d(self.inlet_list[:],self.obst_list[:]))
        
        self.outlet_list = np.setxor1d(self.outlet_list[:],
            np.intersect1d(self.outlet_list[:],self.solid_list[:]))
        self.outlet_list = np.setxor1d(self.outlet_list[:],
            np.intersect1d(self.outlet_list[:],self.obst_list[:]))

        self.obst_list = np.setxor1d(self.obst_list[:],
            np.intersect1d(self.obst_list[:],self.solid_list[:]))
            
        
        self.ndType = np.zeros((int(self.nnodes),),dtype=np.int32)
        # node types:
        # regular fluid node: 0
        # solid node: 1
        # velocity zm (inlet) node: 2
        # pressure zp (outlet) node: 3
        self.ndType[np.union1d(self.obst_list,self.solid_list).astype(np.int32)] = 1
        self.inlet_list = np.array(self.inlet_list).astype(np.int32);
        self.outlet_list = np.array(self.outlet_list).astype(np.int32);
        self.ndType[self.inlet_list] = 2
        self.ndType[self.outlet_list] = 3
       
    def write_mat_file(self, geom_filename):
        """
          generate the mat file to interface with genInput.py.  Needs to save
          Lx_p, Ly_p, Lz_p, Lo, Ny_divs, rho_p, nu_p, and ndType.
          note that the snl and obst_list need to be combined into one list 
        """
        mat_dict = {}
        mat_dict['Lx_p'] = self.Lx_p
        mat_dict['Ly_p'] = self.Ly_p
        mat_dict['Lz_p'] = self.Lz_p
        mat_dict['Lo'] = self.obst.get_Lo()
        mat_dict['Ny_divs'] = self.N_divs
        mat_dict['rho_p'] = self.rho_p
        mat_dict['nu_p'] = self.nu_p

        mat_dict['ndType'] = list(self.ndType[:])
        mat_dict['pRef_idx'] = self.pRef_indx;
        
        
        # add a list of subset nodes
        
        mat_dict['ssNds']=self.get_subset_nodeSet()
        
        if geom_filename[-4:] != ".mat":
            geom_filename+=".mat"

        scipy.io.savemat(geom_filename,mat_dict,appendmat=False)

    def add_subset(self,ss):
        """
        add a subset object to the list of subsets associated with this fluid channel
        
        """
        
        self.subsetList.append(ss)
        
    def get_subset_nodeSet(self):
        """
        convert list of subset objects to set of fluid channel nodes
        contained within the subsets.
        """
        nodeSet = np.empty(0)
        for ss in self.subsetList:
            nodeSet = np.union1d(nodeSet,ss.get_members(self.x[:],self.y[:],self.z[:]))
            
        return list(nodeSet.flatten())
    def set_pRef_indx(self,Xref,Yref,Zref):
        """
        find the node index within the fluid channel that can
        serve as the pressure reference.
        
        By default this pressure reference will be placed at 95% of 
        the channel length (in the Z-direction) at the X-and Y- center 
        of the channel.  Save this pRef index in the geometry file and
        write it into the params file for post-processing"
        """
        xInd = np.where(np.abs(self.x - Xref)<self.dx/1.5);
        yInd = np.where(np.abs(self.y - Yref)<self.dx/1.5);
        zInd = np.where(np.abs(self.z - Zref)<self.dx/1.5);
        
        ndList = np.intersect1d(xInd,yInd);
        ndList = np.intersect1d(zInd,ndList)
        
        if len(ndList>0):
            self.pRef_indx = ndList[0]
        else:
            self.pRef_indx = 0
            print("Warning: failed to select pRef node!")
            
        
        
    
    def write_bc_vtk(self):
        """
         write node lists to properly formatted VTK files
        """
        print("Creating boundary condition arrays")
        obst_array = np.zeros(int(self.nnodes))
        obst_array[list(self.obst_list)] = 100.

        #print type(self.inlet_list)
        inlet_array = np.zeros(int(self.nnodes))
        inlet_array[list(self.inlet_list)] = 200.

        outlet_array = np.zeros(int(self.nnodes))
        outlet_array[list(self.outlet_list)] = 300.

        solid_array = np.zeros(int(self.nnodes))
        solid_array[list(self.solid_list)] = 500.
        
        dims = [int(self.Nx), int(self.Ny), int(self.Nz)]
        origin = [0., 0., 0.]
        dx = self.x[1] - self.x[0]
        spacing = [dx, dx, dx] #uniform lattice
        
        print("Writing boundary conditions to VTK files")
        writeVTK(inlet_array,'inlet','inlet.vtk',dims,origin,spacing)
        writeVTK(outlet_array,'outlet','outlet.vtk',dims,origin,spacing)
        writeVTK(obst_array,'obst','obst.vtk',dims,origin,spacing)
        writeVTK(solid_array,'solid','solid.vtk',dims,origin,spacing)
        
        if len(self.subsetList) > 0:
            print("Writing subset lists to VTK files")
            
            ss_idx = 0;
            for ss in self.subsetList:
                print('Writing subset %d to file'%ss_idx)
                ssNodes = ss.get_members(self.x[:],self.y[:],self.z[:])
                ss_array = np.zeros(int(self.nnodes))
                ss_array[list(ssNodes)] = ss_idx*1000. + 1000.
                ssName = 'subSpace'+str(ss_idx)
                ssFile = ssName+'.vtk'
                writeVTK(ss_array,ssName,ssFile,dims,origin,spacing)
                ss_idx+=1
            


     # must have geometry set first
    def set_channel_walls(self,walls=['left','right','top','bottom']): 
        """
         set up to 4 walls as solid walls for the simulation
        """
        solid_list_a = np.empty(0).flatten()
        solid_list_b = np.empty(0).flatten()
        solid_list_c = np.empty(0).flatten()
        solid_list_d = np.empty(0).flatten()

        for w in walls:
            if w=='right':
                solid_list_a = np.array(np.where((self.x==0.))).flatten()
            elif w=='left':
                solid_list_b = np.array(np.where((self.x > (self.Lx_p-self.dx/2.)))).flatten()
            elif w=='top':
                solid_list_d = np.array(np.where((self.y > (self.Ly_p-self.dx/2.)))).flatten()
            elif w=='bottom':
                solid_list_c = np.array(np.where((self.y == 0.))).flatten()

        solid_list = np.array(np.union1d(solid_list_a,solid_list_b)); 
        solid_list = np.array(np.union1d(solid_list,solid_list_c))
        self.solid_list = np.array(np.union1d(solid_list,solid_list_d))

