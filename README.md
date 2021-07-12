# tLBM
Not "The Lego Batman Movie"; just another LBM repository.

Like so many repositories, this one is under development.  When things are running it will be a 3D lattice Boltzmann based computational fluid dynamics solver.  The problems that can be solved: simulate viscous fluid flow through a parallelepiped channel with a specified velocity inlet and specified pressure outlet.  Some flexibility is included with the upper/lower/right/left walls of the box to allow for logical periodicity. (e.g. simulate flow over an infinite array of cylinders)  

Geometry definition, domain partitioning, and general input data management is done with Python 3. 

Output data file processing and packaging (currently, into HDF5-formatted) binary data files for visualization (currently, with ParaView) will also be done with Python 3.

The main computational pieces will be written in C++ with MPI libraries for distributed parallelism.  Branches enabled with CUDA and/or OpenCL will (over time) be developed and made available when ready.

This is a research code aimed at gaining greater insight into fluid flows through CFD, but also a platform for which students can experiment with distributed parallel HPC.  Please keep your expectations low regarding results.

Requirements:
1. CMake
2. MPI
3. HDF5

Python Requirements for pre- and post-processing:
1. numpy
2. scipy
3. h5py
4. pymetis
5. pyvista

