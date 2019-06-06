# tLBM
Not "The Lego Batman Movie"; just another LBM repository.

Like so many repositories, this one is under development.  When things are running it will be a 3D lattice Boltzmann based computational fluid dynamics solver.  The problems that can be solved: simulate viscouse fluid flow through a parallelepiped channel with a specified velocity inlet and specified pressure outlet.  Some flexibility is included with the upper/lower/right/left walls of the box to allow for logical periodicity. (e.g. simulate flow over an infinite array of cylinders)  

Geometry definition, domain partitioning, and general input data management is done with Python.  

Output data file processing and pacakging (currently, into HDF5-formatted) binary data files for visualization (currently, with ParaView) will also be done with Python.

The main computational pieces will be written in C++ with MPI libraries for distributed parallelism.  Branches enabled with CUDA and/or OpenCL will (over time) be developed and made available when ready.

This is a research code aimed at gaining greater insight into fluid flows through CFD, but also a platform for which students can experiment with distributed parallel HPC.  Please keep your expectations low regarding results.
