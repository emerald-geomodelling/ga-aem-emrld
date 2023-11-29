#!/bin/tcsh

#Standalone
#galeisbstdem.exe galeisbstdem.con

#Use 4 MPI processes
#mpirun -np 4 galeisbstdem.exe galeisbstdem.con

#Use 3 OpenMP threads
galeisbstdem.exe galeisbstdem.con 3

