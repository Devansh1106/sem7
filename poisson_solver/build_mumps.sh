#!/bin/bash
MUMPS_DIR=/home/devansh/MUMPS_5.5.1
mpicc mumps_solver.c -I$MUMPS_DIR/include -L$MUMPS_DIR/lib -ldmumps -lmumps_common -lpord -lblas -o mumps_solver
