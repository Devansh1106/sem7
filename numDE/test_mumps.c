#include "mpi.h"
#include "dmumps_c.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    DMUMPS_STRUC_C id;
    memset(&id, 0, sizeof(id));
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);

    id.job = -1;  // init
    dmumps_c(&id);

    int n = 5;
    int nnz = 7;

    int irn[7] = {1,2,3,4,5,2,3};
    int jcn[7] = {1,2,3,4,5,1,2};
    double a[7] = {10,10,10,10,10,1,1};
    double rhs[5] = {1,1,1,1,1};

    if (myid == 0) {
        id.n = n;
        id.nz = nnz;
        id.irn = irn;
        id.jcn = jcn;
        id.a = a;
        id.rhs = rhs;
    }

    id.job = 6;
    dmumps_c(&id);

    id.job = -2;
    dmumps_c(&id);

    MPI_Finalize();
    return 0;
}
