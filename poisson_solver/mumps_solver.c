#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "dmumps_c.h"

int main(int argc, char **argv)
{
    MUMPS_INT n;
    MUMPS_INT8 len_a;
    DMUMPS_STRUC_C id;
    int* irn;
	int* jcn;
	double* rhs;
	double* a;
    int size, error = 0, n1;   // for MPI
    MUMPS_INT myid, ierr;
    double start_time = 0.0, end_time = 0.0;
    #define JOB_INIT -1
    #define JOB_END -2
    #define USE_COMM_WORLD -987654

    ierr = MPI_Init(&argc, &argv);		//MPI initialization
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0){
        FILE *fp = fopen("../numDE/cp.txt", "r");
        if (!fp) { perror("../numDE/cp.txt"); MPI_Abort(MPI_COMM_WORLD, 1); }

        // 1. Read n and nnz
        fscanf(fp, "%d %d", &n1, &len_a);
        n = (n1*n1);
        
        // ------------------------ TEST PROBLEM -------------------------------
        // irn = malloc(2 * sizeof(int));
        // jcn = malloc(2 * sizeof(int));
        // a = malloc(2 * sizeof(double));
        // rhs = malloc(2 * sizeof(double));
        // irn[0] = 1; irn[1] = 2;
        // jcn[0] = 1; jcn[1] = 2;
        // a[0] = 1; a[1] = 2;
        // rhs[0] = 1; rhs[1] = 4;
        // n = 2; len_a = 2; n1 = 2;

        // ------------------------- TEST PROBLEM ENDS -------------------------

        // 2. Allocate
        irn = malloc(len_a * sizeof(int));
        jcn = malloc(len_a * sizeof(int));
        a = malloc(len_a * sizeof(double));
        rhs = malloc(n * sizeof(double)); // n is already n1*n1

        if (!irn || !jcn || !a || !rhs) {
            fprintf(stderr, "Host Error: Memory allocation failed.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int k = 0;
        while (k < len_a && fscanf(fp, "%d %d %lf", &irn[k], &jcn[k], &a[k]) == 3){
            k++;
        }
        printf("n, len_a: %d %d\n", n, len_a);
        int i = 0;
        while (i < n && fscanf(fp, "%lf\n", &rhs[i]) == 1){
            i++;
        }
            
        fclose(fp);
        printf("File reading done!");
    } // myid == 0

    //record the start time
	if (myid == 0){
        start_time = MPI_Wtime();
    }

    id.comm_fortran = USE_COMM_WORLD;
	id.par = 1; id.sym = 0;
	id.job = JOB_INIT;
    dmumps_c(&id);
  
    
    // define the problem on host
	if(myid == 0)
	{
        id.n = n; id.nnz = len_a;
        id.irn = irn; id.jcn = jcn;
		id.a = a; id.rhs = rhs;
	}
	#define ICNTL(I) icntl[(I)-1] 
	id.ICNTL(1) = 6; id.ICNTL(2) = 6; id.ICNTL(3) = 6; id.ICNTL(4) = 0;//Supressing error msgs
    id.ICNTL(18) = 0; id.ICNTL(5) = 0;

    
	//Call the MUMPS package (analyze , factorization and solve)
    id.job = 6;             // combined 1,2,3 jobs
    dmumps_c(&id);

	if (id.infog[0] < 0)
	{
		printf("(PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
		myid, id.infog[0], id.infog[1]);
		error = 1;
	}

	// Terminate instance
	id.job = JOB_END;
	dmumps_c(&id);

    if (myid  == 0)
    {   
        end_time = MPI_Wtime();
        printf("Time Taken = %lf\n", (end_time - start_time));
        FILE *fp1 = fopen("sol.txt", "w");
        if (!fp1) { perror("sol.txt"); exit(1); }
    
        printf("sol:\n");
        for (int j = 0; j < n1; j++) {
            for (int i = 0; i < n1; i++) {
                fprintf(fp1, "%.17e ", rhs[i * n1 + j]);
            }
            fprintf(fp1, "\n");
        }

        fclose(fp1);
        printf("Solution is in sol.txt file!");
        free(irn);
        free(jcn);
        free(a);
        free(rhs);
    }
    ierr = MPI_Finalize();
    return 0;

}