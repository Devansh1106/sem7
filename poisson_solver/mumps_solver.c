#include <stdio.h>
#include <stdlib.h>
// #include <string.h>
#include "mpi.h"
#include "dmumps_c.h"

int main(int argc, char **argv)
{
    MUMPS_INT n;
    MUMPS_INT8 len_a;
    DMUMPS_STRUC_C id;
    // memset(&id, 0, sizeof(id));
    MUMPS_INT* irn;
	MUMPS_INT* jcn;
	double* rhs;
	double* a;
    int size, error = 0;   // for MPI
    MUMPS_INT myid, ierr;
    int _n, _len_a;
    double start_time = 0.0, end_time = 0.0;
    #define JOB_INIT -1
    #define JOB_END -2
    #define USE_COMM_WORLD -987654

    ierr = MPI_Init(&argc, &argv);		//MPI initialization
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    // printf("%d\n", myid);
    // ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
    // id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);

    // if (myid == 0){                                      // myid == 0 should read the file
                                                            // and distribute it to every rank
        FILE *fp = fopen("matrix.txt", "r");
        if (!fp) { perror("matrix.txt"); exit(1); }

        // 1. Read n and nnz
        fscanf(fp, "%d %d", &_n, &_len_a);
        // printf("%d %d ", n ,len_a);

        // 2. Allocate
        irn = malloc(_len_a * sizeof(int));
        jcn = malloc(_len_a * sizeof(int));
        a = malloc(_len_a * sizeof(double));
        rhs = malloc(_n *_n * sizeof(double));
        // printf(sizeof(rhs));

        // 3. Read matrix entries
        for (int k = 0; k < _len_a; k++){
            fscanf(fp, "%d %d %lf", &irn[k], &jcn[k], &a[k]);
            // printf("%d\n", irn[k]);
            // printf("%d\n", jcn[k]);
            // printf("%lf\n", a[k]);
        }
            // printf("%d\n", rhs[k]);

        // 4. Read RHS
        // for (int i = 0; i != EOF; i++){
        //     fscanf(fp, "%lf", &rhs[i]);
        //     printf("%lf\n", rhs[i]);
        // }
        int i = 0;
        while (fscanf(fp, "%lf ", &rhs[i]) == 1){
            // printf("%lf\n", rhs[i]);
            i++;
        }

        // while (fscanf(fp, " %lf,", &val) == 1) { 
        //     rhs

        fclose(fp);
        if (myid == 0) {printf("File reading done!");}
        // printf("%d\n", irn[1]);
        // printf("%d\n", jcn[1]);
        // printf("%d\n", a[1]);
        // printf("%lf\n", rhs[3]);

    // }
    // else{
    //     irn = NULL;
    //     jcn = NULL;
    //     a = NULL;
    //     rhs = NULL;
    // }

	if (myid == 0){
        // for (int i = 0; i < 16; i++){
            // printf("%lf\n", rhs[i]);
        // }
        start_time = MPI_Wtime();
    }		//record the start time

    id.comm_fortran = USE_COMM_WORLD;
	id.par = 1; id.sym = 0;
	id.job = JOB_INIT;
    dmumps_c(&id);
    
    n = _n * _n;
    len_a = _len_a;
    // printf("%d %d\n", n, len_a);
	// if (myid == 0){
        //     for (int i = 0; i < n; i++){
            //         printf("%lf\n", rhs[i]);
            //     }
            //     // start_time = MPI_Wtime();
            // }	
    
    
    // define the problem on host
	if(myid == 0)
	{
        // // printf("%lf\n", irn[7]);
        // for (int i = 0; i < n; i++){
        //         printf("%lf\n", rhs[i]);
        //     }
        id.n = n; id.nnz = len_a; id.irn=irn; id.jcn = jcn;         // in 2D system size is n x n !!!
		id.a = a; id.rhs = rhs;
        // if (myid == 0){ printf("%lf\n", id.rhs[4]);}
        // printf("%lf\n", id.rhs[4]);
	}
	#define ICNTL(I) icntl[(I)-1] 
    // id.ICNTL(21) = 0;
	id.ICNTL(1) = 6; id.ICNTL(2) = 6; id.ICNTL(3) = 6; id.ICNTL(4) = 0;//Supressing error msgs

	// if (myid == 0){
    //     for (int i = 0; i < n; i++){
    //         printf("%d\n", i);
    //         printf("%lf\n", id.rhs[i]);
    //     }
    //     // start_time = MPI_Wtime();
    // }	
	//Call the MUMPS package (analyze , factorization and solve)
    id.job = 6;
    dmumps_c(&id);
	// id.job = 1;         // combine job=1,2,3
	// dmumps_c(&id);
    // id.job = 2;         // combine job=1,2,3
	// dmumps_c(&id);
    // id.job = 3;         // combine job=1,2,3
	// dmumps_c(&id);

    // if (myid == 0){ printf("%lf\n", id.rhs[4]);}


    // error printing
	if (id.infog[0] < 0)
	{
		printf("(PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
		myid, id.infog[0], id.infog[1]);
		error = 1;
	}
    // else{
    //     if (myid == 0){printf("%lf\n", rhs[1]);}
    // }

	// Terminate instance
	id.job = JOB_END;
	dmumps_c(&id);

    // if (myid == 0)
    // {
    //     if (!error){
            // stop_time = MPI_Wtime();
            // return rhs, start_time - stop_time;
    //     }
    //     else{
    //         printf("Error!, not able to return rhs.")
    //         stop_time = MPI_Wtime();
    //         return -1, start_time - stop_time;
    //     }
    // }
    if (myid  == 0)
    {
        end_time = MPI_Wtime();
        printf("Time Taken = %lf\n", (end_time - start_time));
        // printf("%lf", &rhs[1]);
        FILE *fp1 = fopen("sol.txt", "w");
        if (!fp1) { perror("sol.txt"); exit(1); }
    
        // // 4. Write sol
        // for (int i = 0; i < n*n; i++){
        //     fprintf(fp1, "%lf ", rhs[i]);
        // }
        // Writing 2D solution back
        // FILE *fp1 = fopen("sol.txt", "w");
        for (int i = 0; i < n; i++){
            printf("%lf ", rhs[i]);
        }
        for (int j = 0; j < _n; j++) {
            for (int i = 0; i < _n; i++) {
                fprintf(fp1, "%lf ", rhs[i * _n + j]);
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

// double* say_hi(double *a){
//     printf("hi, Solution");
//     for(size_t i = 0; i < 2; i++)
//         a[i] = a[i] + 1;
//     return a;
// }




// // Writing 2D solution back
// FILE *fp1 = fopen("sol.txt", "w");
// for (int j = 0; j < n; j++) {
//     for (int i = 0; i < n; i++) {
//         fprintf(fp1, "%lf ", rhs[i * n + j]);
//     }
//     fprintf(fp1, "\n");
// }
// fclose(fp1);
