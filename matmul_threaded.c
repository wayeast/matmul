#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>


#define LG 50000
#define SM 25
#define ROWS 5
#define COLS 4
#define M 7
#define N 4
#define P 9
#define NO_THREADS 3


/************************************************
 * Auxiliary function and struct declarations
 * **********************************************
 */
double *allocatearray(int);
double **allocatetable(int, int);
void rand_pop_array(double *, int);
void rand_pop_table(double **, int, int);
void work_array(int);
void work_table(int, int);
void sum_with_threads();
void thread_matmul(int, int, int, int);
void *thread_sum_array(void *);
void *thread_mult_tables(void *);
int **assign_rows(int, int);
void free_table(double **, int);
void print_table(double **, int, int);
void free_row_assts(int **, int);

struct matmul_args {
    // arguments for multiplications of matrices of
    // dimensions (mxn)(nxp)
    double **A_mat;
    double **B_mat;
    double **C_mat;
    int *C_rows;
    int n;
    int p;
    int thread_id;
};

struct sum_array_args {
    double *array;
    double *sum;
    int no_elems;
    int thread_id;
};
/***************************************************
 * End auxiliary declarations
 * *************************************************
*/


/*
 * *************************************
 * Main routine
 * *************************************
 */
int main(int argc, char **argv) {
    if (argc != 5) {
        printf("Usage: <m_value> <n_value> <p_value> <numThreads>\n");
        return 1;
    }
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int p = atoi(argv[3]);
    int numThreads = atoi(argv[4]);
    //work_array(SM);
    //work_table(ROWS, COLS);
    //sum_with_threads();
    thread_matmul(m, n, p, numThreads);
    return 0;
}


/*
 * *************************************************
 * Accessory function definitions
 * *************************************************
 */

void thread_matmul(int m, int n, int p, int no_threads) {
    int i, j;

    // create and initialize A, B, C
    printf("Creating tables A, B, C\n");
    double **A, **B, **C;
    A = allocatetable(m, n);
    B = allocatetable(n, p);
    rand_pop_table(A, m, n);
    rand_pop_table(B, n, p);
    C = allocatetable(m, p);
    for (i=0; i<m; i++) {
        for (j=0; j<p; j++) {
            C[i][j] = 0.0;
        }
    }

    // assign rows to (eventual) threads
    printf("Creating thread row assignments\n");
    int **row_nos;
    row_nos = assign_rows(m, no_threads);

    // create thread arguments array
    printf("Creating thread arguments array\n");
    struct matmul_args *thread_args;
    thread_args = (struct matmul_args *) malloc((unsigned) no_threads * sizeof(struct matmul_args));
    for (i=0; i<no_threads; i++) {
        thread_args[i].A_mat = A;
        thread_args[i].B_mat = B;
        thread_args[i].C_mat = C;
        thread_args[i].C_rows = row_nos[i];
        thread_args[i].n = n;
        thread_args[i].p = p;
        thread_args[i].thread_id = i;
    }

    // create threads on stack
    pthread_t threads[no_threads];
    for (i=0; i<no_threads; i++) {
        printf("Creating thread %d\n", i);
        pthread_create(&threads[i], NULL, thread_mult_tables, (void *) &thread_args[i]);
    }
    for (i=0; i<no_threads; i++) {
        printf("Joining thread %d\n", i);
        pthread_join(threads[i], NULL);
    }

    printf("\nTable A:\n");
    print_table(A, m, n);
    printf("Table B:\n");
    print_table(B, n, p);
    printf("Table C:\n");
    print_table(C, m, p);

    printf("Cleaning up\n");
    free_table(A, m);
    free_table(B, n);
    free_table(C, m);
    free_row_assts(row_nos, no_threads);
    free(thread_args);
}  // end thread_matmul

void print_table(double **t, int rows, int cols){
    int i, j;
    for (i=0; i<rows; i++){
        for (j=0; j<cols; j++) {
            printf("%-11.8f", t[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void free_table(double **t, int rows) {
    int i;
    for (i=0; i<rows; i++) {
        free(t[i]);
    }
    free(t);
    return;
}

void free_row_assts(int **rat, int threads){
    int i;
    for (i=0; i<threads; i++) free(rat[i]);
    free(rat);
    return;
}

int **assign_rows(int no_rows, int no_threads) {
    // when called, assts should be an unassigned pointer
    int i, j;
    int **assts;

    // create array of row counts for each thread
    printf("In assign_rows, creating row count array on stack.\n");
    int rows_per_thread[no_threads];
    for (i=0; i<no_threads; i++) rows_per_thread[i] = no_rows / no_threads;
    for (i=0; i<no_rows%no_threads; i++) rows_per_thread[i]++;
    for (i=0; i<no_threads; i++) printf("%d\n", rows_per_thread[i]);

    // create an assignment table of pointers: no_threads rows x
    // rows_per_thread[i] + 1 columns.  Initialize all elements to -1
    printf("In assign_rows, initializing assignments table.\n");
    assts = (int **) malloc((unsigned) no_threads * sizeof(int*));
    for (i=0; i<no_threads; i++) assts[i] = (int *) malloc((unsigned) (rows_per_thread[i] + 1) * sizeof(int));
    for (i=0; i<no_threads; i++) {
        for (j=0; j<rows_per_thread[i] + 1; j++){
            //printf("Initializing assts[%d][%d] to -1\n", i, j);  // testing
            assts[i][j] = -1;
        }
    }
    // populate assignment table and attach row index value to each assts[i][j]
    // pointer
    printf("In assign_rows, populating row assignment table.\n");
    for (i=0, j=-1; i<no_rows; i++) {
        if (i%no_threads == 0) j++;
        assts[i%no_threads][j] = i;
    }
    printf("Exiting assign_rows.\n");
    /*
     * For testing
    for (i=0; i<no_threads; i++){
        printf("Thread %d: ", i);
        for (j=0; j<rows_per_thread[i]; j++){
            printf("%-5d", assts[i][j]);
        }
        printf("\n");
    }
    */
    return assts;
}

void *thread_mult_tables(void *argument) {
    struct matmul_args *args = (struct matmul_args *) argument;
    int i, k, row, col;
    int *C_row;

    i=0;
    while (args->C_rows[i] >= 0) {
        row = args->C_rows[i];
        printf("Thread %d calculating row %d\n", args->thread_id, row);
        for (col=0; col<args->p; col++) {
            for (k=0; k<args->n; k++) {
                args->C_mat[row][col] += args->A_mat[row][k] * args->B_mat[k][col];
            }
        }
        i++;

    }
}

void sum_with_threads() {
    int i;

    // create and populate arrays to sum
    double *arrays[NO_THREADS];
    for (i=0; i<NO_THREADS; i++){
        arrays[i] = allocatearray(LG);
        rand_pop_array(arrays[i], LG);
    }
    printf("Finished creating data arrays.\n");

    // prepare thread arguments
    double *sums;
    sums = (double *) malloc((unsigned) NO_THREADS * sizeof(double));
    struct sum_array_args *thread_args;
    thread_args = (struct sum_array_args *) malloc((unsigned) NO_THREADS * sizeof(struct sum_array_args));
    for (i=0; i<NO_THREADS; i++) {
        thread_args[i].array = arrays[i];
        thread_args[i].sum = &sums[i];
        thread_args[i].no_elems = LG;
        thread_args[i].thread_id = i;
    }

    // launch threads
    pthread_t threads[NO_THREADS];
    for (i=0; i<NO_THREADS; i++) {
        printf("Creating thread %d...\n", i);
        pthread_create(&threads[i], NULL, thread_sum_array, (void *) &thread_args[i]);
    }
    // join threads
    for (i=0; i<NO_THREADS; i++) {
        pthread_join(threads[i], NULL);
        printf("Joined thread %d\n", i);
    }
    
    // free allocated memory
    for (i=0; i<NO_THREADS; i++) {
        free(arrays[i]);
    }
    free(sums);
    free(thread_args);

    printf("Made it to the end.\n");
    return;
}

void *thread_sum_array(void *argument) {
    printf("Entering thread_sum_array.\n");
    int i;
    struct sum_array_args *args = (struct sum_array_args *) argument;

    // initialize sum to 0.0
    *(args->sum) = 0.0;
    // sum elements of array
    for (i=0; i<args->no_elems; i++) {
        *(args->sum) += (args->array)[i];
    }
    printf("Thread %d exiting after summing %f\n", args->thread_id, *(args->sum));

}

void work_table(int rows, int cols){
    printf("Entering work_table.  Working with table size %dx%d\n", rows, cols);
    double **tab;
    int r, c;
    tab = allocatetable(rows, cols);
    rand_pop_table(tab, rows, cols);

    for (r=0; r<rows; r++){
        for (c=0; c<cols; c++){
            printf("Row %d, Col %d : %.6f\n", r, c, tab[r][c]);
        }
    }
    printf("Freeing table...\n");
    for (r=0; r<rows; r++) {
        free(tab[r]);
    }
    free(tab);
    return;
}

void work_array(int size) {
    int i;
    printf("Entering work_array.  Working with array size %d.\n", size);

    double *array;
    array = allocatearray(size);
    rand_pop_array(array, size);

    for (i=0; i<size; i++){
        printf("Element %d : %.6f\n", i, array[i]);
    }

    printf("Freeing array and exiting work_array...\n");
    free(array);
    return;
}

double *allocatearray(int elements) {
    printf("Entering allocatearray...\n");
    double *ptr;
    ptr = (double *) malloc((unsigned) elements * sizeof(double));
    if (ptr != 0) printf("Allocation successful.  Exiting allocatearray.\n");
    return ptr;
}

double **allocatetable(int rows, int cols){
    printf("Entering allocatetable.\n");
    double **ptr;
    int i;
    ptr = (double **) malloc((unsigned) rows * sizeof(double *));
    if (ptr != 0) {
        for (i=0; i<rows; i++) {
            ptr[i] = (double *) malloc((unsigned) cols * sizeof(double));
            if (ptr[i] == 0) {
                printf("Error allocating row pointer %d\n", i);
                return NULL;
            }
        }
    }
    else {
        printf("Error allocating table pointer.\n");
        return NULL;
    }
    printf("Allocated %dx%d table.\n", rows, cols);
    return ptr;
}

void rand_pop_array(double *a, int s) {
    printf("Populating array of size %d\n", s);
    int i;
    for (i=0; i<s; i++){
        a[i] = (double) rand() / (double) RAND_MAX;
    }
}

void rand_pop_table(double **t, int rows, int cols) {
    printf("Populating table size %dx%d\n", rows, cols);
    int r, c;
    for (r=0; r<rows; r++){
        for (c=0; c<cols; c++) {
            t[r][c] = (double) rand() / (double) RAND_MAX;
        }
    }
}
