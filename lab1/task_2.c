#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define DEBUG

const int kRoot = 0;

double NextDouble() {
    return (double)rand() / RAND_MAX * 2.0 - 1.0;
}

void InitVector(double *vec, const int size) {
    for (int i = 0; i < size; ++i) {
        vec[i] = NextDouble();
    }
}

void PrintMat(double *mat, const int rows, const int cols) {
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            printf("%.3f\t", mat[row * cols + col]);
        }
        printf("\n");
    }
    printf("##################################################################################\n");
}

void PrintVec(double *vec, const int size) {
    for (int row = 0; row < size; ++row) {
        printf("%.3f\t", vec[row]);
    }
    printf("\n##################################################################################\n");
}

double MatVecInRowMode(double *mat, double *vec, double *res, const int rows, const int cols,
                     const int my_rank, const int comm_sz) {
    int rows_per_process[comm_sz];
    int base = rows / comm_sz;
    int rem = rows % comm_sz;
    for (int process = 0; process < comm_sz; ++process) {
        rows_per_process[process] = base + (process < rem ? 1 : 0);
    }

    int local_rows = rows_per_process[my_rank];
    double *local_mat = calloc(local_rows * cols, sizeof(double));

    int sendcount[comm_sz];
    for (int process = 0; process < comm_sz; ++process) {
        sendcount[process] = rows_per_process[process] * cols;
    }

    int displs[comm_sz]; // смещения (в количестве элементов)
    displs[0] = 0;
    for (int process = 1; process < comm_sz; ++process) {
        displs[process] = displs[process - 1] + sendcount[process - 1];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    MPI_Scatterv(mat, sendcount, displs, MPI_DOUBLE, local_mat, sendcount[my_rank], MPI_DOUBLE,
                 kRoot, MPI_COMM_WORLD);
    if (my_rank == kRoot) {
        MPI_Bcast(vec, cols, MPI_DOUBLE, kRoot, MPI_COMM_WORLD);
    } else {
        vec = calloc(cols, sizeof(double));
        MPI_Bcast(vec, cols, MPI_DOUBLE, kRoot, MPI_COMM_WORLD);
    }

    double *local_res = calloc(local_rows, sizeof(double));
    for (int row = 0; row < local_rows; ++row) {
        double sum = 0.0;
        double *mat_row = local_mat + (uint64_t)row * (uint64_t)cols;
        for (int col = 0; col < cols; ++col) {
            sum += mat_row[col] * vec[col];
        }
        local_res[row] = sum;
    }

    int displs_row[comm_sz]; // смещения (в количестве строк)
    displs_row[0] = 0;
    for (int process = 1; process < comm_sz; ++process) {
        displs_row[process] = displs_row[process - 1] + rows_per_process[process - 1];
    }
    MPI_Gatherv(local_res, local_rows, MPI_DOUBLE, res, rows_per_process, displs_row, MPI_DOUBLE, kRoot, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    free(local_mat);
    free(local_res);
    if (my_rank != kRoot) {
        free(vec);
    }

    double local_time = t_end - t_start;
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return max_time;
}

double MatVecInColMode(double *mat, double *vec, double *res, const int rows, const int cols,
                     const int my_rank, const int comm_sz) {
    int base = cols / comm_sz;
    int rem = cols % comm_sz;
    int cols_per_process[comm_sz];
    for (int process = 0; process < comm_sz; ++process) {
        cols_per_process[process] = base + (process < rem ? 1 : 0);
    }
    int displs_col[comm_sz];
    displs_col[0] = 0;
    for (int process = 1; process < comm_sz; ++process) {
        displs_col[process] = displs_col[process - 1] + cols_per_process[process - 1];
    }
    int local_cols = cols_per_process[my_rank];

    double *sendbuf = NULL;
    int *sendcounts = NULL, *displs = NULL;
    if (my_rank == kRoot) {
        sendcounts = calloc(comm_sz, sizeof(int));
        displs = calloc(comm_sz, sizeof(int));
        for (int process = 0; process < comm_sz; ++process) {
            sendcounts[process] = cols_per_process[process] * rows;
        }
        displs[0] = 0;
        for (int process = 1; process < comm_sz; ++process) {
            displs[process] = displs[process - 1] + sendcounts[process - 1];
        }

        sendbuf = calloc(rows * cols, sizeof(double)); // храним столбцы как строки
        for (int process = 0; process < comm_sz; ++process) {
            int column_start = displs_col[process];
            int curr_local_cols = cols_per_process[process];
            double *dest = sendbuf + displs[process];
            for (int row = 0; row < rows; ++row) {
                const double *rowptr = mat + (uint64_t)row * cols + column_start;
                memcpy(dest + (uint64_t)row * curr_local_cols, rowptr, sizeof(double) * curr_local_cols);
            }
        }
#ifdef DEBUG
        PrintMat(sendbuf, cols, rows);
#endif
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    double *local_mat = calloc(rows * local_cols, sizeof(double));
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_DOUBLE, local_mat, rows * local_cols, MPI_DOUBLE, kRoot, MPI_COMM_WORLD);

    double *local_vec = calloc(local_cols, sizeof(double));
    MPI_Scatterv(vec, cols_per_process, displs_col, MPI_DOUBLE, local_vec, local_cols, MPI_DOUBLE, kRoot, MPI_COMM_WORLD);

    double *local_res = calloc(rows, sizeof(double));
    for (int row = 0; row < rows; ++row) {
        double sum = 0.0;
        double *mat_row = local_mat + (uint64_t)row * local_cols;
        for (int col = 0; col < local_cols; ++col) {
            sum += mat_row[col] * local_vec[col];
        }
        local_res[row] = sum;
    }

    MPI_Reduce(local_res, res, rows, MPI_DOUBLE, MPI_SUM, kRoot, MPI_COMM_WORLD);
    double t_end = MPI_Wtime();
    double total_time = t_end - t_start;

    if (my_rank == 0) {
        free(sendbuf);
        free(sendcounts); 
        free(displs);
    }
    free(local_mat);
    free(local_vec);
    free(local_res);

    double max_time = 0.0;
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, kRoot, MPI_COMM_WORLD);
    return max_time;
}

double MatVecInBlockMode(double *mat, double *vec, double *res, const int rows, const int cols,
                     const int my_rank, const int comm_sz) {
}

int main(int argc, char *argv[]) {

    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    enum Mode { ROW, COL, BLOCK } mode = ROW;
    const int equal = 0;
    if (argc >= 2) {
        if (strcmp(argv[1], "row") == equal) {
            mode = ROW;
        } else if (strcmp(argv[1], "col") == equal) {
            mode = COL;
        } else if (strcmp(argv[1], "block") == equal) {
            mode = BLOCK;
        } else {
            if (my_rank == kRoot) {
                fprintf(stderr, "Invalid argument: incorrect mode. Mod can only have following values: row, col or block.\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    } else {
        if (my_rank == kRoot) {
            fprintf(stderr, "Invalid argument: mode is required. Mod can only have following values: row, col or block.\n");
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    double *vec = NULL, *mat = NULL;
    int rows = rand() % 10, cols = rand() % 10;
    MPI_Bcast(&rows, 1, MPI_INT, kRoot, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, kRoot, MPI_COMM_WORLD);
    double *res = calloc(rows, sizeof(double));
    if (my_rank == kRoot) {
        vec = calloc(cols, sizeof(double));
        mat = calloc(rows * cols, sizeof(double));
        InitVector(vec, cols);
        InitVector(mat, rows * cols);
    }

    double total_time = 0;
    switch (mode) {
        case ROW:
            total_time = MatVecInRowMode(mat, vec, res, rows, cols, my_rank, comm_sz);
            break;
        case COL:
            total_time = MatVecInColMode(mat, vec, res, rows, cols, my_rank, comm_sz);
            break;
        case BLOCK:
            total_time = MatVecInBlockMode(mat, vec, res, rows, cols, my_rank, comm_sz);
            break;
    }

    if (my_rank == kRoot) {
#ifdef DEBUG
        PrintMat(mat, rows, cols);
        PrintVec(vec, cols);
        PrintVec(res, rows);
#endif
        printf("Total time = %f seconds\n", total_time);
    }

    MPI_Finalize();

    free(res);
    if (my_rank == kRoot) {
        free(vec);
        free(mat);
    }

    return 0;
}
