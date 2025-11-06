#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>

// #define DEBUG

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
    printf(
        "\n##################################################################################\n");
}

double MatVecInRowMode(double *mat, double *vec, double *res, const int rows, const int cols,
                       const int my_rank, const int comm_sz) {

    int rows_per_process[comm_sz];
    int displs_row[comm_sz]; // смещения в строках
    int base = rows / comm_sz;
    int rem = rows % comm_sz;
    displs_row[0] = 0;
    for (int p = 0; p < comm_sz; ++p) {
        rows_per_process[p] = base + (p < rem ? 1 : 0);
        if (p > 0) {
            displs_row[p] = displs_row[p - 1] + rows_per_process[p - 1];
        }
    }

    int local_rows = rows_per_process[my_rank];
    double *local_mat = calloc(local_rows * cols, sizeof(double));

    int sendcounts_elem[comm_sz];
    int displs_elem[comm_sz];
    for (int p = 0; p < comm_sz; ++p) {
        sendcounts_elem[p] = rows_per_process[p] * cols;
    }
    displs_elem[0] = 0;
    for (int p = 1; p < comm_sz; ++p) {
        displs_elem[p] = displs_elem[p - 1] + sendcounts_elem[p - 1];
    }

    double* vec_bcast = vec;
    if (my_rank != kRoot) {
        vec_bcast = calloc(cols, sizeof(double));
    }
    
    double *local_res = calloc(local_rows, sizeof(double));
    
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    MPI_Request requests[2];
    
    MPI_Iscatterv(mat, sendcounts_elem, displs_elem, MPI_DOUBLE,
                  local_mat, sendcounts_elem[my_rank], MPI_DOUBLE,
                  kRoot, MPI_COMM_WORLD, &requests[0]);

    MPI_Ibcast(vec_bcast, cols, MPI_DOUBLE, kRoot, MPI_COMM_WORLD, &requests[1]);
    
    MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

    for (int row = 0; row < local_rows; ++row) {
        double sum = 0.0;
        double *mat_row = local_mat + (uint64_t)row * cols;
        for (int col = 0; col < cols; ++col) {
            sum += mat_row[col] * vec_bcast[col];
        }
        local_res[row] = sum;
    }

    MPI_Gatherv(local_res, local_rows, MPI_DOUBLE, 
                res, rows_per_process, displs_row, MPI_DOUBLE,
                kRoot, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    free(local_mat);
    free(local_res);
    if (my_rank != kRoot) {
        free(vec_bcast);
    }

    double local_time = t_end - t_start;
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, kRoot, MPI_COMM_WORLD);
    return max_time;
}

double MatVecInColMode(double *mat, double *vec, double *res, const int rows, const int cols,
                       const int my_rank, const int comm_sz) {
    int base = cols / comm_sz;
    int rem = cols % comm_sz;
    int *cols_per_process = malloc(comm_sz * sizeof(int));
    int *displs_col = malloc(comm_sz * sizeof(int));
    displs_col[0] = 0;
    for (int p = 0; p < comm_sz; ++p) {
        cols_per_process[p] = base + (p < rem ? 1 : 0);
        if (p > 0) {
            displs_col[p] = displs_col[p - 1] + cols_per_process[p - 1];
        }
    }
    int local_cols = cols_per_process[my_rank];
    int start_col = displs_col[my_rank];

    size_t local_mat_elems = (size_t)rows * local_cols;
    double *local_mat;
    if (local_mat_elems == 0) {
        local_mat = malloc(sizeof(double));
    } else {
        local_mat = calloc(local_mat_elems, sizeof(double));
    }

    double *local_vec;
    if (local_cols == 0) local_vec = malloc(sizeof(double));
    else local_vec = calloc(local_cols, sizeof(double));

    double *local_res;
    if (rows == 0) local_res = malloc(sizeof(double));
    else local_res = calloc(rows, sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    if (my_rank == kRoot) {
        MPI_Request *send_requests = NULL;
        MPI_Datatype *col_block_types = NULL;
        int request_count = 0;

        if (comm_sz > 1) {
            send_requests = malloc((comm_sz - 1) * sizeof(MPI_Request));
            col_block_types = malloc((comm_sz - 1) * sizeof(MPI_Datatype));
        }

        for (int p = 0; p < comm_sz; ++p) {
            if (p == kRoot) continue;

            int gsizes[2] = {rows, cols};
            int subsizes[2] = {rows, cols_per_process[p]};
            int starts[2] = {0, displs_col[p]};

            MPI_Type_create_subarray(2, gsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &col_block_types[request_count]);
            MPI_Type_commit(&col_block_types[request_count]);

            MPI_Isend(mat, 1, col_block_types[request_count], p, 0, MPI_COMM_WORLD, &send_requests[request_count]);
            request_count++;
        }

        for (int i = 0; i < rows; ++i) {
            memcpy(local_mat + (uint64_t)i * local_cols,
                   mat + (uint64_t)i * cols + start_col,
                   (size_t)local_cols * sizeof(double));
        }

        if (request_count > 0) {
            MPI_Waitall(request_count, send_requests, MPI_STATUSES_IGNORE);
            for (int i = 0; i < request_count; ++i) {
                MPI_Type_free(&col_block_types[i]);
            }
            free(send_requests);
            free(col_block_types);
        }

    } else {
        MPI_Recv(local_mat, rows * local_cols, MPI_DOUBLE, kRoot, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Scatterv(vec, cols_per_process, displs_col, MPI_DOUBLE,
                 local_vec, local_cols, MPI_DOUBLE,
                 kRoot, MPI_COMM_WORLD);

    for (int row = 0; row < rows; ++row) {
        double sum = 0.0;
        double *mat_row = local_mat + (uint64_t)row * local_cols;
        for (int col = 0; col < local_cols; ++col) {
            sum += mat_row[col] * local_vec[col];
        }
        local_res[row] = sum;
    }

    MPI_Reduce(local_res, res, rows, MPI_DOUBLE, MPI_SUM, kRoot, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    free(local_mat);
    free(local_vec);
    free(local_res);
    free(cols_per_process);
    free(displs_col);

    double local_time = t_end - t_start;
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, kRoot, MPI_COMM_WORLD);
    return max_time;
}

double MatVecInBlockMode(double *mat, double *vec, double *res, const int rows, const int cols,
                         const int my_rank, const int comm_sz) {

    int dims[2] = {0, 0};
    MPI_Dims_create(comm_sz, 2, dims);
    int p_row = dims[0];
    int p_col = dims[1];

    MPI_Comm grid_comm;
    int periods[2] = {0, 0};
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &grid_comm);

    int grid_rank = -1;
    int grid_coords[2] = {0, 0};
    MPI_Comm_rank(grid_comm, &grid_rank);
    MPI_Cart_coords(grid_comm, grid_rank, 2, grid_coords);
    int my_row = grid_coords[0];
    int my_col = grid_coords[1];

    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();

    int base_rows = rows / p_row;
    int rem_rows = rows % p_row;
    int local_rows = base_rows + (my_row < rem_rows ? 1 : 0);
    int start_row = 0;
    for (int i = 0; i < my_row; ++i) {
        start_row += base_rows + (i < rem_rows ? 1 : 0);
    }

    int base_cols = cols / p_col;
    int rem_cols = cols % p_col;
    int local_cols = base_cols + (my_col < rem_cols ? 1 : 0);
    int start_col = 0;
    for (int j = 0; j < my_col; ++j) {
        start_col += base_cols + (j < rem_cols ? 1 : 0);
    }

    size_t local_mat_elems = (size_t)local_rows * local_cols;
    double *local_mat;
    if (local_mat_elems == 0) local_mat = malloc(sizeof(double)); else local_mat = calloc(local_mat_elems, sizeof(double));
    double *local_vec;
    if (local_cols == 0) local_vec = malloc(sizeof(double)); else local_vec = calloc(local_cols, sizeof(double));

    if (my_rank == kRoot) {
        MPI_Request *requests = NULL;
        MPI_Datatype *block_types = NULL;
        int request_count = 0;

        if (comm_sz > 1) {
            requests = malloc((comm_sz - 1) * sizeof(MPI_Request));
            block_types = malloc((comm_sz - 1) * sizeof(MPI_Datatype));
        }

        for (int i = 0; i < p_row; ++i) {
            for (int j = 0; j < p_col; ++j) {
                int dest_coords[2] = {i, j};
                int dest_rank;
                MPI_Cart_rank(grid_comm, dest_coords, &dest_rank);

                if (dest_rank == kRoot) continue;

                int block_rows = base_rows + (i < rem_rows ? 1 : 0);
                int block_start_row = 0;
                for (int r = 0; r < i; ++r) {
                    block_start_row += base_rows + (r < rem_rows ? 1 : 0);
                }

                int block_cols = base_cols + (j < rem_cols ? 1 : 0);
                int block_start_col = 0;
                for (int c = 0; c < j; ++c) {
                    block_start_col += base_cols + (c < rem_cols ? 1 : 0);
                }

                int gsizes[2] = {rows, cols};
                int subsizes[2] = {block_rows, block_cols};
                int starts[2] = {block_start_row, block_start_col};

                MPI_Type_create_subarray(2, gsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &block_types[request_count]);
                MPI_Type_commit(&block_types[request_count]);

                MPI_Isend(mat, 1, block_types[request_count], dest_rank, 0, grid_comm, &requests[request_count]);
                request_count++;
            }
        }

        for (int r = 0; r < local_rows; ++r) {
            memcpy(local_mat + r * local_cols, mat + (start_row + r) * cols + start_col, (size_t)local_cols * sizeof(double));
        }

        if (request_count > 0) {
            MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
            for (int i = 0; i < request_count; ++i) {
                MPI_Type_free(&block_types[i]);
            }
            free(requests);
            free(block_types);
        }

    } else {
        MPI_Recv(local_mat, local_rows * local_cols, MPI_DOUBLE, kRoot, 0, grid_comm, MPI_STATUS_IGNORE);
    }

    MPI_Comm col_comm;
    MPI_Comm_split(grid_comm, my_col, my_row, &col_comm);

    MPI_Comm leaders_comm = MPI_COMM_NULL;
    int color = (my_row == 0) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(grid_comm, color, my_col, &leaders_comm);

    if (leaders_comm != MPI_COMM_NULL) {
        int leaders_rank;
        MPI_Comm_rank(leaders_comm, &leaders_rank);

        int *sendcounts = NULL;
        int *displs = NULL;

        if (leaders_rank == 0) {
            sendcounts = malloc(p_col * sizeof(int));
            displs = malloc(p_col * sizeof(int));
            int current_displ = 0;
            for (int j = 0; j < p_col; ++j) {
                sendcounts[j] = base_cols + (j < rem_cols ? 1 : 0);
                displs[j] = current_displ;
                current_displ += sendcounts[j];
            }
        }

        MPI_Scatterv(vec, sendcounts, displs, MPI_DOUBLE,
                     local_vec, local_cols, MPI_DOUBLE,
                     0, leaders_comm);

        if (leaders_rank == 0) {
            free(sendcounts);
            free(displs);
        }
        MPI_Comm_free(&leaders_comm);
    }

    MPI_Bcast(local_vec, local_cols, MPI_DOUBLE, 0, col_comm);

    double *local_res = calloc(local_rows ? local_rows : 1, sizeof(double));
    for (int i = 0; i < local_rows; ++i) {
        for (int j = 0; j < local_cols; ++j) {
            local_res[i] += local_mat[i * local_cols + j] * local_vec[j];
        }
    }

    MPI_Comm row_comm;
    MPI_Comm_split(grid_comm, my_row, my_col, &row_comm);

    double *row_res = NULL;
    if (my_col == 0) {
        row_res = calloc(local_rows ? local_rows : 1, sizeof(double));
    }

    MPI_Reduce(local_res, row_res, local_rows, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    if (my_col == 0) {
        int *recvcounts = NULL;
        int *displs = NULL;

        if (my_rank == kRoot) {
            recvcounts = malloc(p_row * sizeof(int));
            displs = malloc(p_row * sizeof(int));
            displs[0] = 0;
            for (int i = 0; i < p_row; ++i) {
                recvcounts[i] = base_rows + (i < rem_rows ? 1 : 0);
                if (i > 0) {
                    displs[i] = displs[i - 1] + recvcounts[i - 1];
                }
            }
        }

        MPI_Gatherv(row_res, local_rows, MPI_DOUBLE,
                    res, recvcounts, displs, MPI_DOUBLE,
                    0, col_comm);

        if (my_rank == kRoot) {
            free(recvcounts);
            free(displs);
        }
        free(row_res);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    free(local_mat);
    free(local_vec);
    free(local_res);

    if (row_comm != MPI_COMM_NULL) MPI_Comm_free(&row_comm);
    if (col_comm != MPI_COMM_NULL) MPI_Comm_free(&col_comm);
    if (grid_comm != MPI_COMM_NULL) MPI_Comm_free(&grid_comm);

    double local_time = t_end - t_start;
    double max_time = 0.0;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, kRoot, MPI_COMM_WORLD);
    return max_time;
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
                fprintf(stderr,
                        "Invalid argument: incorrect mode. Mod can only have following values: "
                        "row, col or block.\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    } else {
        if (my_rank == kRoot) {
            fprintf(stderr,
                    "Invalid argument: mode is required. Mod can only have following values: row, "
                    "col or block.\n");
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    srand(time(NULL));
    int rows = rand() % 10 + 1, cols = rand() % 10 + 1;
    if (argc >= 3) {
        for (uint32_t i = 0; i < strlen(argv[2]); ++i) {
            if (!isdigit(argv[2][i])) {
                if (my_rank == kRoot) {
                    fprintf(stderr,
                            "Invalid argument: expect positive integer number - dimension of the "
                            "matrix or rows and cols.\n");
                }
                MPI_Finalize();
                return EXIT_FAILURE;
            }
        }

        char *endptr = NULL;
        rows = cols = atoi(argv[2]);
        if (rows == 0) {
            if (my_rank == kRoot) {
                fprintf(stderr,
                        "Invalid argument: expect positive integer number - dimension of the "
                        "matrix or rows and cols.\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    }
    if (argc >= 4) {
        for (uint32_t i = 0; i < strlen(argv[3]); ++i) {
            if (!isdigit(argv[3][i])) {
                if (my_rank == kRoot) {
                    fprintf(stderr,
                            "Invalid argument: expect positive integer number - dimension of the "
                            "matrix or rows and cols.\n");
                }
                MPI_Finalize();
                return EXIT_FAILURE;
            }
        }

        char *endptr = NULL;
        cols = atoi(argv[3]);
        if (rows == 0) {
            if (my_rank == kRoot) {
                fprintf(stderr,
                        "Invalid argument: expect positive integer number - dimension of the "
                        "matrix or rows and cols.\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    }

    double *vec = NULL, *mat = NULL;
    MPI_Bcast(&rows, 1, MPI_INT, kRoot, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, kRoot, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

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
