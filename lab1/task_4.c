#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define __USE_XOPEN
#include <math.h>
#include <string.h>

#define MAX_ITER 1000
#define TOLERANCE 1e-6

double f_function(double x, double y) {
    return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            printf("Usage: %s <grid_size>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int N = atoi(argv[1]);
    double c = 0.0;
    double h = 1.0 / (N + 1);

    if (N % size != 0) {
        if (rank == 0) {
            printf("Grid size must be divisible by the number of processes\n");
        }
        MPI_Finalize();
        return 1;
    }

    // Распределение данных
    int local_N = N / size;
    int start_row_global = rank * local_N;
    int end_row_global = start_row_global + local_N - 1;

    // Выделение памяти (+2 для граничных условий)
    double **u = (double **)malloc((local_N + 2) * sizeof(double *));
    double **u_new = (double **)malloc((local_N + 2) * sizeof(double *));
    for (int i = 0; i < local_N + 2; i++) {
        u[i] = (double *)malloc((N + 2) * sizeof(double));
        u_new[i] = (double *)malloc((N + 2) * sizeof(double));
    }

    // Инициализация - все точки включая границы
    for (int i = 0; i < local_N + 2; i++) {
        for (int j = 0; j < N + 2; j++) {
            u[i][j] = c;
            u_new[i][j] = c;
        }
    }

    // Буферы для обмена
    double *send_top = (double *)malloc((N + 2) * sizeof(double));
    double *send_bottom = (double *)malloc((N + 2) * sizeof(double));
    double *recv_top = (double *)malloc((N + 2) * sizeof(double));
    double *recv_bottom = (double *)malloc((N + 2) * sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // Основной итерационный цикл
    for (int iter = 0; iter < MAX_ITER; iter++) {
        
        // Обмен граничными слоями
        if (rank > 0) {
            MPI_Sendrecv(u[1], N + 2, MPI_DOUBLE, rank - 1, 0,
                        recv_top, N + 2, MPI_DOUBLE, rank - 1, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            memcpy(u[0], recv_top, (N + 2) * sizeof(double));
        }
        
        if (rank < size - 1) {
            MPI_Sendrecv(u[local_N], N + 2, MPI_DOUBLE, rank + 1, 0,
                        recv_bottom, N + 2, MPI_DOUBLE, rank + 1, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            memcpy(u[local_N + 1], recv_bottom, (N + 2) * sizeof(double));
        }

        // Волновая схема вычислений
        for (int s = 2; s <= N + local_N; s++) {
            for (int i_local = 1; i_local <= local_N; i_local++) {
                int i_global = start_row_global + i_local;
                int j_global = s - i_global;
                
                if (j_global >= 1 && j_global <= N) {
                    double x = j_global * h;
                    double y = i_global * h;
                    double f_val = f_function(x, y);
                    
                    // Метод Гаусса-Зейделя
                    u_new[i_local][j_global] = 0.25 * (
                        u_new[i_local-1][j_global] + u[i_local+1][j_global] +
                        u_new[i_local][j_global-1] + u[i_local][j_global+1] -
                        h * h * f_val
                    );
                }
            }
        }

        // Проверка сходимости
        double local_max_diff = 0.0;
        for (int i = 1; i <= local_N; i++) {
            for (int j = 1; j <= N; j++) {
                double diff = fabs(u_new[i][j] - u[i][j]);
                if (diff > local_max_diff) {
                    local_max_diff = diff;
                }
                u[i][j] = u_new[i][j]; // Обновление значений
            }
        }

        double global_max_diff;
        MPI_Allreduce(&local_max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        if (global_max_diff < TOLERANCE) {
            if (rank == 0) {
                printf("Convergence reached after %d iterations\n", iter + 1);
            }
            break;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Grid: %d, Processes: %d, Time: %.6f seconds\n", 
               N, size, end_time - start_time);
    }

    // Освобождение памяти
    for (int i = 0; i < local_N + 2; i++) {
        free(u[i]);
        free(u_new[i]);
    }
    free(u);
    free(u_new);
    free(send_top);
    free(send_bottom);
    free(recv_top);
    free(recv_bottom);

    MPI_Finalize();
    return 0;
}