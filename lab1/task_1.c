#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>

double NextDouble() {
    return (double)rand() / RAND_MAX * 2.0 - 1.0;
}

int main(int argc, char *argv[]) {

    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    const int root = 0;

    uint32_t count_points = 200000000;  // 200'000'000
    if (argc >= 2) {
        for (uint32_t i = 0; i < strlen(argv[1]); ++i) {
            if (!isdigit(argv[1][i])) {
                if (my_rank == root) {
                    fprintf(stderr,
                            "Invalid argument: expect positive integer number - count points.\n");
                }
                MPI_Finalize();
                return EXIT_FAILURE;
            }
        }

        char *endptr = NULL;
        count_points = strtoul(argv[1], NULL, 10);
        if (count_points == 0) {
            if (my_rank == root) {
                fprintf(stderr,
                        "Invalid argument: expect positive integer number - count points.\n");
            }
            MPI_Finalize();
            return EXIT_FAILURE;
        }
    }

    uint32_t point_per_process = count_points / comm_sz;
    if (my_rank == root) {
        point_per_process += count_points % comm_sz;
    }

    uint32_t points_in_circle = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    double t_start = MPI_Wtime();
    for (uint32_t i = 0; i < point_per_process; ++i) {
        double x = NextDouble(), y = NextDouble();
        if (x * x + y * y <= 1.0) {
            ++points_in_circle;
        }
    }
    uint32_t total_points_in_circle = 0;
    MPI_Reduce(&points_in_circle, &total_points_in_circle, 1, MPI_UINT32_T, MPI_SUM, root,
               MPI_COMM_WORLD);
    double t_end = MPI_Wtime();

    if (my_rank == root) {
        double total_time = t_end - t_start;
        printf("Points in circle: %d\n", total_points_in_circle);
        printf("Total time = %f seconds\n", total_time);
        printf("pi = %.15Lf\n",
               4.0l * (long double)total_points_in_circle / (long double)count_points);
    }

    MPI_Finalize();

    return 0;
}
