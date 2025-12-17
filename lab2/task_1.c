// компиляция
// gcc task_1.c -o task_1 -fopenmp -lm

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>

struct Point {
    double x;
    double y;
};

int is_in_mandelbrot(double complex c, int max_iter) {
    double complex z = 0.0 + 0.0 * I;
    for (int i = 0; i < max_iter; ++i) {
        z = z * z + c;
        if ((creal(z) * creal(z) + cimag(z) * cimag(z)) > 4.0) {
            return i + 1;
        }
    }
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Использование: %s <nthreads> <npoints_side>\n", argv[0]);
        return 1;
    }

    int nthreads = atoi(argv[1]);
    int npoints_side = atoi(argv[2]);
    int max_iter = 1000;

    if (nthreads <= 0 || npoints_side <= 0) {
        fprintf(stderr, "Количество потоков и точек должно быть положительным числом.\n");
        return 1;
    }

    omp_set_num_threads(nthreads);

    double x_min = -2.0, x_max = 1.0;
    double y_min = -1.5, y_max = 1.5;

    double dx = (x_max - x_min) / npoints_side;
    double dy = (y_max - y_min) / npoints_side;

    int capacity = npoints_side * npoints_side;
    struct Point* mandelbrot_points = (struct Point*)malloc(capacity * sizeof(struct Point));
    if (mandelbrot_points == NULL) {
        fprintf(stderr, "Ошибка выделения памяти.\n");
        return 1;
    }
    int point_count = 0;

    double start_time, end_time;

    start_time = omp_get_wtime();

    #pragma omp parallel
    {
        struct Point* local_points = (struct Point*)malloc(capacity / nthreads * sizeof(struct Point) + 1);
        int local_point_count = 0;

        #pragma omp for
        for (int i = 0; i < npoints_side; ++i) {
            for (int j = 0; j < npoints_side; ++j) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double complex c = x + y * I;

                if (is_in_mandelbrot(c, max_iter) == 0) {
                    local_points[local_point_count].x = x;
                    local_points[local_point_count].y = y;
                    local_point_count++;
                }
            }
        }

        #pragma omp critical
        {
            for (int k = 0; k < local_point_count; ++k) {
                mandelbrot_points[point_count] = local_points[k];
                point_count++;
            }
        }
        free(local_points);
    }
    
    end_time = omp_get_wtime();

    printf("Время вычислений: %.4f секунд\n", end_time - start_time);


    FILE* outfile = fopen("mandelbrot_points.csv", "w");
    if (outfile == NULL) {
        fprintf(stderr, "Не удалось открыть файл для записи.\n");
        free(mandelbrot_points);
        return 1;
    }

    fprintf(outfile, "real,imaginary\n");
    for (int i = 0; i < point_count; ++i) {
        fprintf(outfile, "%f,%f\n", mandelbrot_points[i].x, mandelbrot_points[i].y);
    }

    fclose(outfile);
    free(mandelbrot_points);

    printf("Найдено %d точек. Координаты сохранены в mandelbrot_points.csv\n", point_count);

    return 0;
}