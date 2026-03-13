#include <stdio.h>
#include <stdlib.h>
#include "hamiltonian.h"

double potential(double x) {
    return 0.0;
}

double coeff_a(double x) {
    return (x < 0.5) ? 1.0 : 5.0;
}

double **mass(double *x, int N) {
    double **m = alloc_matrix(N);
    int i;

    // i = 0
    m[0][0] = (x[1] - 0.0) / 3.0;
    m[0][1] = (x[0] - 0.0) / 6.0;

    // i = 1, ..., N-2
    for (i = 1; i < N-1; ++i) {
        m[i][i-1] = (x[i  ] - x[i-1]) / 6.0;
        m[i][i  ] = (x[i+1] - x[i-1]) / 3.0;
        m[i][i+1] = (x[i+1] - x[i  ]) / 6.0;
    }

    // i = N-1
    m[N-1][N-2] = (x[N-1] - x[N-2]) / 6.0;
    m[N-1][N-1] = (1.0    - x[N-2]) / 3.0;

    return m;
}

double **hamiltonian(double *x, int N) {
    double **h = alloc_matrix(N);
    int i;

    // i = 0
    h[0][0] =  coeff_a( x[0]         / 2.0) /  x[0]
              +coeff_a((x[1] + x[0]) / 2.0) / (x[1] - x[0]);
    h[0][1] = -coeff_a((x[1] + x[0]) / 2.0) / (x[1] - x[0]);

    h[0][0] +=  potential( x[0]         / 2.0) *  x[0]         / 3.0
               +potential((x[1] + x[0]) / 2.0) * (x[1] - x[0]) / 3.0;
    h[0][1] +=  potential((x[1] + x[0]) / 2.0) * (x[1] - x[0]) / 6.0;

    // i = 1, ..., N-2
    for (i = 1; i < N-1; ++i) {
        h[i][i-1] = -coeff_a((x[i  ] + x[i-1]) / 2.0) / (x[i  ] - x[i-1]);
        h[i][i  ] =  coeff_a((x[i  ] + x[i-1]) / 2.0) / (x[i  ] - x[i-1])
                    +coeff_a((x[i+1] + x[i  ]) / 2.0) / (x[i+1] - x[i  ]);
        h[i][i+1] = -coeff_a((x[i+1] + x[i  ]) / 2.0) / (x[i+1] - x[i  ]);

        h[i][i-1] +=  potential((x[i  ] + x[i-1]) / 2.0) * (x[i  ] - x[i-1]) / 6.0;
        h[i][i  ] +=  potential((x[i  ] + x[i-1]) / 2.0) * (x[i  ] - x[i-1]) / 3.0
                     +potential((x[i+1] + x[i  ]) / 2.0) * (x[i+1] - x[i  ]) / 3.0;
        h[i][i+1] +=  potential((x[i+1] + x[i  ]) / 2.0) * (x[i+1] - x[i  ]) / 6.0;

    }

    // i = N-1
    h[N-1][N-2] = -coeff_a((x[N-1] + x[N-2]) / 2.0) / (x[N-1] - x[N-2]);
    h[N-1][N-1] =  coeff_a((x[N-1] + x[N-2]) / 2.0) / (x[N-1] - x[N-2])
                  +coeff_a((1.0    + x[N-1]) / 2.0) / (1.0    - x[N-1]);

    h[N-1][N-2] +=  potential((x[N-1] + x[N-2]) / 2.0) * (x[N-1] - x[N-2]) / 6.0;
    h[N-1][N-1] +=  potential((x[N-1] + x[N-2]) / 2.0) * (x[N-1] - x[N-2]) / 3.0
                   +potential((1.0    + x[N-1]) / 2.0) * (1.0    - x[N-1]) / 3.0;

    return h;
}

double **alloc_matrix(int N) {
    double **p = malloc(N * sizeof(double *));
    int i;
    for (i = 0; i < N; ++i) p[i] = calloc(N, sizeof(double));
    return p;
}

void fprint_matrix(FILE *f, double **m, int N) {
    int i, j;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            fprintf(f, "%.16f ", m[i][j]);
        }
        fprintf(f, "\n");
    }
}

void free_matrix(double **m, int N) {
    int i;
    for (i = 0; i < N; ++i) {
        free(m[i]);
    }
    free(m);
}

double *read_nodes(const char *filename, int *N) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("file open");
        exit(1);
    }

    if (fscanf(f, "%d", N) != 1) {
        perror("must start with number of nodes");
        exit(1);
    }

    double *x = malloc(*N * sizeof(double));
    int i;
    for (i = 0; i < *N; ++i) {
        if (fscanf(f, "%lf", &x[i]) != 1) {
            perror("incorrect number of nodes");
            exit(1);
        }
        if (!(0.0 < x[i] && x[i] < 1.0 && (i == 0 || x[i] > x[i-1]))) {
            perror("mesh must be in increasing order, between 0.0 and 1.0");
            exit(1);
        }
    }

    fclose(f);
    return x;
}

int main() {
    int N;
    double *x = read_nodes("nodes.txt", &N);

    printf("Read %d nodes\n", N);

    double **m = mass(x, N);
    double **h = hamiltonian(x, N);

    FILE *m_file = fopen("M.txt", "w");
    FILE *h_file = fopen("H.txt", "w");

    fprint_matrix(m_file, m, N);
    fprint_matrix(h_file, h, N);

    free_matrix(m, N);
    free_matrix(h, N);
    return 0;
}
