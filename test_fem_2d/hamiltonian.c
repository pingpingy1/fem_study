#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hamiltonian.h"

double potential(struct Point x) {
    return 50.0 * x.x + 100.0 * x.y; 
}

double coeff_a(struct Point x) {
    return 1.0;
}

// ====================================

struct Point *read_nodes(const char *filename, int *N_node, int *N_bdry) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("file open");
        exit(1);
    }

    fscanf(f, "%d %d", N_node, N_bdry); 

    struct Point *x = malloc((*N_node + *N_bdry) * sizeof(struct Point));
    int i;
    for (i = 0; i < *N_node + *N_bdry; ++i) {
        if (fscanf(f, "%lf %lf", &x[i].x, &x[i].y) != 2) {
            perror("Incorrect format");
            exit(1);
        }
    }

    fclose(f);
    return x;
}

struct Triangle *read_triangles(const char *filename, int *N_tri) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("file open");
        exit(1);
    }

    fscanf(f, "%d", N_tri);

    struct Triangle *t = malloc(*N_tri * sizeof(struct Triangle));

    int i, v1, v2, v3;
    for (i = 0; i < *N_tri; ++i) {
        if (fscanf(f, "%d %d %d", &v1, &v2, &v3) != 3) {
            perror("Invalid .ele format");
            exit(1);
        }
        
        t[i].v[0] = v1 - 1;
        t[i].v[1] = v2 - 1;
        t[i].v[2] = v3 - 1;
    }

    fclose(f);
    return t;
}

// ====================================

void fprint_matrix(FILE *f, double **m, int N) {
    int i, j;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            fprintf(f, "%.16f ", m[i][j]);
        }
        fprintf(f, "\n");
    }
}

double **alloc_matrix(int N) {
    double **p = malloc(N * sizeof(double *));
    int i;
    for (i = 0; i < N; ++i) {
        p[i] = calloc(N, sizeof(double));
    };
    return p;
}

void free_matrix(double **m, int N) {
    int i;
    for (i = 0; i < N; ++i) {
        free(m[i]);
    }
    free(m);
}

void fprint_triangles(const char *filename, struct Triangle *t, int N_tri, int N_node) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("file open");
        exit(1);
    }

    fprintf(f, "%d\n", N_tri);

    int i, j;
    for (i = 0; i < N_tri; ++i) {
        for (j = 0; j < 3; ++j) {
            fprintf(
                f,
                "%c%d ",
                (t[i].v[j] >= N_node) ? 'b' : 'n' ,
                t[i].v[j] + 1 - (t[i].v[j] >= N_node) * N_node
            );
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

// ====================================

double area(struct Triangle t, struct Point *nodes) {
    double a = 0.0;
    int i;
    for (i = 0; i < 3; ++i) {
        a += nodes[t.v[i]].x * (nodes[t.v[(i+1)%3]].y - nodes[t.v[(i+2)%3]].y);
    }
    return 0.5 * fabs(a);
}

struct Point center(struct Triangle t, struct Point *nodes) {
    struct Point c = {
        (nodes[t.v[0]].x + nodes[t.v[1]].x + nodes[t.v[2]].x) / 3.0,
        (nodes[t.v[0]].y + nodes[t.v[1]].y + nodes[t.v[2]].y) / 3.0
    };
    return c;
}

// ====================================

double **mass(
    struct Triangle *t,
    struct Point *nodes,
    int N_tri,
    int N_node
) {
    double **m = alloc_matrix(N_node);
    int i, j, k;
    double A;

    for (i = 0; i < N_tri; ++i) {
        A = area(t[i], nodes);
        
        // For each pair of indices, if both are less than N_node,
        // add \int \phi_i(x) \phi_j(x) to m.
        // i == j -> A / 6.0
        // i != j -> A / 12.0
        for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
                if (t[i].v[j] < N_node && t[i].v[k] < N_node) {
                    m[t[i].v[j]][t[i].v[k]] += A / ((j == k) ? 6.0 : 12.0);
                }
            }
        }
    }

    return m;
}

double **hamiltonian(struct Triangle *t, struct Point *nodes, int N_tri, int N_node) {
    double **h = alloc_matrix(N_node);
    int i, j, k;
    double A, bj, bk, cj, ck;
    struct Point c;

    for (i = 0; i < N_tri; ++i) {
        A = area(t[i], nodes);
        c = center(t[i], nodes);

        // For each pair of indices, if both are less than N_node,
        // add coeff_a(center) \int \grad\phi_i \cdot \grad\phi_j to h
        // With b_i = x_j - x_k and c_i = y_j - y_k,
        // \int \grad\phi_i \cdot \grad\phi_j = (b_i b_j + c_i c_j) / 4A
        // Also, add potential term: V(center) * \int\phi_i\phi_j
        for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
                if (t[i].v[j] < N_node && t[i].v[k] < N_node) {
                    bj = nodes[t[i].v[(j+1)%3]].x - nodes[t[i].v[(j+2)%3]].x;
                    bk = nodes[t[i].v[(k+1)%3]].x - nodes[t[i].v[(k+2)%3]].x;
                    
                    cj = nodes[t[i].v[(j+1)%3]].y - nodes[t[i].v[(j+2)%3]].y;
                    ck = nodes[t[i].v[(k+1)%3]].y - nodes[t[i].v[(k+2)%3]].y;

                    h[t[i].v[j]][t[i].v[k]] +=
                        coeff_a(c) * (bj * bk + cj * ck) / (4.0 * A) +
                        potential(c) * A / ((j == k) ? 6.0 : 12.0);
                }
            }
        }
    }
    return h;
}

// ====================================

int main() {
    int N_node, N_bdry, N_tri;

    struct Point *x = read_nodes("nodes.txt", &N_node, &N_bdry);
    printf("Read %d nodes and %d boundaries\n", N_node, N_bdry);

    struct Triangle *t = read_triangles("triangles.txt", &N_tri);
    printf("Read %d triangles\n", N_tri);
    
    printf("First triangle:\n");
    printf("%8.2f %8.2f\n", x[t[0].v[1]].x, x[t[0].v[1]].y);
    printf("%8.2f %8.2f\n", x[t[0].v[2]].x, x[t[0].v[2]].y);
    printf("%8.2f %8.2f\n", x[t[0].v[3]].x, x[t[0].v[3]].y);
    printf("Area: %8.2f\n", area(t[0], x));

    double **m = mass(t, x, N_tri, N_node);
    printf("Mass calculation successful\n");

    FILE *m_file = fopen("M.txt", "w");
    printf("Mass file open successful\n");

    fprint_matrix(m_file, m, N_node);
    printf("Mass print successful\n");

    free_matrix(m, N_node);
    printf("Mass free successful\n");
    
    double **h = hamiltonian(t, x, N_tri, N_node);
    printf("Hamiltonian calculation successful\n");

    FILE *h_file = fopen("H.txt", "w");
    printf("Hamiltonian file open successful\n");

    fprint_matrix(h_file, h, N_node);
    printf("Hamiltonian print successful\n");

    free_matrix(h, N_node);

    return 0;
}
