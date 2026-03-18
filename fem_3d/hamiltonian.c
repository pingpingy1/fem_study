#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "hamiltonian.h"

double potential(struct Node p) {
    return 0.0;
}

double coeff_a(struct Node p) {
    return 1.0;
}

// =============================

struct CSR_Matrix create_matrix(int N_node) {
    struct CSR_Matrix m;
    m.N_node = N_node;
    m.rows = malloc(N_node * sizeof(struct CSR_Row));

    int i;
    for (i = 0; i < N_node; ++i) {
        m.rows[i].cnt = 0;
        m.rows[i].capacity = 1;
        m.rows[i].entries = malloc(sizeof(struct CSR_Entry));
    }

    return m;
}

void add_entry(struct CSR_Matrix *m, int row, int col, double val) {
    assert(row >= 0 && row < m->N_node);

    if (m->rows[row].cnt == m->rows[row].capacity) {
        m->rows[row].capacity *= 2;
        m->rows[row].entries = realloc(
            m->rows[row].entries,
            m->rows[row].capacity * sizeof(struct CSR_Entry)
        );
        assert(m->rows[row].entries != NULL);
    }
    m->rows[row].entries[m->rows[row].cnt++] = (struct CSR_Entry) {col, val};
}

void free_matrix(struct CSR_Matrix m) {
    int i, j;
    for (i = 0; i < m.N_node; ++i) {
        free(m.rows[i].entries);
    }
    free(m.rows);
}

// =============================

struct Node *read_nodes(const char *filename, int *N_node, int *N_bdry) {
    FILE *f = fopen(filename, "r");
    assert(f);

    int N_tot, dim, attr, bdry_markers;
    int read_res = fscanf(f, "%d %d %d %d", &N_tot, &dim, &attr, &bdry_markers);
    assert(read_res == 4);
    assert(dim == 3);
    assert(attr == 0);
    assert(bdry_markers == 1);

    *N_node = 0;
    *N_bdry = 0;
    struct Node *nodes = malloc(N_tot * sizeof(struct Node));
    int i;

    for (i = 0; i < N_tot; ++i) {
        int idx, bdry;
        double x, y, z;

        read_res = fscanf(f, "%d %lf %lf %lf %d", &idx, &x, &y, &z, &bdry);
        assert(read_res == 5);

        nodes[i] = (struct Node) {x, y, z, !!bdry};
        if (bdry) *N_bdry += 1;
        else *N_node += 1;
    }

    fclose(f);
    return nodes;
}

struct Tetra *read_tetras(const char *filename, int *N_tetra, int N_tot) {
    FILE *f = fopen(filename, "r");
    assert(f);

    int nodes_per_tet, attr;
    int read_res = fscanf(f, "%d %d %d", N_tetra, &nodes_per_tet, &attr);
    assert(read_res == 3);
    assert(nodes_per_tet == 4);
    assert(attr == 0);

    struct Tetra *tetras = malloc(*N_tetra * sizeof(struct Tetra));
    int i;
    
    for (i = 0; i < *N_tetra; ++i) {
        int idx;
        
        read_res = fscanf(f, "%d %d %d %d %d",
            &idx,
            &tetras[i].v[0],
            &tetras[i].v[1],
            &tetras[i].v[2],
            &tetras[i].v[3]
        );
        assert(read_res == 5);
        
        int j;
        for (j = 0; j < 4; ++j) {
            assert(tetras[i].v[j] > 0 && tetras[i].v[j] <= N_tot);
            tetras[i].v[j] -= 1;
        }
    }

    fclose(f);
    return tetras;
}

struct Mesh read_mesh(const char *node_file, const char *tetra_file) {
    struct Mesh mesh;
    mesh.nodes = read_nodes(node_file, &mesh.N_node, &mesh.N_bdry);
    mesh.N_tot = mesh.N_node + mesh.N_bdry;
    mesh.tetras = read_tetras(tetra_file, &mesh.N_tetra, mesh.N_tot);

    mesh.node_to_idx = malloc(mesh.N_tot * sizeof(int));
    int i, idx = 0;
    for (i = 0; i < mesh.N_tot; ++i) {
        if (mesh.nodes[i].is_bdry) {
            mesh.node_to_idx[i] = -1;
        } else {
            mesh.node_to_idx[i] = idx++;
        }
    }

    return mesh;
}

void free_mesh(struct Mesh mesh) {
    free(mesh.nodes);
    free(mesh.tetras);
    free(mesh.node_to_idx);
}

// =============================

double volume(struct Tetra t, struct Node *nodes) {
    double ax = nodes[t.v[1]].x - nodes[t.v[0]].x,
           ay = nodes[t.v[1]].y - nodes[t.v[0]].y,
           az = nodes[t.v[1]].z - nodes[t.v[0]].z,
           bx = nodes[t.v[2]].x - nodes[t.v[0]].x,
           by = nodes[t.v[2]].y - nodes[t.v[0]].y,
           bz = nodes[t.v[2]].z - nodes[t.v[0]].z,
           cx = nodes[t.v[3]].x - nodes[t.v[0]].x,
           cy = nodes[t.v[3]].y - nodes[t.v[0]].y,
           cz = nodes[t.v[3]].z - nodes[t.v[0]].z;
    double sum =
        ax * by * cz +
        ay * bz * cx +
        az * bx * cy -
        ax * bz * cy -
        ay * bx * cz -
        az * by * cx;
    return fabs(sum) / 6.0;
}

struct Node center(struct Tetra t, struct Node *nodes) {
    return (struct Node) {
        (nodes[t.v[0]].x + nodes[t.v[1]].x + nodes[t.v[2]].x + nodes[t.v[3]].x) / 4.0,
        (nodes[t.v[0]].y + nodes[t.v[1]].y + nodes[t.v[2]].y + nodes[t.v[3]].y) / 4.0,
        (nodes[t.v[0]].z + nodes[t.v[1]].z + nodes[t.v[2]].z + nodes[t.v[3]].z) / 4.0,
        0
    };
}

// =============================

struct CSR_Matrix mass(struct Mesh mesh) {
    // Creates the mass matrix:
    // M_jk = \int dV \phi_j \phi_k
    // For adjacent j and k,
    // \int_tetra dV \phi_j \phi_k = (j == k) ? 0.6V : 0.3V

    struct CSR_Matrix m = create_matrix(mesh.N_node);

    int i;
    for (i = 0; i < mesh.N_tetra; ++i) {
        double V = volume(mesh.tetras[i], mesh.nodes);
        int j, k;
        for (j = 0; j < 4; ++j) {
            for (k = 0; k < 4; ++k) {
                if (!mesh.nodes[mesh.tetras[i].v[j]].is_bdry &&
                    !mesh.nodes[mesh.tetras[i].v[k]].is_bdry) {
                    add_entry(
                        &m,
                        mesh.node_to_idx[mesh.tetras[i].v[j]],
                        mesh.node_to_idx[mesh.tetras[i].v[k]],
                        (j == k) ? V / 10.0 : V / 20.0
                    );
                }
            }
        }
    }

    return m;
}

void _compute_grad_coeff(double *b, double *c, double *d, struct Mesh mesh, int ti) {
    struct Node *n = mesh.nodes;
    struct Tetra t = mesh.tetras[ti];
    int i, j, k, l;
    
    for (i = 0; i < 4; ++i) {
        // Solve for \nabla \phi_i = (b_i, c_i, d_i) with Gaussian elimination
        double M[4][5];
        for (j = 0; j < 4; ++j) {
            int vi = t.v[j];
            M[j][0] = 1.0;
            M[j][1] = n[vi].x;
            M[j][2] = n[vi].y;
            M[j][3] = n[vi].z;
            M[j][4] = (j == i) ? 1.0 : 0.0;
        }

        for (j = 0; j < 4; ++j) {
            int pivot_row = j;
            double max_val = fabs(M[j][j]);

            for (k = j + 1; k < 4; ++k) {
                if (fabs(M[k][j]) > max_val) {
                    max_val = fabs(M[k][j]);
                    pivot_row = k;
                }
            }

            if (pivot_row != j) {
                for (k = 0; k < 5; ++k) {
                    double tmp = M[j][k];
                    M[j][k] = M[pivot_row][k];
                    M[pivot_row][k] = tmp;
                }
            }

            double pivot = M[j][j];
            assert(fabs(pivot) > 1e-14);

            for (k = j; k < 5; ++k) {
                M[j][k] /= pivot;
            }

            for (k = 0; k < 4; ++k) {
                if (k == j) continue;
                double factor = M[k][j];
                for (l = j; l < 5; ++l) {
                    M[k][l] -= factor * M[j][l];
                }
            }
        }

        b[i] = M[1][4];
        c[i] = M[2][4];
        d[i] = M[3][4];
    }
}

struct CSR_Matrix hamiltonian(struct Mesh mesh) {
    // Creates the hamiltonian matrix
    // H_jk = a(center) \int dV \nabla\phi_j \cdot \nabla\phi_k
    //      + V(center) \int dV \phi_j \phi_k
    // For adjacent j and k,
    // \int_tetra dV \nabla\phi_j \cdot \nabla\phi_k
    // = (b_j b_k + c_j c_k + d_j d_k) / 36V
    // where b_i = y_j (z_k - z_l) + y_k (z_l - z_j) + y_l (z_j - z_k)
    //       c_i = z_j (x_k - x_l) + z_k (x_l - x_j) + z_l (x_j - x_k)
    //       c_i = x_j (y_k - y_l) + x_k (y_l - y_j) + x_l (y_j - y_k)
    
    struct CSR_Matrix h = create_matrix(mesh.N_node);

    int i;
    for (i = 0; i < mesh.N_tetra; ++i) {
        double V = volume(mesh.tetras[i], mesh.nodes);
        int j, k;

        struct Node cen = center(mesh.tetras[i], mesh.nodes);
        double b[4], c[4], d[4];
        _compute_grad_coeff(b, c, d, mesh, i);

        for (j = 0; j < 4; ++j) {
            for (k = 0; k < 4; ++k) {
                if (!mesh.nodes[mesh.tetras[i].v[j]].is_bdry &&
                    !mesh.nodes[mesh.tetras[i].v[k]].is_bdry) {
                    add_entry(
                        &h,
                        mesh.node_to_idx[mesh.tetras[i].v[j]],
                        mesh.node_to_idx[mesh.tetras[i].v[k]],
                        coeff_a(cen) * (b[j]*b[k] + c[j]*c[k] + d[j]*d[k]) * V +
                        potential(cen) * ((j == k) ? V * 0.6 : V * 0.3)
                    );

                }
            }
        }
    }

    return h;
}

// =============================

int cmp_entry(const void *a, const void *b) {
    const struct CSR_Entry *ea = a;
    const struct CSR_Entry *eb = b;
    return ea->col - eb-> col;
}

void sort_row(struct CSR_Row *r) {
    qsort(r->entries, r->cnt, sizeof(struct CSR_Entry), cmp_entry);

    int i, j = 0;
    for (i = 1; i < r->cnt; ++i) {
        if (r->entries[i].col == r->entries[j].col) {
            r->entries[j].val += r->entries[i].val;
        } else {
            r->entries[++j] = r->entries[i];
        }
    }
    r->cnt = j + 1;
}

void sort_matrix(struct CSR_Matrix *m) {
    int i;
    for (i = 0; i < m->N_node; ++i) {
        sort_row(&m->rows[i]);
    }
}

void fwrite_matrix(const char *filename, struct CSR_Matrix m) {
    // Assuming m has been sorted!
    int *row_ptr = malloc((m.N_node + 1) * sizeof(int));
    int i, total = 0;
    row_ptr[0] = 0;
    for (i = 0; i < m.N_node; ++i) {
        total += m.rows[i].cnt;
        row_ptr[i+1] = total;
    }

    int *col_ind = malloc(total * sizeof(int));
    double *values = malloc(total * sizeof(double));

    int j, idx = 0;
    for (i = 0; i < m.N_node; ++i) {
        for (j = 0; j < m.rows[i].cnt; ++j) {
            col_ind[idx] = m.rows[i].entries[j].col;
            values[idx] = m.rows[i].entries[j].val;
            idx += 1;
        }
    }

    FILE *f = fopen(filename, "w");
    fprintf(f, "%d %d\n", m.N_node, total);

    for (i = 0; i < m.N_node + 1; ++i)
        fprintf(f, "%d ", row_ptr[i]);
    fprintf(f, "\n");

    for (i = 0; i < total; ++i)
        fprintf(f, "%d ", col_ind[i]);
    fprintf(f, "\n");

    for (i = 0; i < total; ++i)
        fprintf(f, "%.16e ", values[i]);
    
    fclose(f);

    free(row_ptr);
    free(col_ind);
    free(values);
}

void print_matrix(struct CSR_Matrix m) {
    printf("Matrix with %d nodes:\n", m.N_node);
    int i;
    for (i = 0; i < m.N_node; ++i) {
        printf("Row %d:\n", i);
        int j;
        for (j = 0; j < m.rows[i].cnt; ++j) {
            printf("Col %d -> Val %6.4f\n", m.rows[i].entries[j].col, m.rows[i].entries[j].val);
        }
    }
}

// =============================

int main() {
    struct Mesh mesh = read_mesh("nodes.1.node", "nodes.1.ele");
    printf("Read %d internal and %d boundary nodes\n", mesh.N_node, mesh.N_bdry);
    printf("First node: (%6.4f, %6.4f, %6.4f) (%s)\n",
            mesh.nodes[0].x, mesh.nodes[0].y, mesh.nodes[0].z,
            (mesh.nodes[0].is_bdry) ? "Boundary" : "Interior"
          );
    printf("Read %d tetrahedra\n", mesh.N_tetra);
    printf("First tetrahedron: (%d, %d, %d, %d)\n",
            mesh.tetras[0].v[0] + 1,
            mesh.tetras[0].v[1] + 1,
            mesh.tetras[0].v[2] + 1,
            mesh.tetras[0].v[3] + 1
    );
    printf("Volume of first tetrahedron: %6.4f\n", volume(mesh.tetras[0], mesh.nodes));
    
    struct Node cen = center(mesh.tetras[0], mesh.nodes);
    printf("Center of first tetrahedron: (%6.4f, %6.4f, %6.4f)\n", cen.x, cen.y, cen.z);

    double b[4], c[4], d[4];
    _compute_grad_coeff(b, c, d, mesh, 0);
    printf("Gradient components:\n");
    int i;
    for (i = 0; i < 4; ++i) {
        printf("\\nabla\\phi_%d = (%6.4f, %6.4f, %6.4f)\n", i, b[i], c[i], d[i]);
    }

    struct CSR_Matrix m = mass(mesh);
    sort_matrix(&m);
    // print_matrix(m);
    fwrite_matrix("mass.csr", m);

    struct CSR_Matrix h = hamiltonian(mesh);
    sort_matrix(&h);
    // print_matrix(h);
    fwrite_matrix("hamiltonian.csr", h);

    free_mesh(mesh);
    free_matrix(m);
    free_matrix(h);
    return 0;
}
