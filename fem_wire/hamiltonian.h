struct Node {
    double x, y, z;
    int is_bdry;
    int match;
};

struct Tetra {
    int v[4];
};

struct Mesh {
    struct Node *nodes;
    struct Tetra *tetras;

    int N_tot;
    int N_bdry;
    int N_node;
    int N_tetra;

    int *node_to_idx;
};

struct CSR_Entry {
    int col;
    double val;
};

struct CSR_Row {
    int cnt;
    int capacity;
    struct CSR_Entry *entries;
};

struct CSR_Matrix {
    int N_node;
    struct CSR_Row *rows;
};

double COEFF_A = 1.0;

double potential(struct Node p);

struct CSR_Matrix create_matrix(int N_node);
void add_entry(struct CSR_Matrix *m, int row, int col, double val);
void free_matrix(struct CSR_Matrix m);

struct Node *read_nodes(const char *filename, int *N_node, int *N_bdry);
struct Tetra *read_tetras(const char *filename, int *N_tetra, int N_tot);
struct Mesh read_mesh(const char *node_file, const char *tetra_file);
void free_mesh(struct Mesh mesh);

double volume(struct Tetra t, struct Node *nodes);
struct Node center(struct Tetra t, struct Node *nodes);

struct CSR_Matrix shape(struct Mesh mesh);
void _compute_grad_coeff(double *b, double *c, double *d, struct Mesh mesh, int i);
struct CSR_Matrix hamiltonian(struct Mesh mesh);

int cmp_entry(const void *a, const void *b);
void sort_row(struct CSR_Row *r);
void sort_matrix(struct CSR_Matrix *m);
void fwrite_matrix(const char *filename, struct CSR_Matrix m);
void print_matrix(struct CSR_Matrix m);
