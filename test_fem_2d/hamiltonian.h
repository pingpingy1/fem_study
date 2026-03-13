struct Point {
    double x;
    double y;
};

struct Triangle {
    int v[3];
};

double potential(struct Point x);
double coeff_a(struct Point x);

double **alloc_matrix(int N);
void fprint_matrix(FILE *f, double **m, int N);
void free_matrix(double **m, int N);

struct Point *read_nodes(const char *filename, int *N_node, int *N_bdry);
struct Triangle *read_triangles(const char *filename, int *N_tri);
void fprint_triangles(const char *filename, struct Triangle *t, int N_tri, int N_node);

double area(struct Triangle t, struct Point *nodes);
struct Point center(struct Triangle t, struct Point *nodes);

double **mass(struct Triangle *t, struct Point *nodes, int N_tri, int N_node);
double **hamiltonian(struct Triangle *t, struct Point *nodes, int N_tri, int N_node);
