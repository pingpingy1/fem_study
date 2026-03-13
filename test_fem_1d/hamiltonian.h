double potential(double x);
double coeff_a(double x);

double **mass(double *x, int N);
double **hamiltonian(double *x, int N);

double **alloc_matrix(int N);
void fprint_matrix(FILE *f, double **m, int N);
void free_matrix(double **m, int N);
double *read_nodes(const char *filename, int *N);
