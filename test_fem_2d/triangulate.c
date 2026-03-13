#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG 0

struct Point {
    double x;
    double y;
};

struct Edge {
    int v[2];
};

struct Triangle {
    int v[3];
};

struct Point *read_nodes(const char *filename, int *N_node) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("file open");
        exit(1);
    }

    int N_bdry;
    fscanf(f, "%d %d", N_node, &N_bdry);
    *N_node += N_bdry;

    struct Point *nodes = malloc(*N_node * sizeof(struct Point));
    int i;
    for (i = 0; i < *N_node; ++i)
        fscanf(f, "%lf %lf", &nodes[i].x, &nodes[i].y);

    fclose(f);
    return nodes;
}

void fprint_triangles(const char *filename, struct Triangle *t, int N_tri) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("file open");
        exit(1);
    }

    fprintf(f, "%d\n", N_tri);
    
    int i;
    for (i = 0; i < N_tri; ++i)
        fprintf(f, "%d %d %d\n", t[i].v[0] + 1, t[i].v[1] + 1, t[i].v[2] + 1);

    fclose(f);
}

struct Point *supertriangle(struct Point *nodes, int N_node) {
    double minX = nodes[0].x, maxX = nodes[0].x,
           minY = nodes[0].y, maxY = nodes[0].y;
    int i;
    for (i = 1; i < N_node; ++i) {
        minX = fmin(minX, nodes[i].x);
        maxX = fmax(maxX, nodes[i].x);
        minY = fmin(minY, nodes[i].y);
        maxY = fmax(maxY, nodes[i].y);
    }

    double dx = maxX - minX, dy = maxY - minY;
    struct Point *a = malloc(3 * sizeof(struct Point));
    a[0] = (struct Point) {minX - dx,     minY - dy    };
    a[1] = (struct Point) {minX - dx,     maxY + dy * 3};
    a[2] = (struct Point) {maxX + dx * 3, minY - dy    };

    return a;
}

int in_circumcircle(
    struct Point p,
    struct Triangle t,
    struct Point *nodes,
    struct Point *super
) {
    struct Point a, b, c;
    if (t.v[0] >= 0) a = nodes[t.v[0]]; else a = super[-t.v[0] - 1];
    if (t.v[1] >= 0) b = nodes[t.v[1]]; else b = super[-t.v[1] - 1];
    if (t.v[2] >= 0) c = nodes[t.v[2]]; else c = super[-t.v[2] - 1];
    
    double orient =
        (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
 
    double ax = a.x - p.x,
           ay = a.y - p.y,
           bx = b.x - p.x,
           by = b.y - p.y,
           cx = c.x - p.x,
           cy = c.y - p.y;
    
    double det =
        (ax * ax + ay * ay) * (bx * cy - cx * by) -
        (bx * bx + by * by) * (ax * cy - cx * ay) +
        (cx * cx + cy * cy) * (ax * by - bx * ay);

    if (DEBUG) printf("%6.2f %6.2f\n", orient, det);
    if (orient > 0) return det > 0;
    return det < 0;
}

void remove_bad_triangles(struct Triangle *t, int *bad, int *N_tri, int N_bad) {
    int i, j, write, remove;

    write = 0;
    for (i = 0; i < *N_tri; ++i) {
        remove = 0;
        for (j = 0; j < N_bad; ++j) {
            if (i == bad[j]) {
                remove = 1;
                break;
            }
        }

        if (!remove)
            t[write++] = t[i];
    }
    *N_tri = write;
}

void add_point(
    struct Triangle *t,
    struct Point *super,
    struct Point *nodes,
    int to_add,
    int N_node,
    int *N_tri
) {
    int bad[*N_tri];
    int i, N_bad = 0;
    for (i = 0; i < *N_tri; ++i) {
        if (in_circumcircle(nodes[to_add], t[i], nodes, super)) {
            bad[N_bad++] = i;
            if (DEBUG)
                printf("Node %d in triangle %d: {%d, %d, %d}\n", to_add, i, t[i].v[0], t[i].v[1], t[i].v[2]);
        }
    }

    struct Edge *polygon = malloc((10 * *N_tri) * sizeof(struct Edge));
    int N_poly = 0;
    for (i = 0; i < N_bad; ++i) {
        int j;
        for (j = 0; j < 3; ++j) {
            struct Edge edge = {t[bad[i]].v[j], t[bad[i]].v[(j+1)%3]};
            int unique = 1;
            int k;
            for (k = 0; k < N_bad; ++k) {
                if (i == k) continue;

                int l;
                for (l = 0; l < 3; ++l) {
                    if ((edge.v[0] == t[bad[k]].v[l] && edge.v[1] == t[bad[k]].v[(l+1)%3]) ||
                        (edge.v[1] == t[bad[k]].v[l] && edge.v[0] == t[bad[k]].v[(l+1)%3])) {
                        unique = 0;
                        break;
                    }
                }
                if (!unique) break;
            }
            if (unique)
                polygon[N_poly++] = edge;
        }
    }
    if (DEBUG) printf("Starting removal...\n");
    remove_bad_triangles(t, bad, N_tri, N_bad);

    if (DEBUG) {
        printf("After removal:\n");
        for (i = 0; i < *N_tri; ++i)
            printf("Triangle %d: {%d, %d, %d}\n", i, t[i].v[0], t[i].v[1], t[i].v[2]);
        printf("\n");
    }

    for (i = 0; i < N_poly; ++i) {
        t[*N_tri].v[0] = polygon[i].v[0];
        t[*N_tri].v[1] = polygon[i].v[1];
        t[*N_tri].v[2] = to_add;
        ++*N_tri;
    }

    if (DEBUG) {
        printf("After adding:\n");
        for (i = 0; i < *N_tri; ++i)
            printf("Triangle %d: {%d, %d, %d}\n", i, t[i].v[0], t[i].v[1], t[i].v[2]);
        printf("\n\n");
    }
    free(polygon);
}

void remove_super(struct Triangle *t, int *N_tri) {
    int i;
    for (i = *N_tri - 1; i >= 0; --i) {
        if (t[i].v[0] < 0 || t[i].v[1] < 0 || t[i].v[2] < 0) {
            int j;
            --*N_tri;
            for (j = i; j < *N_tri; ++j)
                t[j] = t[j + 1];
        }
    }
}

struct Triangle *bowyer_watson(struct Point *nodes, int N_node, int *N_tri) {
    struct Triangle *t = malloc((100 * N_node + 100) * sizeof(struct Triangle));
    struct Point *super = supertriangle(nodes, N_node);
    t[0] = (struct Triangle) {{-1, -2, -3}};
    *N_tri = 1;

    int i;
    for (i = 0; i < N_node; ++i) {
        if (DEBUG) {
            printf("Adding node %d:\n", i);
            int j;
            for (j = 0; j < *N_tri; ++j)
                printf("Triangle %d: {%d, %d, %d}\n", j, t[j].v[0], t[j].v[1], t[j].v[2]);
            printf("\n");
        }
        add_point(t, super, nodes, i, N_node, N_tri);
    }

    remove_super(t, N_tri);

    free(super);
    return t;
}

int main(void) {
    int N_node;
    struct Point *nodes = read_nodes("nodes.txt", &N_node);
    printf("Read %d nodes\n", N_node);

    int N_tri;
    struct Triangle *t = bowyer_watson(nodes, N_node, &N_tri);
    printf("Created %d triangles\n", N_tri);

    fprint_triangles("triangles.txt", t, N_tri);

    free(nodes);
    free(t);
    return 0;
}
