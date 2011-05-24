typedef struct {
    int offset;
    mpz_t entries[];
} row_t;

/* An economical representation of a sparse-ish symmetric matrix */
typedef struct {
    int n; /* Number of rows in matrix */
    int w; /* Number of elements in each row */
    mpz_t zero; /* Always zero */
    row_t *rows[];
} matrix_t;

matrix_t *maze_matrix(int width, int height);
void matrix_free(matrix_t *m);

void matrix_bareiss(matrix_t *m);
void matrix_print_det(FILE *out, matrix_t *m);
