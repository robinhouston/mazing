typedef struct {
    mpz_t ov; /* original value */
    mpz_t bv; /* value after running Bareiss algorithm */
} ent_t;

typedef struct {
    int offset;
    ent_t entries[];
} row_t;

/* A symmetric band matrix, optimised for progressive determinant computations */
typedef struct {
    int n; /* Number of allocated rows in matrix */
    int w; /* Number of elements in each row */
    int nr; /* Number of active rows (<= number of allocated rows) */
    int det_start; /* Which element to start computing a sub-determinant */
    int min_changed; /* Min index of changed element; n if nothing changed */
    ent_t zero; /* Always zero: used for out-of-band entries */
    row_t *rows[];
} matrix_t;

typedef char direction;
#define DIR_N 0x1
#define DIR_E 0x2
#define DIR_S 0x4
#define DIR_W 0x8

typedef struct {
    int width;
    int height;
    direction conn[]; /* has (width * height) elements */
} maze_t;

maze_t *maze_by_index(int width, int height, mpz_t index);
void maze_free(maze_t *maze);
void maze_print(maze_t *maze);
