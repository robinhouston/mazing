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
