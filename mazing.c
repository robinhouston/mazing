#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include <gmp.h>
#include "mazing.h"


/** Small macro-like utility functions **/

inline static int min(int x, int y) { return x < y ? x : y; }
inline static int max(int x, int y) { return x < y ? y : x; }
inline static void swap_ints(int *x, int *y)
{
    int t = *x; *x = *y; *y = t;
}


/** Chain functions **

Is there a standard name for this data structure? I’m calling it
a chain. It represents an equivalence relation on an initial segment
of the natural numbers. The 'chain' itself is just a big-enough array
of ints.
*/

/* Create a discrete chain of length 'n' */
inline static int *chain_init(int n)
{
    int *chain = malloc(sizeof(int) * n);
    for (int i=0; i < n; i++)
        chain[i] = i; /* initially cells and nodes are the same */
    return chain;
}

/* Free a chain. */
inline static void chain_free(int *chain) { free(chain); }

/* The minimum element of the equivalence class
   containing 'index'. */
static int chain_root(int *chain, int index)
{
    if (chain[index] != index)
        return chain[index] = chain_root(chain, chain[index]);
    return index;
}

/* Force equivalence between a and b,
   i.e. combine their equivalence classes. */
static void chain_link(int *chain, int a, int b)
{
    int x = chain_root(chain, a);
    int y = chain_root(chain, b);
    if (y < x) swap_ints(&x, &y);
    if (x == y) return;
    chain[y] = x;
}


/** Maze representation (for output) **/

/* Create a maze of 'width'x'height cells */
maze_t *maze_init(int width, int height)
{
    int n = width * height;
    maze_t *maze = malloc(sizeof(maze_t) + sizeof(direction) * n);
    
    maze->width = width;
    maze->height = height;
    memset(maze->conn, 0, n);
    
    return maze;
}

/* Free the maze */
void maze_free(maze_t *maze)
{
    free(maze);
}

/* Print the maze to stdout, in ascii art style */
void maze_print(maze_t *maze)
{
    int w = maze->width, h = maze->height;
    
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
            printf( (maze->conn[w*y + x] & DIR_N) == 0 ? "+---" : "+   " );
        printf("+\n|");
        
        for (int x = 0; x < w; x++)
        {
            printf("   ");
            printf( (maze->conn[w*y + x] & DIR_E) == 0 ? "|" : " " );
        }
        printf("\n");
    }
    
    for (int x = 0; x < w; x++)
        printf("+---");
    printf("+");
    
    printf("\n\n");
}

/** Matrix functions **

matrix_t represents a symmetric band matrix with (big) integer entries.
Actually each cell contains two integers: the actual value in that cell
is in the 'ov' member of ent_t, and the 'bv' member is used for the
determinant computations described in the "Determinant computation" section.

Each row contains the half-band ending with the row’s diagonal entry.
The 'w' member is the number of entries in a typical row, i.e.
1 + the half-bandwidth. The row also has an 'offset', which is
the column number of its first entry. This is redundant, but
makes the code simpler.
*/

/* Allocate a zero matrix with 'num_rows' rows (and columns),
and a row length of 'row_length', i.e. a half-bandwidth of
(row_length - 1). */
matrix_t *matrix_init(int num_rows, int row_length, int det_start)
{
    matrix_t *m = malloc(sizeof(matrix_t) + sizeof(row_t *) * num_rows);
    
    m->n = m->nr = num_rows;
    m->w = row_length;
    m->det_start = det_start;
    mpz_init(m->zero.ov);
    mpz_init(m->zero.bv);
    
    for (int i = 0; i < num_rows; i++)
    {
        int this_row_len = min(i+1, row_length);
        row_t *row = m->rows[i] = malloc(sizeof(row_t) + sizeof(ent_t) * this_row_len);
        row->offset = i+1 - this_row_len;
        
        for (int j=0; j < this_row_len; j++)
        {
            ent_t *e = &row->entries[j];
            mpz_init(e->ov);
            mpz_init(e->bv);
        }
    }
    
    return m;
}

/* Free the matrix. */
void matrix_free(matrix_t *m)
{
    for (int i=0; i < m->n; i++)
    {
        int this_row_len = min(i+1, m->w);
        row_t *row = m->rows[i];
        
        for (int j=0; j < this_row_len; j++)
        {
            mpz_clear(row->entries[j].ov);
            mpz_clear(row->entries[j].bv);
        }
        
        free(row);
    }
    
    mpz_clear(m->zero.ov);
    mpz_clear(m->zero.bv);
    free(m);
}

/* Get a pointer to the specified entry of the row */
inline static ent_t *ent_r(row_t *row, int j)
{
    return &row->entries[j - row->offset];
}
/* Get a pointer to the (i,j)th entry, assuming i >= j */
inline static ent_t *ent(matrix_t *m, int i, int j)
{
    row_t *row = m->rows[i];
    if (j < row->offset)
        return &m->zero;
    return ent_r(row, j);
}
/* Get a pointer to the (i,j)th entry */
inline static ent_t *ent_eo(matrix_t *m, int i, int j)
{
    return (i < j) ? ent(m,j,i) : ent(m,i,j);
}

/* The Laplacian matrix for a 'width' x 'height' grid. */
matrix_t *grid_matrix(int width, int height)
{
    int n = width * height;
    matrix_t *m = matrix_init(n, width + 1, 1);
    
    for (int i = 0; i < n; i++)
    {
        row_t *row = m->rows[i];
        int r = i / width;
        int c = i % width;
        bool first_row = (r == 0);
        bool last_row = (r == height - 1);
        bool first_col = (c == 0);
        bool last_col = (c == width - 1);
        
        int num_neighbours = (int) !first_row + (int) !last_row
                           + (int) !first_col + (int) !last_col;
        
        mpz_set_si(ent_r(row, i)->ov, num_neighbours);
        if (!first_row) mpz_set_si(row->entries[i - width - row->offset].ov, -1);
        if (!first_col) mpz_set_si(row->entries[i - 1 - row->offset].ov, -1);
    }
    
    return m;
}

/* Print the matrix to stdout: useful for debugging */
static void matrix_print(matrix_t *m, char *name)
{
    printf("=== %s ===\n", name);
    for (int i = 0; i < m->nr; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            ent_t *e = ent(m, i, j);
            gmp_printf("%2Zd,%2Zd ", e->ov, e->bv);
        }
        fputc('\n', stdout);
    }
    fputc('\n', stdout);
}


/** Determinant computation **

We use the exact determinant algorithm of Bareiss, avoiding
unnecessary calculations with two tricks:

1. The algorithm itself takes advantage of the band matrix
   structure, so that instead of needing O(n^3) operations,
   as the ordinary Bareiss algorithm does, it needs only
   O(n * w^2) where w is the band width. (The symmetry also
   saves a constant factor of 2^3, of course.)

2. Since we need determinants of similar matrices in succession,
   it is wasteful to recalculate the whole Bareiss matrix every time:
   instead we just recompute the part that has changed.

*/

/* Compute the Bareiss matrix, initially */
static void det_init(matrix_t *m)
{
    int n = m->nr, w = m->w;
    
    /* Copy the original matrix to the Bareiss matrix */
    for (int i=0; i < n; i++)
    {
        int this_row_len = min(i+1, w);
        row_t *row = m->rows[i];
        
        for (int j=0; j < this_row_len; j++)
        {
            ent_t *e = &row->entries[j];
            mpz_set(e->bv, e->ov);
        }
    }
    
    /* Now run the Bareiss algorithm */
    ent_t *mkk=0, *mkk_prev=0;
    for (int k=m->det_start; k < n - 1; k++)
    {
        mkk_prev = mkk;
        mkk = ent(m,k,k);
        for (int i = k+1; i < min(n, k+w); i++)
        {
            row_t *row_i = m->rows[i];
            ent_t *mik = ent(m,i,k);
            for (int j = max(k+1, row_i->offset); j <= i; j++)
            {
                ent_t *mjk = ent(m,j,k);
                ent_t *mij = ent_r(row_i,j);
                
                mpz_mul(mij->bv, mij->bv, mkk->bv);
                mpz_submul(mij->bv, mik->bv, mjk->bv);
                if (mkk_prev) mpz_divexact(mij->bv, mij->bv, mkk_prev->bv);
            }
        }
        if (k+w < n) {
            int i = k+w;
            row_t *row_i = m->rows[i];
            for (int j = max(k+1, row_i->offset); j <= i; j++)
            {
                ent_t *mij = ent_r(row_i,j);
                mpz_mul(mij->bv, mij->bv, mkk->bv);
            }
        }
    }
    
    /* Nothing has changed since we last recalculated */
    m->min_changed = m->n;
}

inline static void det_changed(matrix_t *m, int i, int j)
{
    int x = min(i,j);
    if (x < m->min_changed) m->min_changed = x;
}

/* Update the Bareiss matrix to account for changes to the underlying matrix */
static void det_update(matrix_t *m)
{
    for (int i = m->min_changed; i < m->nr; i++)
        for (int j = m->min_changed; j <= i; j++) {
            ent_t *e = ent(m,i,j);
            mpz_set(e->bv, e->ov);
        }
    
    ent_t *mkk_prev;
    if (m->min_changed >= m->w)
    {
        int k = m->min_changed - m->w;
        ent_t *mkk = ent(m,k,k);
        
        for (int i = m->min_changed; i < m->nr; i++)
            for (int j = m->min_changed; j <= i; j++)
            {
                mpz_t *mij = &ent(m,i,j)->bv;
                mpz_mul(*mij, *mij, mkk->bv);
            }
        
        mkk_prev = mkk;
    }
    else
        mkk_prev = 0;
    
    for (int k = max(m->det_start, m->min_changed - m->w); k < m->nr; k++)
    {
        ent_t *mkk = ent(m,k,k);
        for (int i = max(m->min_changed, k+1); i < m->nr; i++)
        {
            ent_t *mik = ent(m,i,k);
            for (int j = max(m->min_changed, k+1); j <= i; j++)
            {
                ent_t *mij = ent(m,i,j);
                ent_t *mjk = ent(m,j,k);
                
                mpz_mul(mij->bv, mij->bv, mkk->bv);
                mpz_submul(mij->bv, mik->bv, mjk->bv);
                if (mkk_prev) mpz_divexact(mij->bv, mij->bv, mkk_prev->bv);
            }
        }
        mkk_prev = mkk;
    }
    
    /* Nothing has changed since we last recalculated */
    m->min_changed = m->n;
}

/** Maze finding **

*/

/* Decide which branch of the tree to descend down */
static bool try_edge(matrix_t *m, mpz_t *index, int *node_chain,
    int from_cell, int to_cell)
{
    int n_i = chain_root(node_chain, to_cell);
    int n_j = chain_root(node_chain, from_cell);
    
    if (n_i < n_j) {
        /* Swap them round so that n_i >= n_j */
        swap_ints(&n_i, &n_j);
        swap_ints(&to_cell, &from_cell);
    }
    
    ent_t *m_ii = ent(m, n_i, n_i);
    ent_t *m_ij = ent(m, n_i, n_j);
    
    if (mpz_cmp_ui(m_ij->ov, 0) >= 0)
    {
        // from_cell is already connected to to_cell
        return false;
    }
    
    ent_t *m_jj = ent(m, n_j, n_j);
    mpz_t *count_wo_edge = &ent(m, m->nr - 1, m->nr - 1)->bv;
    
    /* How many mazes are there without this edge? */
    mpz_sub_ui(m_ii->ov, m_ii->ov, 1);
    mpz_sub_ui(m_jj->ov, m_jj->ov, 1);
    mpz_add_ui(m_ij->ov, m_ij->ov, 1);
    det_changed(m, n_j, n_i);
    det_update(m);
    
    if (mpz_cmp(*index, *count_wo_edge) < 0)
    {
        /* Don’t include the edge */
        return false;
    }
    else
    {
        /* Do include it */
        int start_node = max(0, n_j - m->w + 1);
        int end_node = min(m->n, n_i + m->w);
        
        mpz_add(m_jj->ov, m_jj->ov, m_ii->ov);
        mpz_add(m_jj->ov, m_jj->ov, m_ij->ov);
        det_changed(m, n_j, n_i);
        
        for (int k = start_node; k < end_node; k++)
        {
            mpz_t *n_ik = &ent_eo(m,n_i,k)->ov;
            if (k != n_i && mpz_cmp_ui(*n_ik, 0) != 0) {
                mpz_t *n_jk = &ent_eo(m,n_j,k)->ov;
                mpz_add(*n_jk, *n_jk, *n_ik);
                det_changed(m,n_j,k);
            }
            
            int new_value = (k==n_i ? 1 : 0);
            if (mpz_cmp_ui(*n_ik, new_value) != 0) {
                mpz_set_ui(*n_ik, new_value);
                det_changed(m,n_i,k);
            }
        }
        mpz_sub(*index, *index, *count_wo_edge);
        chain_link(node_chain, n_i, n_j);
        return true;
    }
}

maze_t *maze_by_index(int width, int height, mpz_t index_in)
{
    mpz_t index;
    matrix_t *m = grid_matrix(width, height);
    int n = m->n;
    maze_t *maze = maze_init(width, height);
    int *node_chain = chain_init(n);
    
    mpz_init_set(index, index_in);
    
    det_init(m);
    
    for (int i = n - 1; i > 0; i--)
    {
        m->nr = i + 1;
        
        if (i >= width)
        {
            /* Not on the top row */
            if (try_edge(m, &index, node_chain, i - width, i))
            {
                maze->conn[i-width] |= DIR_S;
                maze->conn[i] |= DIR_N;
            }
        }
        
        if (i % width)
        {
            /* Not in the leftmost column */
            if (try_edge(m, &index, node_chain, i - 1, i))
            {
                maze->conn[i-1] |= DIR_E;
                maze->conn[i] |= DIR_W;
            }
        }
    }
    if (mpz_cmp_ui(index, 0) != 0)
    {
        return 0;
    }
    
    matrix_free(m);
    mpz_clear(index);
    chain_free(node_chain);
    return maze;
}
