#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

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

/** Matrix functions **

matrix_t represents a symmetric band matrix with (big) integer entries.

I was going to explain that in detail, but then I found the wikipedia
page http://en.wikipedia.org/wiki/Band_matrix which reveals that I am
not in fact the first person to have thought of representing matrices
this way.

Each row contains the half-band ending with the diagonal entry.
The 'w' member is the number of entries in a typical row, i.e.
1 + the half-bandwidth. The row also has an 'offset', which is
the column number of its first entry. This is redundant, but
makes the code simpler.
*/

/* Allocate a zero matrix with 'num_rows' rows (and columns),
and a row length of 'row_length', i.e. a half-bandwidth of
(row_length - 1). */
matrix_t *matrix_init(int num_rows, int row_length)
{
    matrix_t *m = malloc(sizeof(matrix_t) + sizeof(row_t *) * num_rows);
    
    m->n = num_rows;
    m->w = row_length;
    mpz_init(m->zero);
    
    for (int i = 0; i < m->n; i++)
    {
        int this_row_len = min(i+1, m->w);
        row_t *row = m->rows[i] = malloc(sizeof(row_t) + sizeof(mpz_t) * this_row_len);
        
        for (int j=0; j < this_row_len; j++)
        {
            mpz_init(row->entries[j]);
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
            mpz_clear(row->entries[j]);
        }
        
        free(row);
    }
    
    mpz_clear(m->zero);
    free(m);
}

/* The Laplacian matrix for a 'width' x 'height' grid. */
matrix_t *grid_matrix(int width, int height)
{
    int n = width * height;
    matrix_t *m = matrix_init(n, width + 1);
    
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
        
        row->offset = first_row ? 0 : i - width;
        
        mpz_set_si(row->entries[i - row->offset], num_neighbours);
        if (!first_row) mpz_set_si(row->entries[i - width - row->offset], -1);
        if (!first_col) mpz_set_si(row->entries[i - 1 - row->offset], -1);
    }
    
    return m;
}

/* Get a pointer to the specified entry of the row */
inline static mpz_t *ent_r(row_t *row, int j)
{
    return &row->entries[j - row->offset];
}
/* Get a pointer to the (i,j)th entry, assuming i >= j */
inline static mpz_t *ent(matrix_t *m, int i, int j)
{
    row_t *row = m->rows[i];
    if (j < row->offset)
        return &m->zero;
    return ent_r(row, j);
}
/* Get a pointer to the (i,j)th entry */
inline static mpz_t *ent_eo(matrix_t *m, int i, int j)
{
    return (i < j) ? ent(m,j,i) : ent(m,i,j);
}

/* Print the matrix to stdout: useful for debugging */
static void matrix_print(matrix_t *m, char *name)
{
    printf("=== %s ===\n", name);
    for (int i = 0; i < m->n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            gmp_printf("%2Zd ", *ent(m, i, j));
        }
        fputc('\n', stdout);
    }
    fputc('\n', stdout);
}

/** Maze finding **

Find a maze by index, descending the tree guided by index
and weighing the branches using what is essentially
the exact determinant algorithm of Bareiss, but avoiding
unnecessary calculations with two tricks:

1. The algorithm itself takes advantage of the band matrix
   structure, so that instead of needing O(n^3) operations,
   as the ordinary Bareiss algorithm does, it needs only
   O(n * w^2) where w is the band width. (The symmetry also
   saves a constant factor of 2^3, of course.)

2. Since we need determinants of similar matrices in succession,
   it is wasteful to recalculate the whole Bareiss matrix every time.
   Instead we just recompute the part that has changed.

*/
void matrix_count_trees(matrix_t *m, int start_node, int end_node)
{
    int n = m->n, w = m->w;
    mpz_t *mkk=0, *mkk_prev=0;
    
    for (int k=start_node; k < end_node - 1; k++)
    {
        mkk_prev = mkk;
        mkk = ent(m,k,k);
        for (int i = k+1; i < min(n, k+w); i++)
        {
            row_t *row_i = m->rows[i];
            mpz_t *mik = ent(m,i,k);
            for (int j = max(k+1, row_i->offset); j <= i; j++)
            {
                mpz_t *mjk = ent(m,j,k);
                mpz_t *mij = ent_r(row_i,j);
                
                mpz_mul(*mij, *mij, *mkk);
                mpz_submul(*mij, *mik, *mjk);
                if (mkk_prev) mpz_divexact(*mij, *mij, *mkk_prev);
            }
        }
        if (k+w < n) {
            int i = k+w;
            row_t *row_i = m->rows[i];
            for (int j = max(k+1, row_i->offset); j <= i; j++)
            {
                mpz_t *mij = ent_r(row_i,j);
                mpz_mul(*mij, *mij, *mkk);
            }
        }
    }
}

maze_t *maze_init(int width, int height)
{
    int n = width * height;
    maze_t *maze = malloc(sizeof(maze_t) + sizeof(direction) * n);
    
    maze->width = width;
    maze->height = height;
    memset(maze->conn, 0, n);
    
    return maze;
}

void maze_free(maze_t *maze)
{
    free(maze);
}

static void recompute(matrix_t *m, matrix_t *im, int from, int to)
{
    mpz_t *mkk_prev;
    
    for (int i = from; i < to; i++)
        for (int j = from; j <= i; j++)
            mpz_set(*ent(m,i,j), *ent(im,i,j));
    
    if (from > m->w)
    {
        int k = from - m->w;
        mpz_t *mkk = ent(m,k,k);
        
        for (int i = from; i < to; i++)
            for (int j = from; j <= i; j++)
            {
                mpz_t *mij = ent(m,i,j);
                mpz_mul(*mij, *mij, *mkk);
            }
        
        mkk_prev = mkk;
    }
    else
        mkk_prev = 0;
    
    for (int k = max(1, from - m->w); k < to; k++)
    {
        mpz_t *mkk = ent(m,k,k);
        for (int i = max(from, k+1); i < to; i++)
        {
            mpz_t *mik = ent(m,i,k);
            for (int j = max(from, k+1); j <= i; j++)
            {
                mpz_t *mij = ent(m,i,j);
                mpz_t *mjk = ent(m,j,k);
                
                mpz_mul(*mij, *mij, *mkk);
                mpz_submul(*mij, *mik, *mjk);
                if (mkk_prev) mpz_divexact(*mij, *mij, *mkk_prev);
            }
        }
        mkk_prev = mkk;
    }
}

static bool try_edge(matrix_t *im, matrix_t *m, mpz_t *index, int *node_chain, int last_node,
    int from_cell, int to_cell)
{
    int n_i = chain_root(node_chain, to_cell);
    int n_j = chain_root(node_chain, from_cell);
    
    if (n_i < n_j) {
        /* Swap them round so that n_i >= n_j */
        swap_ints(&n_i, &n_j);
        swap_ints(&to_cell, &from_cell);
    }
    
    mpz_t *im_ii = ent(im, n_i, n_i);
    mpz_t *im_ij = ent(im, n_i, n_j);
    
    if (mpz_cmp_ui(*im_ij, 0) >= 0)
    {
        printf("There aren't any edges %d->%d (%d->%d)\n", from_cell, to_cell, n_j, n_i);
        return false;
    }
    
    mpz_t *im_jj = ent(im, n_j, n_j);
    mpz_t *count_wo_edge = ent(m,last_node-1,last_node-1);
    
    //gmp_printf("Testing edge %d -> %d\n", from_cell, to_cell);
    
    /* How many mazes are there without this edge? */
    mpz_sub_ui(*im_ii, *im_ii, 1);
    mpz_sub_ui(*im_jj, *im_jj, 1);
    mpz_add_ui(*im_ij, *im_ij, 1);
    
    recompute(m, im, min(n_i, n_j), last_node);
    gmp_printf("without (%d, %d) there are %Zd (index = %Zd)\n", from_cell, to_cell, *count_wo_edge, *index);
    
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
        
        matrix_print(im, "before munging");
        gmp_printf("im[%d][%d] = %Zd\n", n_i, n_j, *im_ij);
        mpz_add(*im_jj, *im_jj, *im_ii);
        mpz_add(*im_jj, *im_jj, *im_ij);
        printf("Zapping row %d\n", n_i);
        for (int k = start_node; k < end_node; k++)
        {
            if (k != n_i)
                mpz_add(*ent_eo(im,n_j,k), *ent_eo(im,n_j,k), *ent_eo(im,n_i,k));
            mpz_set_ui(*ent_eo(im,n_i,k), (k==n_i ? 1 : 0));
        }
        mpz_sub(*index, *index, *count_wo_edge);
        printf("node_chain[%d] := %d\n", to_cell, n_j);
        node_chain[to_cell] = n_j;
        return true;
    }
}

void maze_print(maze_t *maze)
{
    int w = maze->width, h = maze->height;
    
    for (int y = 0; y < h; y++)
    {
        for (int x = 0; x < w; x++)
            printf( (maze->conn[w*y + x] & DIR_N) == 0 ? "+--+" : "+  +" );
        printf("\n");
        
        for (int x = 0; x < w; x++)
        {
            printf( (maze->conn[w*y + x] & DIR_W) == 0 ? "+" : " " );
            printf("  ");
            printf( (maze->conn[w*y + x] & DIR_E) == 0 ? "+" : " " );
        }
        printf("\n");
    }
    
    for (int x = 0; x < w; x++)
        printf("+--+");
    
    printf("\n\n");
}

maze_t *maze_by_index(int width, int height, mpz_t index_in)
{
    mpz_t index;
    matrix_t *im = grid_matrix(width, height);
    matrix_t *m = grid_matrix(width, height);
    int n = m->n;
    maze_t *maze = maze_init(width, height);
    int *node_chain = malloc(sizeof(int) * n);
    
    mpz_init_set(index, index_in);
    
    matrix_count_trees(m, 1, n);
    for (int i=0; i < n; i++)
        node_chain[i] = i; /* initially cells and nodes are the same */
    
    for (int i = n - 1; i > 0; i--)
    {
        if (i >= width)
        {
            /* Not on the top row */
            if (try_edge(im, m, &index, node_chain, i+1, i - width, i))
            {
                maze->conn[i-width] |= DIR_S;
                maze->conn[i] |= DIR_N;
                printf("...including edge %d->%d\n", i - width, i);
            }
            matrix_print(im, "im");
            maze_print(maze);
        }
        
        if (i % width)
        {
            /* Not in the leftmost column */
            if (try_edge(im, m, &index, node_chain, i+1, i - 1, i))
            {
                maze->conn[i-1] |= DIR_E;
                maze->conn[i] |= DIR_W;
                printf("...including edge %d->%d\n", i - 1, i);
            }
            matrix_print(im, "im");
            maze_print(maze);
        }
    }
    if (mpz_cmp_ui(index, 0) != 0)
    {
        gmp_printf("!! index == %Zd\n", index);
    }
    
    matrix_free(m);
    matrix_free(im);
    mpz_clear(index);
    free(node_chain);
    return maze;
}


matrix_t *maze_matrix(int width, int height)
{
    matrix_t *m = grid_matrix(width, height);
    matrix_count_trees(m, 0, m->n);
    return m;
}
