#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include "mazing.h"


#define DEBUG 1

#if DEBUG
#  define LOG(x) printf x
#else
#  define LOG(x)
#endif


matrix_t *matrix_init(int num_rows, int row_length)
{
    matrix_t *m = malloc(sizeof(matrix_t) + sizeof(row_t *) * num_rows);
    
    m->n = num_rows;
    m->w = row_length;
    mpz_init(m->zero);
    
    for (int i = 0; i < m->n; i++)
    {
        row_t *row = m->rows[i] = malloc(sizeof(row_t) + sizeof(mpz_t) * row_length);
        
        for (int j=0; j < m->w; j++)
        {
            mpz_init(row->entries[j]);
        }
    }
    
    return m;
}

void matrix_free(matrix_t *m)
{
    for (int i=0; i < m->n; i++)
    {
        row_t *row = m->rows[i];
        
        for (int j=0; j < m->w; j++)
        {
            mpz_clear(row->entries[j]);
        }
        
        free(row);
    }
    
    mpz_clear(m->zero);
    free(m);
}

matrix_t *maze_matrix(int width, int height)
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

inline static mpz_t *ent(matrix_t *m, int i, int j)
{
    row_t *row = m->rows[i];
    if (j < row->offset)
        return &m->zero;
    return &row->entries[j - row->offset];
}

void matrix_bareiss(matrix_t *m)
{
    int n = m->n;
    mpz_t *mkk=0, *mkk_prev, *mik, *mjk, *mij;
    
    for (int k=0; k < n - 1; k++)
    {
        mkk_prev = mkk;
        mkk = ent(m,k,k);
        for (int i = k+1; i < n; i++)
        {
            mik = ent(m,i,k);
            for (int j = k+1; j <= i; j++)
            {
                mjk = ent(m,j,k);
                mij = ent(m,i,j);
                
                mpz_mul(*mij, *mij, *mkk);
                mpz_submul(*mij, *mik, *mjk);
                if (k > 0) mpz_divexact(*mij, *mij, *mkk_prev);
            }
        }
    }
}

void matrix_print_det(FILE *out, matrix_t *m)
{
    int last = m->n - 2;
    mpz_t *x = ent(m,last,last);
    mpz_out_str(out, 10, *x);
}
