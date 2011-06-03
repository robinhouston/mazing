/* fmc.c - Fast Maze Counter */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "fmc.h"

/** Small helper functions **/

/* The nth triangular number */
inline static int tri(int n)
{
    return n * (n+1) / 2;
}

/* The largest power of two less than or equal to n */
inline static int msb(int n)
{
    int p = 1;
    while (p <= n) p <<= 1;
    return p >> 1;
}


/** Matrices **/

/* A symmetric matrix of integers.
The (i,j)th entry of m, where i >= j, is stored in m->entries[tri(i) + j] */
typedef struct {
    int n;
    mpz_t entries[];
} fmc_matrix;

/* Allocate an 'n'x'n' zero matrix. */
static fmc_matrix *fmc_matrix_init(int n)
{
    int num_ents = tri(n);
    fmc_matrix *m = malloc(sizeof(fmc_matrix) + sizeof(mpz_t) * num_ents);
    m->n = n;
    for (int i = 0; i < num_ents; i++)
        mpz_init(m->entries[i]);
    
    return m;
}

/* Free the matrix. */
static void fmc_matrix_free(fmc_matrix *m)
{
    int num_ents = tri(m->n);
    for (int i = 0; i < num_ents; i++)
        mpz_clear(m->entries[i]);
    free(m);
}

inline static mpz_t *ent(fmc_matrix *m, int i, int j)
{
    if (i < j) {
        int x = i; i = j; j = x;
    }
    
    return &m->entries[tri(i) + j];
}

/*
static void fmc_matrix_print(FILE *out, char *name, fmc_matrix *m)
{
    fprintf(out, "=== %s ===\n", name);
    for (int i = 0; i < m->n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            gmp_fprintf(out, "%Zd ", *ent(m, i, j));
        }
        fputc('\n', out);
    }
    fputc('\n', out);
}
*/

/* Copy the contents of one matrix over another. */
static void fmc_matrix_set(fmc_matrix *dest, fmc_matrix *src)
{
    int n = tri(src->n);
    for (int i=0; i < n; i++)
        mpz_set(dest->entries[i], src->entries[i]);
}

/* result := result - other */
static void fmc_matrix_sub(fmc_matrix *result, fmc_matrix *other)
{
    int n = tri(result->n);
    for (int i=0; i < n; i++)
        mpz_sub(result->entries[i], result->entries[i], other->entries[i]);
}

/* dest := ma x mb */
static void fmc_matrix_mul(fmc_matrix *dest, fmc_matrix *ma, fmc_matrix *mb)
{
    int n = dest->n;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            mpz_t *cell = &dest->entries[tri(i) + j];
            mpz_set_si(*cell, 0);
            for (int k = 0; k < n; k++)
                mpz_addmul(*cell, *ent(ma,i,k), *ent(mb,k,j));
        }
    }
}

/* Multiply 'src' by the tridiagonal matrix that has 4 down the main diagonal
and -1 immediately above and below, and stores the result in 'dest'.

("mobfi" stands for "minus-one-bordered four-times-identity".)

Assumes dest != src.
*/
static void fmc_matrix_mul_mobfi(fmc_matrix *dest, fmc_matrix *src)
{
    int n = dest->n;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int x = tri(i) + j;
            mpz_t *cell = &dest->entries[x];
            mpz_mul_2exp(*cell, src->entries[x], 2);
            if (i > 0)
                mpz_sub(*cell, *cell, *ent(src, i-1, j));
            if (i < n-1)
                mpz_sub(*cell, *cell, *ent(src, i+1, j));
        }
    }
}

/* Perform the Bareiss algorithm. */
static void bareiss(fmc_matrix *m)
{
    int n = m->n;
    mpz_t *mkk=0, *mkk_prev;
    
    for (int k=0; k < n; k++)
    {
        mkk_prev = mkk;
        mkk = ent(m,k,k);
        
        for (int i=k+1; i < n; i++)
        {
            mpz_t *mik = ent(m,i,k);
            for (int j = k+1; j <= i; j++)
            {
                mpz_t *mij = ent(m,i,j);
                mpz_t *mjk = ent(m,j,k);
                
                mpz_mul(*mij, *mij, *mkk);
                mpz_submul(*mij, *mik, *mjk);
                if (k > 0) mpz_divexact(*mij, *mij, *mkk_prev);
            }
        }
    }
}

/* Swap two matrix pointers */
inline static void swap(fmc_matrix **x, fmc_matrix **y)
{
    fmc_matrix *t = *x; *x = *y; *y = t;
}

/* The determinant of the block matrix that represents
the planar dual of the 'width'x'height' grid. 

Because this block matrix has a very simple structure,
we can compute its determinant very efficiently by
expressing the determinant as a recurrence and then
*/
static fmc_matrix *dmf(int width, int height)
{
    int n = width - 1;
    
    fmc_matrix *a = fmc_matrix_init(n);
    fmc_matrix *b = fmc_matrix_init(n);
    fmc_matrix *c = fmc_matrix_init(n);
    
    /* scratch space */
    fmc_matrix *new_a = fmc_matrix_init(n);
    fmc_matrix *new_b = fmc_matrix_init(n);
    fmc_matrix *temp = fmc_matrix_init(n);
    
    for (int i=0; i < n; i++)
    {
        mpz_set_si(*ent(a,i,i), -1);
        mpz_set_si(*ent(c,i,i), +1);
    }
        
    for (int bit = msb(height); bit > 0; bit >>= 1)
    {
        /* a, b := b^2 - a^2, bc - ab */
        fmc_matrix_mul(new_a, b, b);
        fmc_matrix_mul(temp, a, a);
        fmc_matrix_sub(new_a, temp);
        
        fmc_matrix_mul(new_b, b, c);
        fmc_matrix_mul(temp, a, b);
        fmc_matrix_sub(new_b, temp);
        
        swap(&a, &new_a);
        swap(&b, &new_b);
        
        if ((height & bit) > 0)
        {
            /* a, b := b, bM - a */
            fmc_matrix_set(new_a, b);
            fmc_matrix_mul_mobfi(new_b, b);
            fmc_matrix_sub(new_b, a);
            swap(&a, &new_a);
            swap(&b, &new_b);
        }
        
        /* c := bM - a */
        fmc_matrix_mul_mobfi(c, b);
        fmc_matrix_sub(c, a);
    }
    fmc_matrix_free(a);
    fmc_matrix_free(c);
    fmc_matrix_free(new_a);
    fmc_matrix_free(new_b);
    fmc_matrix_free(temp);
    
    return b;
}

/* Fast Maze Counter.

Count the number of mazes on a 'width'x'height grid,
and store the result in 'out'.
*/
void fmc(mpz_t *out, int width, int height)
{
    fmc_matrix *m = dmf(width, height);
    bareiss(m);
    mpz_set(*out, m->entries[tri(m->n) - 1]);
    fmc_matrix_free(m);
}

