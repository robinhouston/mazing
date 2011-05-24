#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <gmp.h>

#include "mazing.h"

int main(int argc, char **argv)
{
    int width, height;
    matrix_t *m;
    
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s width height\n", argv[0]);
        return EX_USAGE;
    }
    
    width = atoi(argv[1]);
    height = atoi(argv[2]);
    
    m = maze_matrix(width, height);
    matrix_bareiss(m);
    
    matrix_print_det(stdout, m);
    printf("\n");
    
    matrix_free(m);
    return 0;
}