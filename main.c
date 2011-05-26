#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <gmp.h>

#include "mazing.h"
#include "fmc.h"

void print_count(int width, int height)
{
    mpz_t count;
    
    mpz_init(count);
    fmc(&count, width, height);
    
    gmp_printf("There are %Zd different mazes on a %dx%d grid.\n",
        count, width, height);
    
    mpz_clear(count);
}

void print_maze(int width, int height, mpz_t index)
{
    maze_t *maze = maze_by_index(width, height, index);
    
    // XXXX
    
    maze_free(maze);
}

int main(int argc, char **argv)
{
    int width, height;
    mpz_t index;
    
    if (argc < 3 || argc > 4)
    {
        fprintf(stderr, "Usage: %s width height [index]\n", argv[0]);
        return EX_USAGE;
    }
    
    width = atoi(argv[1]);
    height = atoi(argv[2]);
    
    if (width <= 0 || height <= 0)
    {
        fprintf(stderr, "Usage: %s width height [index]\n", argv[0]);
        fprintf(stderr, "width and height must be positive\n");
        return EX_USAGE;
    }
    
    if (argc == 3)
    {
        print_count(width, height);
        return 0;
    }
    
    /* Construct a maze by index */
    mpz_init(index);
    gmp_sscanf(argv[3], "%Zd", &index);
    print_maze(width, height, index);
    mpz_clear(index);
    return 0;
}
