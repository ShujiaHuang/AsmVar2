/*
 * Encodes nucleotides (A,C,G,T or a,c,g,t) into an integer
 * */
#include <stdio.h>
#include <assert.h>
#include <string.h>

// define the DEBUG symbol
// #define DEBUG 1

inline unsigned int hashEncode(char *seq, unsigned int hashmer) {

    // Encodes nucleotides (A,C,G,T or a,c,g,t) into an integer
    assert(hashmer <= strlen(seq));

    unsigned int h = 0;
    int c, i = 0;
    for (; i < hashmer; ++i) {
        // Just a simple hash function.
        c = seq[i] & 7;  // a,A->1  c,C->3  g,G->7  t,T->4
        if (c == 7) c = 2;

        h = (h << 2) + (unsigned int)(c & 3);
    }

#ifdef DEBUG
    printf("seq = %s; hashmer = %d; h = %d\n", seq, hashmer, h);
#endif

    return h;
}


