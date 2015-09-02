/*
 * Encodes nucleotides (A,C,G,T or a,c,g,t) into an integer
 * */
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <assert.h>

using namespace std;

// define the DEBUG symbol
//#define DEBUG 1
// BKDR Hash Function (https://www.byvoid.com/blog/string-hash-compare)
unsigned int hash_encode(string seq, unsigned int hashmer) {

    // Encodes nucleotides (A,C,G,T or a,c,g,t) into an integer
    assert(hashmer <= seq.length());
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int h = 0;
    int c, i = 0;
    for (; i < hashmer; ++i) {
        h = h * seed + seq[i];
    }

#ifdef DEBUG
    cout << "seq = " << seq << "; hashmer = " << hashmer << "; h = " << h & 0x7FFFFFFF;
#endif

    return (h & 0x7FFFFFFF);
}

unsigned int *hash_array(string seq, unsigned int hashmer) {

    size_t size = seq.length() - hashmer + 1;
    unsigned int *array = new unsigned int[size];
    for (size_t i(0); i < size; ++i) {
        array[i] = hash_encode(seq.substr(i, hashmer), hashmer);
    }
    return array;
}

void delete_ptr(void *p) {

    delete [] p;
	return;
}

extern "C" {

    unsigned int hashEncode(char *seq, unsigned int hashmer) {
        return hash_encode(seq, hashmer);
    }
    unsigned int *hashEncodeArray(char *seq, unsigned int hashmer) {
        return hash_array(seq, hashmer);
    }
    void deleteptr(void *p) {
        delete_ptr(p);
    }
}
