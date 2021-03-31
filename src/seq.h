#ifndef __SEQ_H
#define __SEQ_H

#include <stdint.h>

/* Define adapter kmer database constants 

     ADAPTER_KMER_SIZE = length of each kmer. With length 10, the chance of a
     random match is around one in a million. Since the data is compressed by
     the number of pixels in the graph, the random hits should be masked.

     ADAPTER_DB_SIZE = The maximum number of kmers possible with a given kmer
     length (4^length). Bit shifting can be used in place of the math::pow
     function since the number of bases we encode is a power of 2, removing the
     need for the math library. To use bit shifting, the ADAPTER_KMER_SIZE needs
     to be multpled by the number of bits each base encoding takes (i.e. 2)

     ADAPTER_DB_MASK = Mask for the kmer hash. Since the db size will always be
     a power of 2, we can subtract 1 to get a mask of all 1s.*/
#define ADAPTER_KMER_SIZE 10
#define ADAPTER_DB_SIZE (1 << (ADAPTER_KMER_SIZE * 2))
#define ADAPTER_DB_MASK (ADAPTER_DB_SIZE - 1)

extern int lookup[];
extern char rev_lookup[];

typedef enum {
    guess = 0,
    phred33 = 33,
    phred64 = 64
} encoding_t;

typedef struct {
    uint64_t scores[91];
    uint64_t avg_scores[91];
    uint64_t content[4];
    uint64_t length_count;
    uint64_t kmer_count;
} base_information;

typedef struct {
    base_information *bases;
    uint64_t max_length;
    uint64_t number_of_sequences;
    int max_score;
    int min_score;
    encoding_t encoding;
    float * avg_score;
    int read_type;
} sequence_data;

int* read_adapters(char *adapters_file);
sequence_data* read_fastq(char *fastq_file, int *kmers, encoding_t encoding);
/* sequence_data* transform(sequence_data* data); */

#endif
