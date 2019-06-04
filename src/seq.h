#ifndef __SEQ_H
#define __SEQ_H

#include <stdint.h>

typedef struct {
    uint64_t scores[91];
    uint64_t content[4];
    uint64_t length_count;
    uint64_t kmer_count;
} base_information;

typedef struct {
    base_information *bases;
    uint64_t max_length;
    uint64_t original_max_length;
    uint64_t number_of_sequences;
} sequence_data;

int* read_adapters(char *adapters_file);
sequence_data* read_fastq(char *fastq_file, int *kmers);
sequence_data* transform(sequence_data* data);

#endif
