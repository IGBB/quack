#include "seq.h"

#include <math.h>
#include <stdio.h>
#include <zlib.h>

#include "kseq.h"

#define unlikely(x) __builtin_expect ((x), 0)
#define likely(x)       __builtin_expect((x),1)


/* Convert ASCII to Integer for A T C and G, grouping A-T and G-C */
/*                A     C           G                                      T */
int lookup[20] = {0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

KSEQ_INIT(gzFile, gzread)

int* read_adapters(char *adapters_file) {
    int kmer_size = 10;
    int array_size = pow(4, kmer_size);
    gzFile fp;
    kseq_t *seq;
    int i, l, index;
    fp = gzopen(adapters_file, "r");
    seq = kseq_init(fp);
    int *kmers = malloc(array_size*sizeof(int));
    memset(kmers, 0, array_size*sizeof(int));
    while ((l = kseq_read(seq)) >= 0) {
        index = 0;
            for (i = 0; i < kmer_size; i++) {
                index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
            }
        for (; i < seq->seq.l; i++) {
            index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
            kmers[index] = 1;
        }

    }
    kseq_destroy(seq);
    gzclose(fp);
    return kmers;
}

sequence_data* read_fastq(char *fastq_file, int *kmers) {
    gzFile fp;
    kseq_t *seq;
    int i, l, index;
    int kmer_size = 10;
    int array_size = pow(4, kmer_size);
    int max_length = 0;
    fp = gzopen(fastq_file, "r");
    seq = kseq_init(fp);
    base_information *bases = NULL;
    sequence_data *to_return = malloc(sizeof(sequence_data));
    int number_of_sequences = 0;

    while ((l = kseq_read(seq)) >= 0) {
        if (unlikely(seq->seq.l > max_length)) {
            bases = realloc(bases, seq->seq.l*sizeof(base_information));
            memset(bases+max_length, 0, (seq->seq.l - max_length)*sizeof(base_information));
            max_length = seq->seq.l;
        }
        for (i = 0; i < seq->seq.l; i++) {
            int base = seq->seq.s[i];
            int offset = lookup[base-65 & ~32];
            bases[i].content[offset]++;
            int quality = seq->qual.s[i]-33;
            bases[i].scores[quality]++;
        }
        index = 0;
        for (i = 0; i < kmer_size; i++) {
            index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
        }
        if (kmers) {
            for (; kmers[index] == 0 && i < seq->seq.l; i++) {
                index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
            }
        }
        if (i < seq->seq.l) {
            bases[i].kmer_count++;
        }

        bases[seq->seq.l-1].length_count++;
        number_of_sequences++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    to_return->bases = bases;
    to_return->max_length = max_length;
    to_return->number_of_sequences = number_of_sequences;
    return to_return;
}

sequence_data* transform(sequence_data* data) {
    int i, j;
    data->original_max_length = data->max_length;
    // binning
    if (data->max_length > 3000) {
        fprintf(stderr, "Binning...\n");
        int bin_size = 100;
        int unbinned;
        int binned = 0;

        for (unbinned = 1; unbinned < data->max_length; unbinned++) {
             if (unbinned%bin_size == 0){
                binned++;
                for (i = 0; i < 4; i++) {
                    data->bases[binned].content[i] = 0;
                }
                for (i = 0; i < 91; i++) {
                    data->bases[binned].scores[i] = 0;
                }
                data->bases[binned].length_count = 0;
             }
             for (i = 0; i < 4; i++) {
                data->bases[binned].content[i] = data->bases[binned].content[i] + data->bases[unbinned].content[i];
             }
             for (i = 0; i < 91; i++) {
                data->bases[binned].scores[i] = data->bases[binned].scores[i] + data->bases[unbinned].scores[i];
             }
            // fprintf(stderr, "%d\n", data->bases[binned].length_count);
            data->bases[binned].length_count = data->bases[binned].length_count + data->bases[unbinned].length_count;
            data->bases[binned].kmer_count = data->bases[binned].kmer_count + data->bases[unbinned].kmer_count;
        }
        data->max_length = binned;
    }

    for (i = 1; i < data->max_length; i++) {
        data->bases[i].kmer_count = data->bases[i-1].kmer_count + data->bases[i].kmer_count;
    }

    // transforming data for drawing (percentages rather than counts)
    for (i = 0; i < data->max_length; i++) {
        int content_sum = 0;
        int score_sum = 0;
        /* for (j = 0; j < 4; j++) { */
        /*     content_sum = content_sum + data->bases[i].content[j]; */
        /* } */
        /* if (content_sum != 0) { */
        /*     for (j = 0; j < 4; j++) { */
        /*         data->bases[i].content[j] = 100*data->bases[i].content[j]/content_sum; */
        /*     } */
        /* } */
        for (j = 0; j < 91; j++) {
            score_sum = score_sum + data->bases[i].scores[j];
        }
        if (score_sum != 0) {
            for (j = 0; j < 91; j++) {
                data->bases[i].scores[j] = 100*data->bases[i].scores[j]/score_sum;
            }
        }
        data->bases[i].length_count = ceil(100*(float)data->bases[i].length_count/data->number_of_sequences);
        data->bases[i].kmer_count = ceil(100*(float)data->bases[i].kmer_count/(float)data->number_of_sequences);

    }
    return data;
}
