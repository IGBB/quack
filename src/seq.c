#include "seq.h"

#include <stdio.h>
#include <zlib.h>

#include "kseq.h"

#define unlikely(x) __builtin_expect ((x), 0)
#define likely(x)       __builtin_expect((x),1)


/* Convert ASCII to Integer for A T C and G, grouping A-T and G-C. Padding
   database to next power of 2 to make it easy to bind index to array bounds */
/*                A     C           G                                      T */
int lookup[32] = {0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* Hash function to encode A,T,C,G to 2-bit integer
     base - 65   = Removes the bottom section of ASCII table. Makes 'A' = 0
     ... & 31    = Converts lowercase letters to uppercase and forces index 
                     within lookup array bounds
     lookup[...] = Translates ascii ATCG to integer
*/
#define encode_base(base) (lookup[(base - 65) & 31])

/* Create next kmer
                                 / Make room for newest base by shift kmer by 2
                                 |             / endcode base to 2-bit integer
                                 |             |                 / Remove oldest base
                                 |             |                 |
 */
#define next_kmer(kmer, base) (((kmer << 2) + encode_base(base)) & ADAPTER_DB_MASK)

KSEQ_INIT(gzFile, gzread)

int* read_adapters(char *adapters_file) {
    gzFile fp;
    kseq_t *seq;
    int i, l, index;

    /* Create kmer hash database, making sure it's zero'd */
    int *kmers = calloc(ADAPTER_DB_SIZE, sizeof(int));
    if(kmers == NULL){
      perror("Can't create adapter database");
      exit(EXIT_FAILURE);
    }
    
    /* Open adapter file */
    fp = gzopen(adapters_file, "r");
    /* Print error and exit if adapter file can't be opened */
    if(fp == NULL){
      perror("Can't open adapter file");
      exit(EXIT_FAILURE);
    }
    
    /* Loop through each adapter sequence, storing each adapter kmer hash in
       database */
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      /* Clear index kmer */
      index = 0;

      /* Fill index with first (kmer-1) bases */
      for (i = 0; i < ADAPTER_KMER_SIZE-1; i++) {
        index = next_kmer(index, seq->seq.s[i]);
      }

      /* Loop through remaining bases of adapter, adding sliding kmer to
         database */
      for (; i < seq->seq.l; i++) {
        index = next_kmer(index, seq->seq.s[i]);
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
    int max_length = 0;
    base_information *bases = NULL;
    sequence_data *to_return = malloc(sizeof(sequence_data));
    int number_of_sequences = 0;


    fp = gzopen(fastq_file, "r");
    /* Print error and exit if adapter file can't be opened */
    if(fp == NULL){
      perror("Can't open fastq file");
      exit(EXIT_FAILURE);
    }

    /* Parse each sequence in fastq */
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {

      /* if current sequence is longer than any sequence previously read,
         increase base information size and set all new base counts to zero.
         Shouldn't happen very much unless input file is sorted by length.*/
      if (unlikely(seq->seq.l > max_length)) {
        bases = realloc(bases, seq->seq.l*sizeof(base_information));
        memset(bases+max_length, 0, (seq->seq.l - max_length)*sizeof(base_information));
        max_length = seq->seq.l;
      }

      /* Save base and quality info 
           encode base to 2-bit integer
           
           There are currently two encodings for quality score: phred + 33 and
           phred + 64. We will assume the encoding is phred33 and correct later
           if the assumption is wrong.
       */
      for (i = 0; i < seq->seq.l; i++) {
        bases[i].content[encode_base(seq->seq.s[i])]++;
        bases[i].scores[seq->qual.s[i] - 33]++;
      }

      /* Search kmer database if it exists */
      if(kmers) {
        /* Clear index kmer */
        index = 0;

        /* Fill kmer index with first bases */
        for (i = 0; i < ADAPTER_KMER_SIZE; i++) {
          index = next_kmer(index, seq->seq.s[i]);
        }

        /* Loop through remaining bases, searching kmer database, stoping if kmer is found */
        for (; !kmers[index] && i < seq->seq.l; i++) {
          index = next_kmer(index, seq->seq.s[i]);
        }

        /* Add to base kmer count if i didn't reach the end of the sequence */
        if (i < seq->seq.l) {
            bases[i].kmer_count++;
        }

      }

      /* record length and number of sequences */
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
