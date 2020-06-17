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

char rev_lookup[] = "ATGC";

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

sequence_data* read_fastq(char *fastq_file, int *kmers, encoding_t encoding) {
    gzFile fp;
    kseq_t *seq;
    int i, j, l, index;
    int max_length = 0;
    base_information *bases = NULL;
    sequence_data *data = malloc(sizeof(sequence_data));
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

    data->bases = bases;
    data->max_length = max_length;
    data->number_of_sequences = number_of_sequences;

    data->avg_score = malloc(sizeof(float) * max_length);

    /* Transform data
     * - Make kmer_count cummulative
     * - Make length_count cummulative
     * - Guess or validate encoding
     * - Calculate average score of each base,
     *      adjusting for length drop off */
    uint64_t kmer_count = 0;
    uint64_t length_count = 0;
    int min = 91;
    int max = 0;


    /*Adjust length_count. REVERSE  */
    for(i = max_length - 1; i >= 0; i--){
        length_count += bases[i].length_count;
        bases[i].length_count = length_count;
    }
   
    for(i = 0; i < max_length; i++){
        /* Adjust kmer_count */
        kmer_count += bases[i].kmer_count;
        bases[i].kmer_count = kmer_count;



        /* Search the scores below min for non-zero */
        for(j = 0; bases[i].scores[j] == 0 && j < min; j++);
        min = j;
        /* Search the scores above max for non-zero */
        for(j = 90; bases[i].scores[j] == 0 && j > max; j--);
        max = j;

    }

    /* Guess/Validate encoding from min-max range
     * NOTE: all scores are alread offset by 33
     - If min offset is below 31, encoding must be phred33
     - If max offset is below 45 (phred64 score 13), encoding is probably phred33 */
    if(min < 31 || max < 45){
        data->encoding = phred33;
    }else{
        data->encoding = phred64;
    }

    /* Fail if encoding is set and doesn't match guess */
    if( encoding != guess && data->encoding != encoding ){
        perror("Can't open fastq file");
        exit(EXIT_FAILURE);
    }

    if(data->encoding == phred64){
        /* Move score data so that 0 is the first possible score. Remove possibility
         * for negative phred score */
        for(i = 0; i < max_length; i++){
            for(j = min; j <= max; j++){
                bases[i].scores[j - 31] = bases[i].scores[j];
                bases[i].scores[j] = 0;
            }
        }

        /* Adjust min and max score accordingly */
        min -= 31;
        max -= 31;
    }

    data->max_score = max;
    data->min_score = min;

    /* Calculate avg_score */

    for(i = 0; i < max_length; i++){
        data->avg_score[i] = 0;
        for(j = min; j <= max; j++)
            data->avg_score[i] += (float)bases[i].scores[j] * j;
        data->avg_score[i] = data->avg_score[i] / data->bases[i].length_count;
    }
    return data;
}

