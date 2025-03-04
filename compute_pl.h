#ifndef COMPUTEPL_H
#define COMPUTEPL_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int pos;       // 1-based position
    float val[2];  // first and last scaled PL
} indel_info_t;

/*
 * Return chromosome length for "chrM" or "chrY" (0 if unknown).
 */
int get_chrom_length(const char *chrom);

/*
 * Main function that allocates arrays and computes results.
 * 
 * log_gl_scores_out -> pointer to a float** array of shape [chrom_length][4]
 * inserts_out, deletions_out -> arrays of indel_info_t
 * n_inserts_out, n_deletions_out -> number of items in the arrays
 * Return 0 on success, non-0 on error.
 */
int compute_log_genotype_likelihood(
    const char *bcf_file_path,
    const char *chrom,
    float ***log_gl_scores_out,
    indel_info_t **inserts_out,
    int *n_inserts_out,
    indel_info_t **deletions_out,
    int *n_deletions_out
);

/*
 * Free the memory allocated by compute_log_genotype_likelihood.
 */
void free_likelihoods(float **log_gl_scores, int chrom_length,
                     indel_info_t *inserts, indel_info_t *deletions);

#ifdef __cplusplus
}
#endif

#endif
