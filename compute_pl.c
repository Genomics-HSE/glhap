#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
// #include <htslib/hts_idx.h>

/*
    Data structure for storing INDEL info.
    We just keep position and two floating values here.
*/
typedef struct {
    int pos;        // 1-based position
    float val[2];   // The two log-scaled PL values (first, last)
} indel_info_t;

/*
    This function:
     1) Determines the chromosome length for "chrM" or "chrY".
     2) Allocates a [chrom_length][4] float array, initialized to -1.0.
     3) Opens a BCF file (and loads its index).
     4) Iterates over variants in the given chromosome:
         - If REF length != ALT length ⇒ treat as INDEL.
            * If ALT is longer, store in `insertions`.
            * Otherwise, store in `deletions`.
         - Otherwise ⇒ treat as SNP, fill in log_gl_scores[pos][A/T/G/C].
     5) At the end, scales all filled PL values by `-PL/10`.
     6) Returns allocated arrays to the caller.

    Note: The function uses 0-based indexing internally for positions,
    but often returns/stores 1-based positions in the insertion/deletion arrays.
*/

int get_chrom_length(const char *chrom) {
    if (strcmp(chrom, "chrM") == 0) return 16569;
    if (strcmp(chrom, "chrY") == 0) return 57227415;
    return 0; // unknown
}

/* Free the arrays allocated by compute_log_genotype_likelihood. */
void free_likelihoods(float **log_gl_scores, int chrom_length,
                     indel_info_t *inserts, indel_info_t *deletions)
{
    if (log_gl_scores) {
        for (int i = 0; i < chrom_length; i++) {
            free(log_gl_scores[i]);
        }
        free(log_gl_scores);
    }
    if (inserts)   free(inserts);
    if (deletions) free(deletions);
}



int compute_log_genotype_likelihood(
    const char *bcf_file_path,
    const char *chrom,
    float ***log_gl_scores_out, // pointer for 2D array [chrom_length][4]
    indel_info_t **inserts_out,
    int *n_inserts_out,
    indel_info_t **deletions_out,
    int *n_deletions_out
)
{
    // 1. Determine chromosome length
    int chrom_length = 0;
    if (strcmp(chrom, "chrM") == 0) {
        chrom_length = 16569;
    } else if (strcmp(chrom, "chrY") == 0) {
        chrom_length = 57227415;
    } else {
        fprintf(stderr, "Unsupported chromosome: %s\n", chrom);
        return 1;
    }

    // 2. Allocate [chrom_length][4], init to -1.0
    float **log_gl_scores = (float **)calloc(chrom_length, sizeof(float *));
    if (!log_gl_scores) {
        fprintf(stderr, "Memory allocation failed for log_gl_scores.\n");
        return 1;
    }

    for (int i = 0; i < chrom_length; i++) {
        log_gl_scores[i] = (float *)calloc(4, sizeof(float));
        if (!log_gl_scores[i]) {
            fprintf(stderr, "Memory allocation failed for log_gl_scores row.\n");
            return 1;
        }
        // Initialize to -1.0
        for (int j = 0; j < 4; j++) {
            log_gl_scores[i][j] = -1.0f;
        }
    }

    // For insertion/deletion arrays, we use dynamic reallocation
    indel_info_t *inserts   = NULL;
    indel_info_t *deletions = NULL;
    int inserts_capacity = 0, n_inserts = 0;
    int deletions_capacity = 0, n_deletions = 0;

    // 3. Open BCF and read header
    htsFile *fp = bcf_open(bcf_file_path, "r");
    if (!fp) {
        fprintf(stderr, "Error opening BCF file: %s\n", bcf_file_path);
        return 1;
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        fprintf(stderr, "Error reading header from BCF.\n");
        bcf_close(fp);
        return 1;
    }

    // Load the index
    hts_idx_t *idx = bcf_index_load(bcf_file_path);
    if (!idx) {
        fprintf(stderr, "Error loading index for BCF: %s\n", bcf_file_path);
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return 1;
    }

    // 4. Create an iterator to fetch records from the given chromosome
    char region_str[256];
    snprintf(region_str, sizeof(region_str), "%s", chrom);
    hts_itr_t *itr = bcf_itr_querys(idx, hdr, region_str);
    if (!itr) {
        fprintf(stderr, "No records found for %s or cannot create iterator.\n", region_str);
        // Not necessarily an error if region just has no variants
        hts_idx_destroy(idx);
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        // Return empty results
        *log_gl_scores_out = log_gl_scores;
        *inserts_out = NULL;      *n_inserts_out = 0;
        *deletions_out = NULL;    *n_deletions_out = 0;
        return 0;
    }

    bcf1_t *rec = bcf_init();
    if (!rec) {
        fprintf(stderr, "Failed to allocate bcf1_t.\n");
        hts_itr_destroy(itr);
        hts_idx_destroy(idx);
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return 1;
    }

    // We'll need a buffer for PL data
    int32_t *pl_arr = NULL;
    int pl_arr_size = 0;
    
    // 5. Iterate over each variant record
    while (bcf_read(fp, hdr, rec)==0) {
        // fprintf(stderr, "SUCCESS");  
        bcf_unpack(rec, BCF_UN_FMT | BCF_UN_STR);
        // fprintf(stderr, "pos - %lld\n", rec->pos);
        // Single-sample assumption
        int n_sample = bcf_hdr_nsamples(hdr);
        if (n_sample < 1) {
            // No samples => skip
            continue;
        }

        // Get the PL format array
        int nPL = bcf_get_format_int32(hdr, rec, "PL", &pl_arr, &pl_arr_size);
        // printf("%d", nPL);
        if (nPL <= 0) {
            // No PL => skip
            continue;
        }

        // Check if it's an INDEL by comparing REF vs first ALT
        // (If there's no ALT, skip).
        if (rec->n_allele < 2) {
            // e.g. monomorphic site with no ALT
            continue;
        }

        char *ref = rec->d.allele[0];
        char *alt = rec->d.allele[1];
        int ref_len = (int)strlen(ref);
        int alt_len = (int)strlen(alt);
        // bcf_get_variant_types(rec)
        int is_indel = (ref_len != alt_len) && (strcmp(alt, "<*>") != 0);

        if (is_indel) {
            // We'll store first_val = PL[0], last_val = PL[nPL-1],
            // scaled => -PL/10
            float first_val = -1.0f;
            if (pl_arr[0] != bcf_int32_missing)
                first_val = -((float)pl_arr[0]) / 10.0f;

            float last_val = -1.0f;
            if (pl_arr[nPL - 1] != bcf_int32_missing)
                last_val = -((float)pl_arr[nPL - 1]) / 10.0f;

            // Insertion if alt is longer; else deletion
            if (alt_len > ref_len) {
                // Expand insertion array if needed
                if (n_inserts == inserts_capacity) {
                    inserts_capacity = (inserts_capacity == 0) ? 8 : inserts_capacity * 2;
                    inserts = (indel_info_t *)realloc(inserts, inserts_capacity * sizeof(indel_info_t));
                    if (!inserts) {
                        fprintf(stderr, "Memory allocation failed for inserts.\n");
                        // Cleanup & return
                        break;
                    }
                }
                inserts[n_inserts].pos = rec->pos + 1; // store 1-based
                inserts[n_inserts].val[0] = first_val;
                inserts[n_inserts].val[1] = last_val;
                n_inserts++;
            } else {
                // Expand deletion array if needed
                if (n_deletions == deletions_capacity) {
                    deletions_capacity = (deletions_capacity == 0) ? 8 : deletions_capacity * 2;
                    deletions = (indel_info_t *)realloc(deletions, deletions_capacity * sizeof(indel_info_t));
                    if (!deletions) {
                        fprintf(stderr, "Memory allocation failed for deletions.\n");
                        break;
                    }
                }
                deletions[n_deletions].pos = rec->pos + 1; // 1-based
                deletions[n_deletions].val[0] = first_val;
                deletions[n_deletions].val[1] = last_val;
                n_deletions++;
            }
            continue; // skip SNP logic
        }

        // Otherwise, treat as SNP: fill the row in log_gl_scores
        int pos0 = rec->pos; // 0-based already
        // fprintf(stderr, "%d", pos0);
        if (pos0 < 0 || pos0 >= chrom_length) {
            continue; // out of range
        }

        // The Python code walks over each allele with an index pattern.
        // We'll do something simpler: just loop over rec->n_allele, applying PL.
        // Then scale all PL values at the end. For now, we store raw PL in the array.
        int allele_pl_index = 0;
        int increment = 2; // the Python snippet had this pattern

        for (int i = 0; i < rec->n_allele; i++) {
            if (allele_pl_index >= nPL) break; // guard

            char *allele = rec->d.allele[i];
            float current_pl = (pl_arr[allele_pl_index] == bcf_int32_missing)
                                   ? -1.0f
                                   : (float)pl_arr[allele_pl_index];

            // A=0, T=1, G=2, C=3
            if (strcmp(allele, "A") == 0) {
                log_gl_scores[pos0][0] = current_pl;
            } else if (strcmp(allele, "T") == 0) {
                log_gl_scores[pos0][1] = current_pl;
            } else if (strcmp(allele, "G") == 0) {
                log_gl_scores[pos0][2] = current_pl;
            } else if (strcmp(allele, "C") == 0) {
                log_gl_scores[pos0][3] = current_pl;
            } else if (strcmp(allele, "<*>") == 0) {
                // fill any -1.0 in that row
                for (int j = 0; j < 4; j++) {
                    if (log_gl_scores[pos0][j] < 0.0f) {
                        log_gl_scores[pos0][j] = current_pl;
                    }
                }
            }

            allele_pl_index += increment;
            increment++;
        }
    }

    // Cleanup iteration
    if (pl_arr) free(pl_arr);
    bcf_destroy(rec);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    // 6. Scale final PL: => -PL/10
    //    If the position is still -1.0, we leave it that way (missing).
    for (int i = 0; i < chrom_length; i++) {
        for (int j = 0; j < 4; j++) {
            if (log_gl_scores[i][j] >= 0.0f) {
                log_gl_scores[i][j] = -log_gl_scores[i][j] / 10.0f;
            }
        }
    }

    // 7. Set outputs
    *log_gl_scores_out = log_gl_scores;
    *inserts_out = inserts;   *n_inserts_out = n_inserts;
    *deletions_out = deletions; *n_deletions_out = n_deletions;

    return 0;
}

// int main(int argc, char *argv[])
// {
//     if (argc < 3) {
//         fprintf(stderr, "Usage: %s <input.bcf> <chrom>\n", argv[0]);
//         return 1;
//     }
//     const char *bcf_file_path = argv[1];
//     const char *chrom = argv[2];

//     float **log_gl_scores = NULL; // [chrom_length][4]
//     indel_info_t *insertions = NULL;
//     indel_info_t *deletions  = NULL;
//     int n_inserts = 0;
//     int n_deletions = 0;

//     int ret = compute_log_genotype_likelihood(
//         bcf_file_path,
//         chrom,
//         &log_gl_scores,
//         &insertions,
//         &n_inserts,
//         &deletions,
//         &n_deletions
//     );

//     if (ret != 0) {
//         fprintf(stderr, "Error computing genotype likelihoods.\n");
//         return 1;
//     }

//     // Example: print a few results
//     printf("Number of insertions: %d\n", n_inserts);
//     for (int i = 0; i < n_inserts && i < 5; i++) {
//         printf("Insertion #%d at pos=%d: [%f, %f]\n",
//                i, insertions[i].pos,
//                insertions[i].val[0], insertions[i].val[1]);
//     }

//     printf("Number of deletions: %d\n", n_deletions);
//     for (int i = 0; i < n_deletions && i < 5; i++) {
//         printf("Deletion #%d at pos=%d: [%f, %f]\n",
//                i, deletions[i].pos,
//                deletions[i].val[0], deletions[i].val[1]);
//     }

//     // Print log_gl_scores for the first 10 positions
//     printf("First 10 positions of log_gl_scores:\n");
//     for (int i = 0; i < 10; i++) {
//         printf("pos=%d: A=%.3f, T=%.3f, G=%.3f, C=%.3f\n",
//                i+1,
//                log_gl_scores[i][0],
//                log_gl_scores[i][1],
//                log_gl_scores[i][2],
//                log_gl_scores[i][3]);
//     }

//     // TODO: free memory. 
//     //  e.g., if we know chrom length (e.g. 16569 for chrM),
//     //  we can safely do:
//     //    for(int i=0; i<16569; i++) free(log_gl_scores[i]);
//     //    free(log_gl_scores);
//     // free(insertions);
//     // free(deletions);

//     return 0;
// }
