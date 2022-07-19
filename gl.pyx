
from tqdm import tqdm
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
cimport cython
import numpy as np
from cython.parallel import prange
from libc.math cimport pow
from scipy.special import binom

def bam2consensus(
        ref_fname, bam_fname, double ac_threshold=0, double af_threshold=0):

    cdef str consensus = ''
    cdef int max_count, total_cgeount
    cdef str allele
    cdef str max_allele
    for record in SeqIO.parse(ref_fname, "fasta"):
        assert record.id == 'chrM'
        # consensus = "N" * len(record)

    with pysam.AlignmentFile(bam_fname, "rb") as bam:
        allele_counter = Counter()
        for pileup_column in tqdm(bam.pileup(), total=16569, desc = 'consensus dna'):
            assert pileup_column.reference_name == 'chrM'
            pos = pileup_column.reference_pos

            allele_counter.clear()
            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del:
                    allele = "-"
                else:
                    allele = pileup_read.alignment.query_sequence[
                        pileup_read.query_position]
                allele_counter[allele] += 1

            max_allele = "N"
            max_count, total_count = 0, 0
            for allele, count in allele_counter.items():
                if count > max_count:
                    max_count = count
                    max_allele = allele
                total_count += count

            assert max_allele in "ACGTN-"
            if (max_count >= ac_threshold and
                max_count / total_count >= af_threshold):
                consensus += max_allele

    return consensus


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def get_MN(int num_reads, int num_genomes,char[:, :] genomes, bam):
    cdef double[:, :] N, M
    cdef int correct, incorrect, k, i, j, pos
    cdef str seq
    N = np.zeros((num_reads, num_genomes))
    M = np.zeros((num_reads, num_genomes))
    i = 0
    j = 0
    ### ЗАПАРАЛЛЕЛИТЬ
    for read in tqdm(bam.fetch('chrM'), total = bam.count(), desc = 'MN tables'):
        seq = read.query_sequence
        pos = read.reference_start
        if not read.is_mapped:
            continue
        for j in range(num_genomes):
            correct = 0
            incorrect = 0
            for k in range(len(seq)):
                if k+pos > 16568:
                    break
                if seq[k] == chr(genomes[j][k+pos]).upper():
                    correct += 1
                else:
                    # print(seq[k], chr(genomes[j][k+pos]).upper())
                    incorrect += 1
            M[i, j] = correct
            N[i, j] = incorrect
        i += 1
    return np.array(M, dtype=np.float64), np.array(N, dtype=np.float64)





@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def get_Z(long[:, ::1] m, long[:, ::1] n, double[::1] p,double eps):
    cdef long num_reads = m.shape[0]
    cdef long num_genomes = m.shape[1] 
    cdef long i, j
    cdef long[:] Z
    cdef double[:] probs
    cdef double s
    Z = np.empty(num_reads, dtype=int)
    for i in tqdm(range(num_reads), desc = 'assignments'):
        s = 0
        probs = np.zeros(num_genomes, dtype = float)
        # with nogil:
        for j in range(num_genomes):
            probs[j] = binom(m[i, j] + n[i, j], m[i, j]) * pow((1 - eps),(m[i, j])) * pow(eps,(n[i, j])) * p[j]
            s += probs[j]
        
        for j in range(num_genomes):
            probs[j] =probs[j] / s
        
        Z[i] = np.random.choice(np.arange(0, num_genomes), p = probs)
    return Z


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def get_eta(long[:] z,int num_genomes):
    cdef long[:] eta
    eta = np.zeros(num_genomes, dtype = int)
    cdef int num_reads
    num_reads = z.shape[0]
    cdef int i
    for i in range(num_reads):
        eta[z[i]] += 1
    return np.array(eta)
