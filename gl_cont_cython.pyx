
from tqdm import tqdm
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from IPython.display import clear_output
cimport cython
import numpy as np
from cython.parallel import prange
from libc.math cimport pow
from scipy.special import binom


def get_num_reads(str bam_fname):
    ''''
    This function calculate mapped reads
    '''
    samfile = pysam.AlignmentFile(bam_fname, "rb")
    num_reads = 0
    for read in samfile.fetch('chrM'):
        if read.is_mapped:
            num_reads += 1
    samfile.close()
    return num_reads


def bam2consensus(bam_fname, double ac_threshold=0, double af_threshold=0):
    '''
    Make consensus fasta
    
    Parameters:
        ref_name : str
        path to a file with reference
        bam_fname : str
        path to a bam file
        ac_threshold : float
        af_threshold : float
    '''
    cdef str consensus = ''
    cdef int max_count, total_count
    cdef str allele
    cdef str max_allele
   
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
            else:
                consensus += 'N'

    return consensus.replace('-','')


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(True)
def get_MN(char[:, :] genomes,str bam_fname, long[:] aln_coords, same, trunc = 0, verbosity = False):
    '''
    Make matrices and base error for the method.
    
    Parameters:
        genomes: char
        vector where on j row is j genome.
        
        bam_fname: str
        bam file
        
        same : dict
        dictionary of positions and bases where all genomes have the same base.
    
    Output:
        M : np.ndarray[float, float]
            matrix where M[i, j] is number of bases where read i has same bases with genome 
        N : np.ndarray[float, float]
            matrix where N[i, j] is number of bases where read i has different bases with genome j

    
    '''
    bam = pysam.AlignmentFile(bam_fname, "rb")
    cdef double[:, :] M, N
    cdef int k, i, j, pos, offset
    offset = 0
    cdef double correct, incorrect, P_cor
    cdef str seq
    cdef long num_reads = get_num_reads(bam_fname)
    cdef int num_genomes = genomes.shape[0]
    cdef double base_incorr = 0
    cdef double base_total = 0
    N = np.zeros((num_reads, num_genomes))
    M = np.zeros((num_reads, num_genomes))
    i = 0
    j = 0
    for read in tqdm(bam.fetch('chrM'), total = bam.count(), desc = 'MN tables'):
        
        if not read.is_mapped:
            continue
            
        seq = read.query_sequence
        pos = read.reference_start
        
        if read.cigartuples[0][0] == 4: #read is soft clipped
            left_trim = read.cigartuples[0][1]
            seq = seq[left_trim:]
            
        # if read.cigartuples[-1][0] == 4: #read is soft clipped
        #     right_trim = read.cigartuples[-1][1]
        #     seq = seq[:-right_trim]
        
        
        oldest_pos = pos
        
        qual = read.query_qualities
        
        ins_pos = 0
        
        for j in range(num_genomes):
            if "I" in read.cigarstring: # Check for indels
                if j == 0:
                    # print('INDEL')
                    M[i, j] = -1
                    N[i, j] = -1
                else:
                    cigar = read.cigartuples

                    for m in range(len(cigar)):
                        if cigar[m][0] == 1:
                            for s in range(m):
                                ins_pos += cigar[s][1]
                            break
                            
                    # if pos + ins_pos < genomes.shape[1] and chr(genomes[s, aln_coords[pos] + ins_pos]) != seq[aln_coords[pos]+ins_pos]:
                    #     M[i, j] = -1
                    #     N[i, j] = -1
            
            
            
            if "D" in read.cigarstring: # Check for indels
                if j == 0:
                    M[i, j] = -1
                    N[i, j] = -1
                else: # пока нуждается в доработке но посмотрим
                    cigar = read.cigartuples
                    ins_pos = 0
                    for m in range(len(cigar)):
                        if cigar[m][0] == 1:
                            for s in range(m):
                                ins_pos += cigar[s][1]
                            # break
                    if chr(genomes[s, aln_coords[pos+ins_pos]]) != '-':
                        M[i, j] = -1
                        N[i, j] = -1
                        

                        
            correct = 0
            incorrect = 0
            offset = 0
            debug_str = ''

            for k in range(trunc,len(seq)-trunc):
                
                if k + pos + offset >= genomes.shape[1]:
                    break
                # if k + pos + offset >= len(aln_coords):
                    # break;
                    
                while chr(genomes[j, aln_coords[pos] + k +offset]).upper() == '-' or chr(genomes[0, aln_coords[pos] + k + offset]).upper() == '-':
                    if chr(genomes[j, aln_coords[pos] + k + offset]).upper() == '-' and chr(genomes[0, aln_coords[pos] + k + offset]).upper() == '-':
                        offset += 1
                    
                    elif chr(genomes[j, aln_coords[pos] + k + offset]).upper() == '-':
                        if 'D' not in read.cigarstring:
                            M[i, j] = -1
                            N[i, j] = -1
                            break
                        else:
                            correct += 1
                            offset += 1 
                        
                    elif chr(genomes[0, aln_coords[pos] + k + offset]).upper() == '-':
                        if 'I' not in read.cigarstring:
                            M[i, j] = -1
                            N[i, j] = -1
                            break
                        else:
                            break
                    if k + pos + offset >= genomes.shape[1]:
                        break
                            
                    # print(pos + k + offset)
                
                if chr(genomes[j, aln_coords[pos] + k + offset]).upper() == 'N' or seq[k] == 'N':
                    # correct += 0.25
                    # incorrect += 0.75
                    # correct +=1
                    pass
                                
                elif seq[k] == chr(genomes[j, aln_coords[pos] + k + offset]).upper(): #means that read has same base with j genome
                    P_cor = 1 - 10**(- qual[k]/10)
                    correct += P_cor
                    incorrect += 1 - P_cor
                    if verbosity:
                        debug_str += f'{seq[k]}, {chr(genomes[j][aln_coords[pos] + k + offset ]).upper()} => +1\n'
                    # correct += 1 #old version for debugging
                    if aln_coords[pos] + k +offset in same:
                        base_total += 1
                        
                else:
                    # print(seq[k], chr(genomes[j][k+pos]).upper()) #means that read has same difference with j genome
                    P_cor = (10**(- qual[k]/10))/3
                    incorrect += 1 - P_cor
                    correct += P_cor
                    if verbosity:
                        debug_str += f'{seq[k]}, {chr(genomes[j][aln_coords[pos] + k + offset]).upper()} => -1\n'
                    # incorrect += 1
                    if aln_coords[pos] + k + offset in same:
                        # print("plus incorrect")
                        base_total += 1
                        base_incorr += 1
                        

            if M[i, j] != -1:
                M[i, j] = correct
                N[i, j] = incorrect

                
                
                
        
        i += 1
        
    bam.close()
    cdef double base_err = base_incorr / base_total
    return np.array(M, dtype=np.float64), np.array(N, dtype=np.float64), base_err

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def get_mc(double[:, ::1] m, double[:, ::1] n, double eps):
    ''''
    Преобразует M N и base_err в матрицу mc,  которая в дальнейшем используется для остальных вычислений
    '''
    cdef double[:,::1] mc
    cdef long num_reads = m.shape[0]
    cdef long num_genomes = m.shape[1] 
    cdef long j
    mc = np.zeros((num_reads, num_genomes))
    for i in range(num_reads):
        for j in range(num_genomes):
            if m[i, j] == -1:
                mc[i, j] = 0
            else:
                mc[i, j] = binom(m[i, j] + n[i, j], m[i, j]) * (1 - eps)**(m[i, j]) * eps**n[i, j]
    return np.asarray(mc)
    




@cython.cdivision(True)
# @cython.wraparound(False)
# @cython.boundscheck(False)
def get_Zi(double[:,::1] mc, double[::1] p,double eps, long i):
    cdef long num_reads = mc.shape[0]
    cdef long num_genomes = mc.shape[1] 
    cdef long j
    cdef long Z
    cdef double[:] probs
    cdef double s
    s = 0
    probs = np.zeros(num_genomes, dtype = float)
    

    for j in range(num_genomes):
        probs[j] = mc[i, j] * p[j]
        s += probs[j]

        
    for j in range(num_genomes):
        probs[j] =probs[j] / s
        
    Z = random_choice(probs)

    return Z


@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def get_eta(long[:] z,int num_genomes):
    '''
    create vector eta
    eta[j] is number of reads that were predicted to be from j-th genome 
    '''
    
    cdef long[:] eta
    eta = np.zeros(num_genomes, dtype = int)
    
    cdef int num_reads
    num_reads = z.shape[0]
    
    cdef int i
    
    for i in range(num_reads):
        eta[z[i]] += 1
        
    return np.array(eta)


cdef extern from "stdlib.h":
    double drand48()
    void srand48(long int seedval)

    
@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def random_choice(double[:] probs):
    '''Returns random number from 0 to n-1 according to probs
        Check if probs do not sum up to 1!'''
    
    cdef double s = 0
    cdef int i
    cdef int l = len(probs)
    
    

    
    cdef double x = drand48()
    cdef double cum_probs = 0
    cdef long n = 0
    while x > cum_probs:
        cum_probs += probs[n]
        n += 1
    n -= 1
    return n
