#!/usr/bin/env python
# coding: utf-8


import argparse

parser = argparse.ArgumentParser(description='Supply reference fasta and bam file')
parser.add_argument('ref',
                    help='reference fasta')
parser.add_argument('bam',
                    help='bam file')
parser.add_argument('contaminants',
                    help='fasta with all possible contaminants', default='contaminants.fa')
parser.add_argument('nIter',
                    help='number of iteration mcmc', default=10000)


args = parser.parse_args()
ref_fname = args.ref
bam_fname = args.bam
genomes_fname = args.contaminants
nIter = int(args.nIter)

# import of all neccessary modules
import os
from collections import Counter
import pysam
import numpy as np
from tqdm import tqdm
from scipy.special import binom
import scipy.stats as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocess import Pool
from gl_cont_cython import *



def get_base_err(bam_fname, ref, aln_pos, same_set):
    bam = pysam.AlignmentFile(bam_fname, "rb")
    correct = 0
    incorrect = 0
    for readId, read in enumerate(bam.fetch('chrM')):
        
        if not read.is_mapped or 'D' in read.cigarstring or 'I' in read.cigarstring:
            continue
            
        seq = read.query_sequence
        pos = read.reference_start
        
        if read.cigartuples[0][0] == 4: #read is soft clipped
            left_trim = read.cigartuples[0][1]
            seq = seq[left_trim:]
                        

                        
        
        offset = 0
        debug_str = ''

        for k in range(len(seq)):
            if aln_pos[pos+k] in same_set:
                if seq[k].upper() == ref[aln_pos[pos+k]]:
                    correct+=1
                else:
                    # print(pos, k, readId)
                    incorrect += 1
    return correct, incorrect, incorrect/(correct + incorrect)
                
        



def preprocess(ref_fname, genomes_fname, bam_fname):
    
    base = bam_fname[:-4]
    
    pysam.index(bam_fname);
    
    os.system(f"samtools view {bam_fname} chrM -o {base+'_mt.bam'}")
    base = base + '_mt'
    
    print('#EXTRACTING MTDNA OK')
    
    # Один из вариантов получения консенсуса, работает не очень хорошо.
    '''consensus = bam2consensus(ref_fname, bam_fname)
    consensus_fa = '>chrM\n'+''.join(consensus) +'\n'
    with open(f'{base}.fa', 'w') as new_genomes:
        new_genomes.write(consensus_fa)'''
    
    
    # os.system(f'samtools consensus -o {base}.fa {bam_fname}')
    
#     os.system(f'sh gatkconsensus.sh {ref_fname} {base}.bam  {base}1.fa')
#     os.system(f'''awk -i '/^>/{{print ">chrM"; next}}{{print}}' {base}1.fa > {base}.fa''')
#     os.system(f'rm {base}1.fa')
    
    # os.system(f'cp genome_0.fa {base}.fa')
    
    os.system(f'bcftools mpileup  -d 2000 -m 3 -C50 -q 30 -EQ 20 -f {ref_fname} {base}.bam | bcftools call -m --ploidy 1 > {base}.vcf')
    os.system(f'perl CnsMaj3_1.pl -i {base}.vcf -o {base}.fa -l 16569 -cov 1 -diff 0.5 -idiff 0.5 -h {base} -callindels no > {base}.cns')
    
    # return 0
    
    print('#CONSENSUS IS READY')
    
    os.system(f'cat {base}.fa {genomes_fname} > {base}_genomes.fa'); # gather consensus and possible contaminants together
    os.system(f'mafft {base}_genomes.fa >  {base}_aligned.fa') # do multiple alignment
    aligned_genomes = f'{base}_aligned.fa'
    print("#ALL GENOMES ARE READY")
    # new_cons = list(SeqIO.parse(f'{base}_aligned.fa', "fasta"))[0]
    # SeqIO.write(new_cons, f'{base}.real.fa', "fasta") #gives you reference after it have been realigned with MAFFT
    os.system(f'bwa index -a bwtsw {base}.fa') #indexing consensus
    os.system(f'samtools faidx {base}.fa')
    os.system(f'rm {base}.dict')
    os.system(f'picard CreateSequenceDictionary R={base}.fa O={base}.dict')
    os.system(f'samtools fastq {bam_fname} > {base}.fq')
    os.system(f'bwa aln -l 1000 -t 10 {base}.fa {base}.fq > {base}_ra.sai')
    os.system(f"bwa samse -r '@RG\\tID:{base}\\tLB:{base}_L1\\tPL:ILLUMINA\\tSM:{base}' {base}.fa {base}_ra.sai {base}.fq | samtools sort -O BAM -o {base}_ra.sort.bam")
    # os.system(f'picard MarkDuplicates I={base}_ra.sort.bam O={base}_ra.sort.rmdup.bam METRICS_FILE=metrics.txt TMP_DIR=temp REMOVE_DUPLICATES=true ASSUME_SORTED=true ') #VALIDATION_STRINGENCY=LENIENT
    
    total_reads = int(pysam.view('-c', bam_fname))
    need_reads = 64000
    proportion = need_reads / total_reads
    # print(proportion)

    # -s {123+proportion}
    os.system(f'samtools view  -O BAM {base}_ra.sort.bam | samtools sort > {base}_ra.sort.rmdup.bam')
    
    pysam.index(f'{base}_ra.sort.rmdup.bam');
    os.system(f'samtools calmd -Erb {base}_ra.sort.rmdup.bam {base}.fa > {base}_ra.final.bam 2>/dev/null');
    bam_final = f'{base}_ra.final.bam'
    os.system(f'samtools index {base}_ra.final.bam')
    os.system(f'rm {base}_ra.sai')
    print("#BAM FILE IS READY")
    return bam_final, aligned_genomes


# In[171]:


def make_genomes_arr(genomes_fname):
    genomes = list()
    for record in SeqIO.parse(genomes_fname, "fasta"):
        genomes.append(str(record.seq))
    genomes_arr = np.array([list(x) for x in genomes], dtype = 'S1')
    return genomes_arr


# In[172]:


def get_same(genomes_arr):
    same_positions = []
    for i in range(genomes_arr.shape[1]):
        if len(np.unique(genomes_arr[:,i])) == 1:
            same_positions.append(i)
    return set(same_positions)


# In[173]:


def get_base_err(bam_fname, same_dict):
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    same_positions = list(same_dict.keys())
    correct = 0
    total = 0
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    for pileupcolumn in tqdm(samfile.pileup("chrM")):
        pos = pileupcolumn.pos
        if pos not in same_positions:
            continue

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                total += 1
                # query position is None if is_del or is_refskip is set.
                nbase =  pileupread.alignment.query_sequence[pileupread.query_position]
                if nbase == same[pos].decode('ascii').upper():
                    correct += 1
    base_err = 1 - correct/total
    samfile.close()
    return base_err


# In[174]:


def get_aln_pos(reference):
    aln_coor = []
    for i in range(len(reference)):
        if reference[i] != '-':
            aln_coor.append(i)
            
    return np.asarray(aln_coor)


# In[175]:


def do_mcmc(n_iterations = 50000, output_file='', n_threads=8, model=0, show_each=10):
    if output_file != '':
        res = open(output_file,'w')
    num_reads, num_genomes  = MC.shape
    print(MC.shape)
    p = np.random.dirichlet([1]*num_genomes)
    # pool = Pool(n_threads)
    for i in tqdm(range(n_iterations) ):
        
        func = lambda x: get_Zi(MC, p, base_err, x)
        
        # Z = np.array(pool.map_async(func, range(num_reads)).get())
        Z = np.array([func(s) for s in range(num_reads) ])
        eta = get_eta(Z, num_genomes)
        if model == 0:
            p0 = np.random.beta(1 + eta[0],1+num_reads-eta[0])
            p_other = np.random.dirichlet(1+ eta[1:])
            p_other *= (1-p0)/p_other.sum()

            p[0] = p0
            p[1:] = p_other
        else:
            p = np.random.dirichlet(1+ eta)
        if output_file != '':
            res.write(f'iteration {i}')
            res.write(str(p[0]))
        if i % show_each == 0:
            # print(p[0], p[1:].sum()) 
            print(p)
    # pool.close()
    if output_file != '':
        res.close()
    return p


# In[176]:


bam, genomes = preprocess(ref_fname, genomes_fname, bam_fname)





# In[178]:


genomes_arr = make_genomes_arr(genomes)




# In[184]:


same = get_same(genomes_arr)


# In[185]:


genomes0 = (''.join( np.array(genomes_arr, dtype = str)[0])).upper()

# In[188]:


aln_coords = get_aln_pos(genomes0)




# In[190]:


M, N, base_err = get_MN(genomes_arr, bam, aln_coords, same)





# In[198]:





# При большой ошибке большое число ридов не картируется, из-за чего точность оценки base_error падаеты


# In[201]:


for i in range(M.shape[0]): # [Предпологаем сильные ошибки]
    for j in range(M.shape[1]):
        if M[i, j] < N[i, j]:
            M[i, j] = -1
            N[i, j] = -1


# In[203]:


print(f'#base error is {base_err}')


MC = get_mc(M, N, base_err)


# In[207]:


idx = np.where((MC[:,0]!=MC[:,1]))[0] # Так можно?
MC = MC[idx]

# In[4]


P = do_mcmc(nIter, n_threads=1, model=1, show_each=10)
