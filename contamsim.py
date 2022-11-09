#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#simple contamination read create


# In[29]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os
import sys
from tqdm import tqdm
import numpy as np


# In[30]:


import argparse


# In[57]:


parser = argparse.ArgumentParser(description='Call out my name. Adepti Xiao. I will be here')
parser.add_argument('--fnames', nargs='*', help='genomes filenames')
parser.add_argument('--proportions', nargs='*',type=float , help='proportions of contaminants')
parser.add_argument('--coverage', default=20, type=float, help='coverage')
parser.add_argument('--err_base',default=0.01, help='coverage')
parser.add_argument('--del_base',default=0., help='coverage')
parser.add_argument('--ins_base',default=0., help='coverage')
parser.add_argument('--output',default='simulated_data.bam',type=str, help='name of final file')


# In[58]:


def make_fastq(genomes_fnames, proportion, err_base, del_base, ins_base, coverage, output):
    
    assert np.sum(proportion) == 1
    
    cat_str = 'samtools cat'
    os.system('rm contaminants.fa')
    os.system('touch contaminants.fa')
    for i, genome_file in enumerate(genomes_fnames):
        os.system(f'cp {genome_file} ./genome_{i}.fa')
        if i != 0:
            os.system(f'cat {genome_file} >> contaminants.fa')
        os.system(f'''sed -i '' "1s/.*/>chrM/" genome_{i}.fa''')
        fname =  f'genome_{i}.fa'
        command_line =         f'simlord        -rr {fname}        -pi {ins_base}        -pd {del_base}        -ps {err_base}        -fl 100        -c {proportion[i] * coverage}        genome_{i}'
        os.system(command_line)
        os.system(f'samtools view -b genome_{i}.sam > genome_{i}.bam')
        os.system(f'rm genome_{i}.sam genome_{i}.fa genome_{i}.fastq')
        cat_str += f' genome_{i}.bam'
    
    os.system(f'{cat_str} | samtools sort > {output}')
        
        
    
        
        
    print('finish')
    


# In[59]:


args = parser.parse_args(sys.argv[1:])
#"--fnames fasta/A.fasta fasta/L2a1a.fasta --proportions 0.7 0.3 --coverage 3".split()


# In[60]:


args.proportions


# In[61]:


make_fastq(args.fnames, args.proportions, args.err_base, args.del_base, args.ins_base, args.coverage, args.output)


# In[ ]:





# In[ ]:




