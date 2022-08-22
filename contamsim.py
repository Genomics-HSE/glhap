#!/usr/bin/env python
# coding: utf-8



#simple contamination read create





import os
import sys


# In[4]:


def make_fastq(*args, size=100000, err_base=0.03, del_base=0.00, ins_base=0.00, total_reads=100000):
    assert(len(args)%2 == 0)
    
    genomes_fnames = []
    proportion = []
    
    for i, x in enumerate(args):
        if i % 2 == 0:
            genomes_fnames.append(x)
        else:
            proportion.append(float(x))
    
    cat_str = 'samtools cat'
    os.system('rm contaminants.fa')
    os.system('touch contaminants.fa')
    for i, genome_file in enumerate(genomes_fnames):
        os.system(f'cp {genome_file} ./genome_{i}.fa')
        if i != 0:
            os.system(f'cat {genome_file} >> contaminants.fa')
        os.system(f'''sed -i '' "1s/.*/>chrM/" genome_{i}.fa''')
        fname =  f'genome_{i}.fa'
        command_line =         f'simlord        -rr {fname}        -pi {ins_base}        -pd {del_base}        -ps {err_base}        -fl 100        -n {int(total_reads*proportion[i])}        genome_{i}'
        os.system(command_line)
        os.system(f'samtools view -b genome_{i}.sam > genome_{i}.bam')
        os.system(f'rm genome_{i}.sam genome_{i}.fastq')
        cat_str += f' genome_{i}.bam'
    
    os.system(f'{cat_str} | samtools sort > simulated_data.bam')
        
        
    
        
        
    print('finish')
    


# In[6]:

#print(sys.argv)
make_fastq(*sys.argv[1:])


# In[ ]:





# In[ ]:




