import subprocess
import pandas as pd
import io
import os
import glhap_app
import random
import matplotlib.pyplot as plt
import numpy as np
import sys


fn = sys.argv[1]
iterations = int(sys.argv[2])

print('FILE IS', fn,  flush=True)


def get_mean_depth(fn):
    stats = pd.read_csv(io.StringIO(subprocess.check_output(['samtools', 'coverage', fn]).decode()), sep="\t").iloc[0]
    if (stats.iloc[0] != 'chrM'):
        return 0
    mean_depth = stats.iloc[6]
    return mean_depth


def downsample(fn, depth, original_depth=None):
    while True:
        if original_depth == None:
            original_depth = get_mean_depth(fn)
        prop = depth/original_depth
        os.system(f'samtools view -b -s {prop + random.randint(1,10000)} {fn} > {fn[:-4]}.cov_{depth}.bam')
        if get_mean_depth(f'{fn[:-4]}.cov_{depth}.bam') != 0:
            break
    return 0


def haplogrep_result(vcf_fn):
    os.system(f'haplogrep classify --input {vcf_fn} --output {vcf_fn[:-7]}.haplogrep.txt --format vcf > /dev/null')
    with open(f'{vcf_fn[:-7]}.haplogrep.txt') as f:
        s =f.readlines()
    s = s[1].replace('"','').split('\t')[1].replace(' ','')
    
    return s


def glhap_results(vcf_fn):
    s = glhap_app.glhap('refchrm.fa', 'array/array.json', f'{vcf_fn}' )
    s = s.splitlines()[1].split('|')[0].replace(' ','')
    return s


def get_glhap_accuracy(fn, depth, iterations=100):
    os.system(f'bcftools mpileup -Oz -f refchrm.fa {fn} > {fn[:-4]}.mpileup.vcf.gz 2>/dev/null')
    os.system(f'bcftools index {fn[:-4]}.mpileup.vcf.gz 2>/dev/null')
    os.system(f'bcftools call -m -Oz {fn[:-4]}.mpileup.vcf.gz > {fn[:-4]}.call.vcf.gz 2>/dev/null')
    
    true_haplogroup = haplogrep_result(f'{fn[:-4]}.call.vcf.gz')
    coinc = 0
    original_depth = get_mean_depth(fn)
    for i in range(iterations):
        downsample(fn, depth, original_depth)
        fn_down = f'{fn[:-4]}.cov_{depth}.bam'
        os.system(f'bcftools mpileup -Oz -f refchrm.fa {fn_down} > {fn_down[:-4]}.mpileup.vcf.gz 2>/dev/null')
        os.system(f'bcftools index {fn_down[:-4]}.mpileup.vcf.gz 2>/dev/null')
        hap = glhap_results(f'{fn_down[:-4]}.mpileup.vcf.gz')
        if hap==true_haplogroup:
            coinc +=1
        os.system(f'rm {fn_down[:-4]}.mpileup.vcf.gz  {fn_down} ')
    os.system(f'rm {fn[:-4]}.call.vcf.gz {fn[:-4]}.mpileup.vcf.gz')
    return coinc / iterations


def get_haplogrep_accuracy(fn, depth, iterations=100):
    os.system(f'bcftools mpileup -Oz -f refchrm.fa {fn} > {fn[:-4]}.mpileup.vcf.gz 2>/dev/null')
    os.system(f'bcftools index {fn[:-4]}.mpileup.vcf.gz 2>/dev/null')
    os.system(f'bcftools call -m -Oz {fn[:-4]}.mpileup.vcf.gz > {fn[:-4]}.call.vcf.gz 2>/dev/null')
    
    true_haplogroup = haplogrep_result(f'{fn[:-4]}.call.vcf.gz')
    coinc = 0
    original_depth = get_mean_depth(fn)
    for i in range(iterations):
        downsample(fn, depth, original_depth)
        fn_down = f'{fn[:-4]}.cov_{depth}.bam'
        os.system(f'bcftools mpileup -Oz -f refchrm.fa {fn_down} > {fn_down[:-4]}.mpileup.vcf.gz 2>/dev/null')
        os.system(f'bcftools index {fn_down[:-4]}.mpileup.vcf.gz 2>/dev/null')
        os.system(f'bcftools call -m -Oz {fn_down[:-4]}.mpileup.vcf.gz > {fn_down[:-4]}.call.vcf.gz 2>/dev/null')
        hap = haplogrep_result(f'{fn_down[:-4]}.call.vcf.gz')
        if hap==true_haplogroup:
            coinc +=1
        os.system(f'rm {fn_down[:-4]}.mpileup.vcf.gz {fn_down[:-4]}.call.vcf.gz {fn_down}')
    os.system(f'rm {fn[:-4]}.call.vcf.gz {fn[:-4]}.mpileup.vcf.gz')
    return coinc / iterations
        
        
        
    


# In[158]:

depths = [0.1 + 0.1 * i for i in range(0,40)]

# In[156]:

def get_accuracy_graph(fn, iterations=10):
    ys_glhap = []
    ys_haplogrep = []
    for d in depths:
        ys_glhap.append(get_glhap_accuracy(fn, d, iterations))
        ys_haplogrep.append(get_haplogrep_accuracy(fn, d, iterations))
    # plt.plot(depths, ys_glhap, label = 'Glhap')
    # plt.plot(depths, ys_haplogrep, label = 'Haplogrep')
    # plt.legend()
    # plt.show()
    np.savetxt(f'glhap/{fn.split('.')[0]}.glhap.txt', [depths, ys_glhap])
    np.savetxt(f'haplogrep/{fn.split('.')[0]}.haplogrep.txt', [depths, ys_haplogrep])
    


# In[157]:


get_accuracy_graph(fn, iterations)


os.system(f"rm  {fn[:-4]}*cov* {fn[:-4]}*mpileup* {fn[:-4]}*call*")
print('FINISH')
