#!/usr/bin/env python
# coding: utf-8

# In[1]:


from anytree import NodeMixin, RenderTree
import numpy as np
import pysam
from pysam import VariantFile, FastaFile
from anytree import find_by_attr, PreOrderIter
import json


# In[2]:


class Node(NodeMixin):
    def __init__(self, name, snps=[],insertion=[],deletion = [], parent=None, children=None):
        self.name = name
        self.snps = snps
        self.insertion = insertion
        self.deletion = deletion
        self.parent = parent
        self.lh = 0
        if children:  # set children only if given
             self.children = children
    
    def __repr__(self):
             return self.name +" "+ str(self.lh)
    def __str__(self):
        return self.name


# In[3]:


def get_snp(node):
    '''
    Returns snps fron node: Node
    
    Paramters:
    node : Node
    node in phylogenetic tree
    
    Returns:
    list of lists of [old_nucleotide,new_nucleotid, position]
    '''
    snp = list()
    atgc = set(['A','T','G','C','a','t','g','c'])
    for i in range(2,len(node)):
        if node[i][0] in atgc and node[i][-1] in atgc:
            snp.append([node[i][0],node[i][-1],int(node[i][1:-1])])
    return snp


# In[4]:


def get_insertion(node):
    ins = set()
    for i in range(2,len(node)):
        if '.' in node[i]:
            dot_pos = node[i].find('.')
            pos = int(node[i][1:dot_pos])
#             size = ''.join(k for k in node[i][dot_pos+1:] if  k.isdigit())
#             insert = ''.join(k for k in node[i][dot_pos+1:] if  k.isalpha())
            ins.add(pos)
    return ins


# In[5]:


def get_deletion(node):
    deletion = set()
    for i in range(2,len(node)):
        if node[i][-1]=='d':
            if node[i][0].isalpha():
                deletion.add(int(node[i][1:-1]))
            else:
                dash_pos = node[i].find('-')
                deletion.add(int(node[i][1:dash_pos]))
    return deletion


# In[6]:


def make_tree(tree,node,pos=0):    
    '''
    Initializate tree structure of phylogenetic tree
    Parameters:
    tree: list
    list from phylotree-parse
    node: Node
    root node
    pos: int
    unnecessary parameter. position of current node in tree list
    
    '''

    
    posit = pos + 1
    i = posit
    while i<len(tree) and tree[i][0] >= tree[pos][0]+1:
        if tree[i][0] == tree[pos][0]+1:
#             print(node.name,tree[i][1])
            snps = get_snp(tree[i])
            insertion = get_insertion(tree[i])
            tmp = Node(tree[i][1],snps,insertion,parent=node)
            
            make_tree(tree,tmp,i)

        i += 1        


# In[7]:


def get_log_monozygous(bcf: VariantFile):
    '''
    PL scores for bcf
    
    '''
    gls = np.full((16569,4),-1)
    for rec in bcf.fetch():
        pos = rec.pos
        pls = rec.samples.values()[0]['PL']
        alt = rec.alleles
        k = 0
        s = 2
        for i in range(len(alt)):
            if alt[i] == 'A':
                gls[pos-1][0] = pls[k]
            elif alt[i] == 'T':
                gls[pos-1][1] = pls[k]
            elif alt[i] == 'G':
                gls[pos-1][2] = pls[k]
            elif alt[i] == 'C':
                gls[pos-1][3] = pls[k]
            elif alt[i] == '<*>':
                for j in range(4):
                    if gls[pos-1][j] == -1:
                        gls[pos-1][j] = pls[k]
            k += s
            s += 1
    return -gls/10


# In[8]:


def call_likelihood(gls,node,ref,insertions,deletions,lh=0):
    '''
    Calculates likelihood for node considering parent node likelihood.
    
    Parameters:
    gls : np.array
    array of snp likelihoods
    node: Node
    haplogroup node
    ref: FastaFile
    reference genome
    lh: float
    genotype likelihood for parent haplogroup
    
    Returns
    
    lh : float
    likelihood of haplogroup
    '''
    snps = node.snps
    for snp in snps:
#       snp = [old, new, pos]
        pos = snp[2]-1
        if snp[1].capitalize() == 'A':
            lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,0]
#             print(1,snp[1].capitalize())
        if snp[1].capitalize() == 'T':
            lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,1]
#             print(2,snp[1].capitalize())
        if snp[1].capitalize() == 'G':
            lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,2]
#             print(3,snp[1].capitalize())
        if snp[1].capitalize() == 'C':
            lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,3]
#             print(4,snp[1].capitalize())
        
    
        for ins in insertions:
            if ins[0] in node.insertion:
#           ins = [pos,lh]
                pos = ins[0]-1
                lh = lh - calculate_pl(gls,ref,pos) + lh[1]

        for delt in deletions:
            if ins[0] in node.deletion:
#           delt = [pos,lh]
                pos = delt[0]-1
                lh = lh - calculate_pl(gls,ref,pos) + delt[1]
    return lh
        
    
    


# In[9]:


def calculate_likelihood(vcf, ref):
    '''
    Calculates genotype likelihood against reference genome
    
    
    Parameters:
    -----------
    vcf : VariantFile
    vcf/bcf file
    
    ref: FastaFile
    reference genome
    --
    Returns
    -------
    lh : float
    likelihood against reference
    '''
    lh = 0
    ref = ref.fetch('chrM')
    gls = get_log_monozygous(vcf)
#     gls[gls<0]=-10**6
    for i in range(len(ref)):
        if ref[i].capitalize() == 'A':
            lh += gls[i,0]
        if ref[i].capitalize() == 'T':
            lh += gls[i,1]
        if ref[i].capitalize() == 'G':
            lh += gls[i,2]
        if ref[i].capitalize() == 'C':
            lh += gls[i,3]
    return lh


# In[10]:


def prunung(node,ref,gls,deletions,insertions):
    '''
    Calculates genotype likelihood for each haplogroup in the tree
    Parameters:
    -----------
    node: Node
    root node in the tree structure
    ref: FastaFile
    reference file
    gls: np.array
    matrix of snp pl scores
    '''
    if node.parent == None:
        node.lh = ref_lh
    else:
        node.lh = call_likelihood(gls,node,ref,insertions,deletions,node.parent.lh)
    for i in node.children:
        prunung(i,ref,gls,deletions,insertions)


# In[11]:


def calculate_pl(gls, ref,pos):
    '''
    Calculate pl score of certain position on reference
    
    Parameters:
    -----------
    gls : np.array
    matrix of snp pl scores
    ref: FastaFile
    reference
    pos: int
    position
    '''
    lh = 0
    ref = ref.fetch('chrM')
#     gls = get_log_monozygous(vcf)
    i = pos
    if ref[i].capitalize() == 'A':
        lh += gls[i,0]
    if ref[i].capitalize() == 'T':
        lh += gls[i,1]
    if ref[i].capitalize() == 'G':
        lh += gls[i,2]
    if ref[i].capitalize() == 'C':
        lh += gls[i,3]
    return lh


# In[12]:


import argparse
parser = argparse.ArgumentParser(description='Calculation of haplogroups')
parser.add_argument('tree', type=str, help='haplogroup tree')
parser.add_argument('ref', type=str, help='reference fasta')
parser.add_argument('vcf', type=str, help='vcf file of mtdna')
args = parser.parse_args()


# In[13]:



with open(args.tree) as f:
    d = json.load(f) # d - это список python
#     print(d)
for i in d:
    i[0] += 1
    d[0][0]=0

bcf_in = VariantFile(args.vcf) 
ref = FastaFile(args.ref)


# In[14]:




# with open('PhyloTree.org-parser/array.json') as f:
#     d = json.load(f) # d - это список python
# #     print(d)
# for i in d:
#     i[0] += 1
#     d[0][0]=0

# bcf_in = VariantFile("PhyloTree.org-parser/out4.vcf") 
# ref = FastaFile('ref.fa')

# # print(d)


# In[15]:


a = Node(d[0][1],[])
make_tree(d,a,0)


# In[16]:


ref_lh = calculate_likelihood(bcf_in,ref)


# In[17]:


# INSERTIONS
insertions = []
for rec in bcf_in.fetch():
    if (rec.info['INDEL']== True):
        if len(rec.ref)== 1 and len(rec.alts[0]) == 2:
            insertions.append([rec.pos + 1, -rec.samples.values()[0]['PL'][-1]/10])


# In[18]:


# Deletions
deletions = []
for rec in bcf_in.fetch():
    if (rec.info['INDEL']== True):
        if len(rec.ref)== 2 and rec.ref[-1]!='N' and len(rec.alts[0]) == 1:
            deletions.append([rec.pos + 1, -rec.samples.values()[0]['PL'][-1]/10])


# In[19]:


gls = get_log_monozygous(bcf_in)
prunung(a,ref,gls,insertions,deletions)


# In[20]:


key = lambda x: x.lh


# In[21]:


S = list()
for i in PreOrderIter(a):
    S.append(i)


# In[22]:


S.sort(key = key,reverse=True)


# In[23]:


for i in range(10):
    print(S[i].name+",lh = ",S[i].lh)


# In[ ]:





# In[ ]:




