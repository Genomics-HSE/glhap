#!/usr/bin/env python
# coding: utf-8


from anytree import NodeMixin, RenderTree
import numpy as np
import pysam
from pysam import VariantFile, FastaFile
from anytree import find_by_attr, PreOrderIter
from anytree.importer import DictImporter
import json
import pickle


class Node(NodeMixin):
    def __init__(self, name, snps=[], snps_back=[], insertion=[],deletion = [], parent=None, children=None):
        self.name = name
        self.snps = snps
        self.snps_back = snps_back
        self.insertion = insertion
        self.deletion = deletion
        self.parent = parent
        self.lh = 0
        if children:  # set children only if given
             self.children = children
    
    def __repr__(self):
             return self.name + " " + str(self.lh) + " " + str(self.snps)
    def __str__(self):
        return self.name

    
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
    for i in range(2, len(node)):
        record = node[i].replace('(', '').replace(')', '')
        if record[0] in atgc and record[-1] in atgc:
            snp.append([record[0], record[-1], int(record[1:-1])])
    return snp


def get_snp_back(node):
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
    for i in range(2, len(node)):
        record = node[i].replace('(', '').replace(')', '')
        if record[0] in atgc and record[-1]=='!' and record[-2] in atgc:
            snp.append([record[0], record[-2], int(record[1:-2])])
    return snp



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


def make_tree(tree,root,pos=0):    
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
    
    queue = [[0, root]]
    while len(queue) != 0:
        pos, node = queue.pop(0)
        i = pos+1
        while i<len(tree) and tree[i][0] > tree[pos][0]:
            if tree[i][0] == tree[pos][0]+1:
    #             print(node.name,tree[i][1])
                snps = get_snp(tree[i])
                snps_back = get_snp_back(tree[i])
                insertion = get_insertion(tree[i])
                
                tmp = Node(tree[i][1],snps, snps_back,insertion,parent=node)
                queue.append([i, tmp])
                # print(tree[i][1])
            i += 1  
    # print(root.size)
    # print(len(tree))
    return root


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
    return - gls/10


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
    snps_back = node.snps_back
    
    base_dict= {'A':0, 'T':1, 'G':2, 'C':3}
    for snp in snps:
#       snp = [old, new, pos]
        pos = snp[2]-1
        # if snp[1]x.capitalize() == 'A':
            # lh = lh - calculate_pl(gls, ref, pos)+ gls[pos,0]
#             print(1,snp[1].capitalize())
        # if snp[1].capitalize() == 'T':
            # lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,1]
#             print(2,snp[1].capitalize())
        # if snp[1].capitalize() == 'G':
            # lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,2]
#             print(3,snp[1].capitalize())
        # if snp[1].capitalize() == 'C':
            # lh = lh - calculate_pl(gls,ref, pos)+ gls[pos,3]
#             print(4,snp[1].capitalize())

        lh = lh - calculate_pl(gls,ref, pos)+ gls[pos, base_dict[snp[1].capitalize()]]


    for snp in snps_back:
#       snp = [old, new, pos]
        pos = snp[2]-1
        lh = lh - gls[pos, base_dict[snp[0].capitalize()]]+ gls[pos, base_dict[snp[1].capitalize()]]

#TODO implementation to include indels.
    return lh


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


def pruning(node,ref,gls,deletions,insertions, ref_lh):
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
        pruning(i,ref,gls,deletions,insertions, ref_lh)


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


def glhap(ref, tree, bcf_in):
    try:
        with open(tree) as f:
            d = json.load(f)
    except OSError:
        return("No such tree file")

    for i in d:
        i[0] += 1
        d[0][0]=0
    a = Node(d[0][1],[])
    a = make_tree(d,a,0)
    
    try:
        ref = FastaFile(ref)
    except OSError:
        return("No such reference file")
    
    try:
        bcf_in = VariantFile(bcf_in) 
    except OSError:
        return("No such vcf.gz file")
    
    ref_lh = calculate_likelihood(bcf_in,ref)
    


    insertions = []
    for rec in bcf_in.fetch():
        if (rec.info['INDEL']== True):
            if len(rec.ref)== 1 and len(rec.alts[0]) == 2:
                insertions.append([rec.pos + 1, -rec.samples.values()[0]['PL'][-1]/10])



    deletions = []
    for rec in bcf_in.fetch():
        if (rec.info['INDEL']== True):
            if len(rec.ref)== 2 and rec.ref[-1]!='N' and len(rec.alts[0]) == 1:
                deletions.append([rec.pos + 1, -rec.samples.values()[0]['PL'][-1]/10])
    gls = get_log_monozygous(bcf_in)
    pruning(a,ref,gls,insertions,deletions, ref_lh)
    key = lambda x: x.lh
    S = list()

    for i in PreOrderIter(a):
        S.append(i)

    S.sort(key = key,reverse=True)
    j = 0
    i= 0 
    out=''
    str1 = ['Haplogroup', 'PL score']
    out+=f"{str1[0]:<15}|\t{str1[1]:<15}\n"
    while i!=10:
        if S[i].name != '-':
            out+=f"{S[i].name:<15}|\t{S[i].lh:<15}\n"
            j+=1
        i+=1
    print(out)
    
    return out