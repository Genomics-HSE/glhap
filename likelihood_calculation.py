from pysam import VariantFile, FastaFile
import numpy as np

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
    return gls




def call_likelihood(gls,node,ref,lh=0):
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
            lh = lh - calculate_pl(gls,ref, pos) + gls[pos,0]
        if snp[1].capitalize() == 'T':
            lh = lh - calculate_pl(gls,ref, pos) + gls[pos,1]
        if snp[1].capitalize() == 'G':
            lh = lh - calculate_pl(gls,ref, pos) + gls[pos,2]
        if snp[1].capitalize() == 'C':
            lh = lh - calculate_pl(gls,ref, pos) + gls[pos,3]
        
    
#         for ins in insertions:
#             if ins[0] in node.insertion:
#             ins = [pos,lh]
#                 pos = ins[0]-1
#                 lh = lh - calculate_pl(gls,ref,pos) + lh[1]

#         for delt in deletions:
#             if ins[0] in node.deletion:
#             delt = [pos,lh]
#                 pos = delt[0]-1
#                 lh = lh - calculate_pl(gls,ref,pos) + delt[1]
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
    
    
def prunung(node, ref, ref_lh, gls):
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
        node.lh = call_likelihood(gls,node,ref,node.parent.lh)
    for i in node.children:
        prunung(i,ref, ref_lh, gls)
        
        
def calculate_pl(gls, ref, pos):
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