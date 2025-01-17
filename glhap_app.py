#!/usr/bin/env python
# coding: utf-8

import time  # measure time of make_tree_from_yFull_json function
from anytree.importer import JsonImporter
import pickle
import json
from anytree.importer import DictImporter
from anytree import find_by_attr, PreOrderIter
from pysam import VariantFile, FastaFile
import pysam
import numpy as np
from anytree import NodeMixin, RenderTree


class Node(NodeMixin):
    def __init__(
        self,
        name,
        snps=[],
        snps_back=[],
        insertion=[],
        deletion=[],
        parent=None,
        children=None,
    ):
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


def make_tree_from_yFull_json(json_fn):
    """
    Make tree from json file

    Parameters:
    json_fn : str
    path to json file

    Returns:
    root : Node
    root node of the tree
    """
    try:
        with open(json_fn) as f:
            json_text = f.read()
    except OSError:
        return "No such tree file"

    importer = JsonImporter()
    root = importer.import_(json_text)
    return root


def get_snp(node: Node):
    """
    Returns snps from node: Node

    Paramters:
    node : Node
    node in phylogenetic tree

    Returns:
    list of lists of [old_nucleotide, new_nucleotid, position]
    """

    snp = list()
    atgc = set(["A", "T", "G", "C", "a", "t", "g", "c"])
    for i in range(2, len(node)):
        record = node[i].replace("(", "").replace(")", "")
        if record[0] in atgc and record[-1] in atgc:
            snp.append([record[0], record[-1], int(record[1:-1])])
    return snp


def get_snp_back(node):
    """
    Returns snps fron node: Node

    Paramters:
    node : Node
    node in phylogenetic tree

    Returns:
    list of lists of [old_nucleotide,new_nucleotid, position]
    """

    snp = list()
    atgc = set(["A", "T", "G", "C", "a", "t", "g", "c"])
    for i in range(2, len(node)):
        record = node[i].replace("(", "").replace(
            ")", "")  # We agree with all variants
        if record[0] in atgc and record[-1] == "!" and record[-2] in atgc:
            snp.append([record[0], record[-2], int(record[1:-2])])
    return snp


def get_insertion(node):
    ins = set()
    for i in range(2, len(node)):
        if "." in node[i]:
            dot_pos = node[i].find(".")
            pos = int(node[i][1:dot_pos])
            #             size = ''.join(k for k in node[i][dot_pos+1:] if  k.isdigit())
            #             insert = ''.join(k for k in node[i][dot_pos+1:] if  k.isalpha())
            ins.add(pos)
    return ins


def get_deletion(node):
    deletion = set()
    for i in range(2, len(node)):
        if node[i][-1] == "d":
            if node[i][0].isalpha():
                deletion.add(int(node[i][1:-1]))
            else:
                dash_pos = node[i].find("-")
                deletion.add(int(node[i][1:dash_pos]))
    return deletion


def make_tree(tree, root):
    """
    Initializate tree structure of phylogenetic tree
    Parameters:
    tree: list
    list from phylotree-parse
    node: Node
    root node
    pos: int
    unnecessary parameter. position of current node in tree list
    """
    queue = [[0, root]]
    while len(queue) != 0:
        pos, node = queue.pop(0)
        i = pos + 1
        while i < len(tree) and tree[i][0] > tree[pos][0]:
            if tree[i][0] == tree[pos][0] + 1:
                snps = get_snp(tree[i])
                snps_back = get_snp_back(tree[i])
                insertion = get_insertion(tree[i])

                tmp = Node(tree[i][1], snps, snps_back, insertion, parent=node)
                queue.append([i, tmp])
            i += 1
    return root


def get_log_monozygous(bcf: VariantFile, chrom="chrM"):
    """
    Calculate the log-scaled PL scores for a given chromosome from a BCF file.

    Parameters:
    bcf (VariantFile): The input BCF file containing variant data.
    chrom (str): The chromosome to fetch data for. Default is 'chrM'.

    Returns:
    numpy.ndarray: A 2D array of log-scaled PL scores for the specified chromosome.

    Notes:
    - For 'chrM', the length of the chromosome is assumed to be 16569.
    - For 'chrY', the start and end positions are 2781490 and 56887900, respectively.
    - The PL scores are divided by 10 and negated to obtain the log-scaled values.
    - The function handles alleles 'A', 'T', 'G', 'C', and '<*>'.
    """
    if chrom == "chrM":
        N = 16569
    elif chrom == "chrY":
        N = 57227415
    gls = np.full((N, 4), -1)  # TODO: CHANGE FOR SOMETHING BETTER
    for i, rec in enumerate(bcf.fetch(chrom)):
        # if i == 5: stop for testing
        # return - gls / 10
        pos = rec.pos - 1
        pls = rec.samples.values()[0]["PL"]
        alt = rec.alleles
        k = 0
        s = 2
        for i in range(len(alt)):
            if alt[i] == "A":
                gls[pos][0] = pls[k]
            elif alt[i] == "T":
                gls[pos][1] = pls[k]
            elif alt[i] == "G":
                gls[pos][2] = pls[k]
            elif alt[i] == "C":
                gls[pos][3] = pls[k]
            elif alt[i] == "<*>":
                for j in range(4):
                    if gls[pos][j] == -1:
                        gls[pos][j] = pls[k]
            k += s
            s += 1
    return -gls / 10


def call_likelihood(gls, node, ref, insertions, deletions, chrom, lh=0):
    """
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
    """
    if chrom == "chrM":
        N = 16569
    elif chrom == "chrY":
        N = 57227415

    snps = node.snps
    snps_back = node.snps_back

    base_dict = {"A": 0, "T": 1, "G": 2, "C": 3}
    for snp in snps:
        #       snp = [old, new, pos]
        pos = snp[2]-1

        lh = (
            lh
            - gls[pos, base_dict[snp[0].capitalize()]]
            + gls[pos, base_dict[snp[1].capitalize()]]
        )

    for snp in snps_back:
        #       snp = [old, new, pos]
        pos = snp[2] - 1
        lh = (
            lh
            - gls[pos, base_dict[snp[0].capitalize()]]
            + gls[pos, base_dict[snp[1].capitalize()]]
        )

    # TODO implementation to include indels.
    return lh


def calculate_likelihood(ref, gls, chrom="chrM"):
    """
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
    """

    lh = 0
    ref = ref.fetch(chrom)

    if chrom == "chrM":
        N = 16569
    elif chrom == "chrY":
        N = 57227415
    # gls = get_log_monozygous(vcf)

    for i in range(N):
        if ref[i].capitalize() == "A":
            lh += gls[i, 0]
        elif ref[i].capitalize() == "T":
            lh += gls[i, 1]
        elif ref[i].capitalize() == "G":
            lh += gls[i, 2]
        elif ref[i].capitalize() == "C":
            lh += gls[i, 3]
    return lh


def pruning(node, ref, gls, deletions, insertions, ref_lh, chrom):
    """
    Calculates genotype likelihood for each haplogroup in the tree
    Parameters:
    -----------
    node: Node
    root node in the tree structure
    ref: FastaFile
    reference file
    gls: np.array
    matrix of snp pl scores
    """

    # global global_it
    # print("global_it = ", global_it)
    # global_it += 1

    if node.parent == None:
        node.lh = ref_lh
    else:
        node.lh = call_likelihood(
            gls, node, ref, insertions, deletions, chrom, node.parent.lh
        )
    for i in node.children:
        pruning(i, ref, gls, deletions, insertions, ref_lh, chrom)


def calculate_pl(gls, ref, pos, chrom="chrM"):
    """
    Calculate pl score of certain position on reference

    Parameters:
    -----------
    gls : np.array
    matrix of snp pl scores
    ref: FastaFile
    reference
    pos: int
    position
    """
    lh = 0
    base_dict = {"A": 0, "T": 1, "G": 2, "C": 3}
    lh = gls[pos, base_dict[ref[pos].capitalize()]]
    return lh


def glhap(ref, tree, chrom, bcf_in):
    if chrom == "chrM":
        try:
            with open(tree) as f:
                d = json.load(f)
        except OSError:
            return "No such tree file"

        for i in d:
            i[0] += 1
            d[0][0] = 0
        a = Node(d[0][1], [])
        a = make_tree(d, a, 0)
    elif chrom == "chrY":
        a = make_tree_from_yFull_json(tree)
    else:
        print("No such chromosome")
        return "No such chromosome"

    print("Tree is loaded")

    try:
        ref = FastaFile(ref)
    except OSError:
        return "No such reference file"

    try:
        bcf_in = VariantFile(bcf_in)
    except OSError:
        return "No such vcf.gz file"

    print("Reference and bcf file are opened")

    gls = get_log_monozygous(bcf_in, chrom)
    print("Matrix of genotype likelihoods is calculated")

    ref_lh = calculate_likelihood(ref, gls, chrom)  # No need anymore
    print("Reference likelihood is calculated")

    deletions = []  # Имеет смысл добавить вставки удаления корректно
    insertions = []
    if chrom == "chrM":
        N = 16569
        for rec in bcf_in.fetch(chrom):
            if rec.info["INDEL"] == True:
                if len(rec.ref) == 1 and len(rec.alts[0]) == 2:
                    insertions.append(
                        [rec.pos + 1, -rec.samples.values()[0]["PL"][-1] / 10]
                    )

        for rec in bcf_in.fetch(chrom):
            if rec.info["INDEL"] == True:
                if len(rec.ref) == 2 and rec.ref[-1] != "N" and len(rec.alts[0]) == 1:
                    deletions.append(
                        [rec.pos + 1, -rec.samples.values()[0]["PL"][-1] / 10]
                    )

    pruning(a, ref.fetch(chrom), gls, deletions, insertions, ref_lh, chrom)
    print("haplogroup likelihoods are calculated")

    def key(x):
        return x.lh

    S = list()

    for i in PreOrderIter(a):
        S.append(i)

    S.sort(key=key, reverse=True)
    j = 0
    i = 0
    out = ""
    str1 = ["Haplogroup", "PL score"]
    out += f"{str1[0]:<15}|\t{str1[1]:<15}\n"
    while i != 10:
        if S[i].name != "-":
            out += f"{S[i].name:<15}|\t{S[i].lh:<15}\n"
            j += 1
        i += 1
    print(out)

    return out


# ref_fname = "reference/Homo_sapiens.GRCh38.dna.chromosome.Y.fa"
# tree_fname = "isogg.json"
# vcf_fname = "Y_inputfiles/vcf/in.vcf.gz"

out = glhap(
    ref_fname, tree_fname, "chrY", vcf_fname
)
print(out)

start = time.time()
tree = make_tree_from_yFull_json("isogg.json")
end = time.time()
elapsed_time = end - start
print(elapsed_time)


bcf_in = VariantFile("Y_inputfiles/vcf/in1.vcf.gz")
ref = FastaFile("reference/Homo_sapiens.GRCh38.dna.chromosome.Y.fa")

start = time.time()
gls = get_log_monozygous(bcf_in, "chrY")
end = time.time()
elapsed_time = end - start
print(elapsed_time)

ref_lh = calculate_likelihood(ref, gls, chrom)
chrom
insertions = []
deletions = []
global_it = 0
pruning(tree, ref.fetch(chrom), gls, deletions, insertions, ref_lh, chrom)

len(list(PreOrderIter(tree)))

S = list()

for i in PreOrderIter(tree):
    S.append(i)

S.sort(key=lambda x: x.lh, reverse=True)
S[0]

# ref_str = ref.fetch('chrY')
# ref_str[7948322]
# (gls != 0).sum()
# ref_lh
# tree.children[1].snps
# gls[7948322]

# gls[2781478]
