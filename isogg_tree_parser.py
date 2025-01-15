import json
from anytree import Node, NodeMixin
from anytree.exporter import JsonExporter
from pprint import pprint
import csv

snp_csv_filename = 'Phylotree/SNP_Index.csv'
tree_in_fn = 'isogg_tree.json'

class YNode(NodeMixin):
    """
    A class representing a node in a Y-chromosome haplogroup tree.
    Attributes:
        name (str): The name of the node.
        snps (list): A list of SNPs (Single Nucleotide Polymorphisms) associated with the node.
        insertion (list): A list of insertions associated with the node.
        deletion (list): A list of deletions associated with the node.
        children_ids (list): A list of IDs of the node's children.
        parent_id (int): The ID of the node's parent.
        lh (int): A placeholder attribute, default is 0.
    Methods:
        __repr__(): Returns a string representation of the node.
        __str__(): Returns the name of the node as a string.
    """
    def __init__(self, name, snps=[],insertion=[], deletion = [],snps_back=[], children_ids=[], parent_id=0, parent=None, children=None):
        self.name = name
        self.snps = snps
        self.insertion = insertion
        self.deletion = deletion
        self.children_ids = children_ids
        self.parent_id = parent_id
        self.snps_back = snps_back
        # self.parent = parent
        self.lh = 0
        # if children:  # set children only if given
             # self.children = children
    
    def __repr__(self):
             return str(self.name) +" "+ str(self.lh)
    def __str__(self):
        return self.name
    

def make_snp_dict(filename):
    nucleotides = set(['A', 'T', 'G', 'C'])
    snp_csv = filename
    snp_dict = {}
    with open(snp_csv) as f:
        csvfile = csv.DictReader(f, delimiter=',')
        for line in csvfile:
            haplogroup = line['Subgroup Name']
            snp_dict.setdefault(haplogroup, [])
            if (len(line['Mutation Info'])==4 and
                line['Mutation Info'][0].upper() in nucleotides and
                line['Mutation Info'][3].upper() in nucleotides):  # проверка что SNP
                snp_dict[haplogroup].append([line['Mutation Info'][0].upper(),
                                             line['Mutation Info'][-1].upper(),
                                             int(line['Build 38 Number'])])

            
    return(snp_dict)
    
    

# Какой мы хотим формат данных
# {name: str, variants: [list of variants in format (old_pos_new)]}, children:[ids of childrens]


def make_tree(json_file: dict, snps: dict):
    haplogroup_dict = dict()
    for haplogroup in json_file:
        if haplogroup not in haplogroup_dict:
            haplogroup_dict[haplogroup] = YNode(name=haplogroup)
        for child in json_file[haplogroup]:
            # print(haplogroup, child, haplogroup_dict[haplogroup])
            haplogroup_dict[child] = YNode(name=child, parent=haplogroup_dict[haplogroup])
            haplogroup_dict[child].parent = haplogroup_dict[haplogroup]
    for haplogroup in haplogroup_dict:
        if haplogroup in snps:
            haplogroup_dict[haplogroup].snps = snps[haplogroup]
    # print(haplogroup_dict)
    return haplogroup_dict["ROOT (Y-Chromosome 'Adam')"]

snp_dict = make_snp_dict(snp_csv_filename)
with open(tree_in_fn) as f:
        tree = json.load(f)
root = make_tree(tree, snp_dict)
exporter = JsonExporter(indent=2, sort_keys=True)
with open('isogg.json', 'w') as f:
        f.write(exporter.export(root))
