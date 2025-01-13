import json
from anytree import Node, NodeMixin
from anytree.exporter import JsonExporter
from pprint import pprint
import csv

snp_csv_filename = 'SNP_Index.csv'
tree_in_fn = 'yfull_ytree_extracted_2024-12-27.json'

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
    snp_csv = filename
    snp_dict = {}
    with open(snp_csv) as f:
        csvfile = csv.DictReader(f, delimiter=',')
        for line in csvfile:
            names = line['Alternate Names'].split()
            if '' in names:
                names.remove('')
            
            names.append(line['Name'])
            for name in names:
                snp_dict[name] = [line['Build 38 Number'],  line['Mutation Info']]
    return(snp_dict)
    
    
snp_dict = make_snp_dict(snp_csv_filename)

# Какой мы хотим формат данных
# {name: str, variants: [list of variants in format (old_pos_new)]}, children:[ids of childrens]


def make_tree_dict(json_file: dict, parent=-1):
    haplogroup_list = []
    idx = json_file['id'] if json_file['id']!= '' else 'root'
    haplogroup_list.append({'name': json_file['id'] if json_file['id']!= '' else 'root' , 'snps_index':  json_file['snps'], "parent":parent})
    if 'children' in json_file:
        for child in json_file['children']:
            haplogroup_list += make_tree_dict(child, parent=idx)
    return(haplogroup_list)

def expand_snps(haplogroup_dict: dict, snps_dict: dict):
    # TODO: Добавить обработку инделов
    haplogroup_dict.setdefault("snps", [])
    nucleotides = set(['A', 'T', 'G', 'C'])
    for snp_codes in haplogroup_dict['snps_index'].split(','):
        for snp_code in snp_codes.strip().split('/'):
            # print(snp_code)
            if snp_code in snps_dict:
                if len(snps_dict[snp_code][1])==4 and snps_dict[snp_code][1][0] in nucleotides and snps_dict[snp_code][1][3] in nucleotides:
                    haplogroup_dict['snps'].append([snps_dict[snp_code][1][0], snps_dict[snp_code][1][3], int(snps_dict[snp_code][0])])
                break
                    


def make_anytree(tree_json_filename, snp_csv_filename):
    with open(tree_in_fn) as f:
        tree = json.load(f)
    haplo_dict = make_tree_dict(tree)   
    snp_dict = make_snp_dict(snp_csv_filename)
    anytree_dict = {}
    for x in haplo_dict:
        expand_snps(x, snp_dict)
        anytree_dict[x['name']] = YNode(name=x['name'], parent_id=x['parent'], snps=x['snps'])
    for x in anytree_dict:
        if x!='root':
            anytree_dict[x].parent = anytree_dict[anytree_dict[x].parent_id]
    root = anytree_dict['root']
    exporter = JsonExporter(indent=2, sort_keys=True)
    with open('yFull.json', 'w') as f:
        f.write(exporter.export(root))

    return(anytree_dict['root'])

make_anytree(tree_in_fn, snp_csv_filename)
    
