from anytree import NodeMixin, RenderTree
from anytree import find_by_attr, PreOrderIter
import json
import numpy as np


#Next class represents nodes from phylotree with all their mutations
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
        elif node[i][0] in atgc and node[i][-1]=='!' and node[i][-1] in atgc:
            snp.append([node[i][0],node[i][-2],int(node[i][1:-2])])
    return snp


def get_insertion(node):
    '''
    Currently unused
    '''
    ins = set()
    for i in range(2,len(node)):
        if '.' in node[i] and 'd' not in node[i] and '!' not in node[i] and '(' not in node[i] and ')' not in node[i]:
            dot_pos = node[i].find('.')
            pos = int(node[i][0:dot_pos])
            size = ''.join(k for k in node[i][dot_pos+1:] if  k.isdigit())
            insert = ''.join(k for k in node[i][dot_pos+1:] if  k.isalpha())
            ins.add(pos)
    return ins


def get_deletion(node):
    '''
    Currently unused
    '''
    deletion = set()
    for i in range(2,len(node)):
        if node[i][-1]=='d':
            if node[i][0].isalpha():
                # print(node)
                deletion.add(int(node[i][1:-1]))
            else:
                dash_pos = node[i].find('-')
                # deletion.add(int(node[i][1:dash_pos]))
    return deletion


def make_tree_(tree,node,pos=0):    
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
            if tree[i][2] == 'reserved':
                i += 1
                continue
            # print(node.name,tree[i][1])
            snps = get_snp(tree[i])
            insertion = get_insertion(tree[i])
            # insertion = set()
            deletion = get_deletion(tree[i])
            # deletion = set()
            tmp = Node(tree[i][1], snps, insertion, deletion, parent=node)
            
            make_tree_(tree,tmp,i)

        i += 1        

    
    
    
    
def make_tree(json_file):
    with open(json_file) as f:
        d = json.load(f) # d - это список python
    for i in d:
        i[0] += 1
    d[0][0]=0
    a = Node(d[0][1],[])
    make_tree_(d,a,0)
    return a
      