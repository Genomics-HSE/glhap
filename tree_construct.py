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
    
    

    