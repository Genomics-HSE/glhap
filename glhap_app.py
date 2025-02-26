#!/usr/bin/env python
# coding: utf-8

import time  # measure time of make_tree_from_yFull_json function
from anytree.importer import JsonImporter
import pickle
import json
from anytree.importer import DictImporter
from anytree.search import find
from anytree import find_by_attr, PreOrderIter, Walker
from pysam import VariantFile, FastaFile
import pysam
import numpy as np
from anytree import NodeMixin, RenderTree
from pprint import pprint
from collections import deque
from typing import Any, List, Set, Dict, Tuple
from anytree import Node, RenderTree
from anytree.exporter import DotExporter
from IPython.display import Image
import os

class PhylogeneticNode(NodeMixin):
    """
    Represents a node in a phylogenetic tree.

    Attributes:
        name (str): Unique identifier or name for the node.
        snps (list): List of forward single nucleotide polymorphisms associated with the node.
        backward_snps (list): List of reverse/compensatory single nucleotide polymorphisms.
        insertions (list): List of insertion mutations.
        deletions (list): List of deletion mutations.
        parent (PhylogeneticNode, optional): Reference to the parent node.
        children (list): List of child nodes.
        likelihood (float): Likelihood score of the node (default is 0).
    """

    def __init__(
        self,
        name: str,
        snps: list = None,
        backward_snps: list = None,
        insertions: list = None,
        deletions: list = None,
        parent: "PhylogeneticNode" = None,
        children: list = None,
    ):
        """
        Initializes a new instance of PhylogeneticNode.

        Args:
            name (str): The name or identifier of the node.
            snps (list, optional): A list of forward SNPs. Defaults to an empty list.
            backward_snps (list, optional): A list of backward SNPs. Defaults to an empty list.
            insertions (list, optional): A list of insertion mutations. Defaults to an empty list.
            deletions (list, optional): A list of deletion mutations. Defaults to an empty list.
            parent (PhylogeneticNode, optional): Reference to the parent node. Defaults to None.
            children (list, optional): List of child nodes. Defaults to an empty list.
        """
        self.name = name
        self.snps = snps if snps is not None else []
        self.backward_snps = backward_snps if backward_snps is not None else []
        self.insertions = insertions if insertions is not None else []
        self.deletions = deletions if deletions is not None else []
        self.parent = parent
        self.likelihood = 0
        self.children = children if children is not None else []

    def __repr__(self) -> str:
        """
        Returns an unambiguous string representation of the node.

        Returns:
            str: A string containing the node's name, likelihood, and SNPs.
        """
        return f"{self.name} {self.likelihood} {self.snps}"

    def __str__(self) -> str:
        """
        Returns a human-readable string representation of the node.

        Returns:
            str: The node's name.
        """
        return self.name



def import_yfull_tree(json_file_path: str) -> Any:
    """
    Create a tree structure from a YFull JSON file.

    This function reads the JSON file specified by `json_file_path` and uses a 
    JsonImporter to parse its content into a tree structure. It returns the 
    root node of the tree.

    Args:
        json_file_path (str): The path to the YFull JSON file.

    Returns:
       Node: The root node of the tree.

    Raises:
        FileNotFoundError: If the file cannot be found or opened.
        Exception: If an error occurs during JSON import.
    """
    try:
        with open(json_file_path, "r", encoding="utf-8") as file:
            json_content = file.read()
    except OSError as error:
        raise FileNotFoundError(f"Unable to open the JSON file: {json_file_path}") from error

    importer = JsonImporter()
    try:
        root_node = importer.import_(json_content)
    except Exception as error:
        raise Exception("An error occurred while importing the JSON content into a tree structure.") from error

    return root_node


from typing import List, Dict, Any

def extract_snps(node: Node, chrom: str = 'chrY', reference: str = 'rCRS',  rcrs_dict: Dict[int, str] = None) -> List[List]:
    """
    Extract SNP (Single Nucleotide Polymorphism) information from a phylogenetic node.

    The function assumes that the node is a sequence where SNP records start from the third element.
    Each SNP record is expected to be formatted with parentheses around the position and nucleotide change,
    e.g. "(A123T)". The function cleans the record, extracts the original nucleotide, position, and
    the new nucleotide, and then replaces the original nucleotide with an alternate value from rcrs_dict
    if available.

    Args:
        node Node: A phylogenetic node where SNP records are stored from index 2 onward.
        rcrs_dict (Dict[int, str]): A dictionary mapping positions to alternative nucleotides.
    
    Returns:
        List[List]: A list of SNP records, each represented as a list:
                    [original_or_replaced_old_nucleotide, new_nucleotide, position].
    """
    snps = []
    valid_nucleotides = {"A", "T", "G", "C", "a", "t", "g", "c"}

    # Process SNP records starting from index 2
    for record in node[2:]:
        cleaned_record = record.replace("(", "").replace(")", "")
        old_nucleotide = cleaned_record[0]
        new_nucleotide = cleaned_record[-1]
        
        # Validate that both nucleotides are valid
        if old_nucleotide in valid_nucleotides and new_nucleotide in valid_nucleotides:
            try:
                position = int(cleaned_record[1:-1])
            except ValueError:
                # Skip records with invalid position formatting
                continue
            
            # Replace old nucleotide if a mapping exists in rcrs_dict
            if reference == 'rCRS':
                updated_old = rcrs_dict.get(position, old_nucleotide)
            snps.append([updated_old, new_nucleotide, position])

    return snps


def extract_backward_snps(node: Any) -> List[List]:
    """

    Args:
        node (Any): A phylogenetic node containing SNP records.

    Returns:
        List[List]: A list of backward SNP records, where each record is a list of the form:
                    [original_nucleotide, backward_nucleotide, position].
    """
    backward_snps = []
    valid_nucleotides = {"A", "T", "G", "C", "a", "t", "g", "c"}

    for record in node[2:]:
        cleaned_record = record.replace("(", "").replace(")", "")
        if (cleaned_record[0] in valid_nucleotides and 
            cleaned_record[-1] == "!" and 
            cleaned_record[-2] in valid_nucleotides):
            try:
                position = int(cleaned_record[1:-2])
            except ValueError:
                # Skip record if position conversion fails
                continue
            backward_snps.append([cleaned_record[0], cleaned_record[-2], position])
    
    return backward_snps


def get_insertion(node):
    """
    Extract insertion positions from a list of mutation records in a phylogenetic node.

    This function scans through the mutation records starting from index 2, looking for records
    that indicate an insertion mutation. An insertion is identified by the presence of a dot (".")
    within the record. The function then extracts the insertion position (an integer found between
    the first character and the dot) and returns a set of unique positions.

    Args:
        mutation_records (List[str]): A list of mutation records (strings) from a phylogenetic node.
                                      Insertion records should contain a dot character.

    Returns:
        Set[int]: A set of insertion positions extracted from the mutation records.
    """
    insertion_positions: Set[int] = set()
    
    # Process mutation records starting from index 2
    for record in node[2:]:
        if "." in record:
            dot_index = record.find(".")
            try:
                position = int(record[1:dot_index])
                insertion_positions.add(position)
            except ValueError:
                # Skip record if position conversion fails
                continue

    return insertion_positions


def get_deletion(node):
    deletion_positions: Set[int] = set()

    # Process mutation records starting from index 2
    for record in node[2:]:
        if record.endswith("d"):
            try:
                if record[0].isalpha():
                    # Format: e.g. "A123d" => position = 123
                    pos_str = record[1:-1]
                else:
                    # Format: e.g. "123-XYZd" => position is taken from characters before the dash
                    dash_index = record.find("-")
                    if dash_index == -1:
                        continue  # Skip if dash is not found in the expected format
                    pos_str = record[1:dash_index]
                
                deletion_positions.add(int(pos_str))
            except ValueError:
                # Skip record if position conversion fails
                continue

    return deletion_positions


def build_mt_phylogenetic_tree(tree_filename: str) -> PhylogeneticNode:
    """
    Build a phylogenetic tree from a JSON file containing tree data.

    This function loads tree data from a JSON file specified by `tree_filename`. The JSON data 
    is expected to be a list, where each element is a record representing a node in the tree. 
    The first element of each record indicates the node's level in the tree (which will be 
    incremented by one for proper tree construction, except for the root which is set to 0), 
    and the second element is the node's identifier (or name).

    In addition, a mapping file 'rsrs_to_rcrs.txt' is loaded to map positions to reference 
    nucleotides. This mapping is used when extracting SNP (Single Nucleotide Polymorphism) data 
    for each node using helper functions (`get_snp`, `get_snp_back`, and `get_insertion`).

    The tree is constructed via a breadth-first traversal of the list, establishing parent-child 
    relationships based on node levels.

    Args:
        tree_filename (str): Path to the JSON file containing the tree data.

    Returns:
        Node: The root node of the constructed phylogenetic tree.

    Raises:
        FileNotFoundError: If the JSON file or the mapping file 'rsrs_to_rcrs.txt' is not found.
        Exception: If an error occurs during file parsing.
    """
    # Load tree data from the JSON file.
    try:
        with open(tree_filename, "r", encoding="utf-8") as file:
            tree_data: List[List[Any]] = json.load(file)
    except OSError as error:
        raise FileNotFoundError(f"Tree file not found: {tree_filename}") from error

    # Load the RSRS-to-rCRS mapping.
    rcrs_mapping: Dict[int, str] = {}
    try:
        with open("rsrs_to_rcrs.txt", "r", encoding="utf-8") as file:
            for line in file:
                parts = line.split()
                rcrs_mapping[int(parts[0])] = parts[2]
    except OSError as error:
        raise FileNotFoundError("Mapping file 'rsrs_to_rcrs.txt' not found") from error


    # Adjust the level information for each record:
    # Increment all node levels by 1, then explicitly set the root level to 0.
    for record in tree_data:
        record[0] += 1
    tree_data[0][0] = 0

    # Initialize the root node using the identifier from the first record.
    root = PhylogeneticNode(tree_data[0][1], [])
    
    # Use a deque for efficient FIFO queue operations during breadth-first construction.
    queue = deque()
    queue.append((0, root))  # Each element is a tuple: (index in tree_data, corresponding Node)

    while queue:
        current_index, current_node = queue.popleft()
        next_index = current_index + 1

        # Process child nodes: nodes with a level greater than the current node.
        while next_index < len(tree_data) and tree_data[next_index][0] > tree_data[current_index][0]:
            # Only direct children have a level exactly one greater.
            if tree_data[next_index][0] == tree_data[current_index][0] + 1:
                # Extract mutation data for the child node using helper functions.
                snps = extract_snps(tree_data[next_index], rcrs_dict=rcrs_mapping)
                backward_snps = extract_backward_snps(tree_data[next_index])
                insertion = get_insertion(tree_data[next_index])
                deletion = get_deletion(tree_data[next_index])
                
                # Create the child node with the extracted mutation data.
                child_node = PhylogeneticNode(tree_data[next_index][1], snps, backward_snps, [], parent=current_node)
                queue.append((next_index, child_node))
            next_index += 1

    return root


def compute_log_genotype_likelihood(
    bcf: VariantFile, chrom: str = "chrM"
) -> Tuple[np.ndarray, Dict[int, List[float]], Dict[int, List[float]]]:
    """
    Calculate log-scaled genotype likelihoods (PL scores) for a specified chromosome from a BCF file.

    This function fetches variant records from the provided BCF file for the given chromosome,
    computes log-scaled PL scores for each variant position, and collects insertion and deletion
    information in separate dictionaries.

    Parameters:
        bcf (VariantFile): Input BCF file containing variant data.
        chrom (str): Chromosome to process. Supported values are "chrM" and "chrY". Default is "chrM".

    Returns:
        Tuple[np.ndarray, Dict[int, List[float]], Dict[int, List[float]]]:
            - A 2D numpy array (shape: [N, 4]) of log-scaled PL scores for the chromosome, where N is the
              chromosome length.
            - A dictionary mapping positions (int) to lists of two float values for insertion events.
            - A dictionary mapping positions (int) to lists of two float values for deletion events.

    Notes:
        - For "chrM", the chromosome length is assumed to be 16569.
        - For "chrY", the chromosome length is set to 57227415.
        - PL scores are divided by 10 and negated to obtain log-scaled values.
        - The function handles alleles "A", "T", "G", "C", and "<*>".
        - For INDEL records, the first and last PL scores (after scaling) are stored.
        - The ordering of PL scores for SNP records assumes that the first score is at index 0 and
          subsequent indices are determined by an increasing offset.
    """
    # Set chromosome length based on specified chromosome.
    if chrom == "chrM":
        chrom_length = 16569
    elif chrom == "chrY":
        chrom_length = 57227415
    else:
        raise ValueError(f"Unsupported chromosome: {chrom}")
    
    # Dictionaries to store insertion and deletion data.
    insertions: Dict[int, List[float]] = {}
    deletions: Dict[int, List[float]] = {}

    # Initialize an array for genotype likelihoods with a default value (-1).
    log_gl_scores = np.full((chrom_length, 4), -1.0, dtype=float)

    # Iterate over variant records for the specified chromosome.
    for rec in bcf.fetch(chrom):
        # Assume a single sample per record.
        sample_data = list(rec.samples.values())[0]
        current_log_pl_scores = sample_data["PL"]

        # Process INDEL records separately.
        if rec.info.get("INDEL", False):
            # Check if it is an insertion (reference allele is shorter than alternate allele).
            if len(rec.ref) < len(rec.alts[0]):
                insertions[rec.pos] = [
                    -log_gl_scores[0] / 10.0,
                    -log_gl_scores[-1] / 10.0,
                ]
            else:
                deletions[rec.pos] = [
                    -log_gl_scores[0] / 10.0,
                    -log_gl_scores[-1] / 10.0,
                ]
            continue

        # Process SNP records.
        pos = rec.pos - 1  # Adjust from 1-based to 0-based indexing.
        alleles = rec.alleles

        # The PL scores ordering is assumed to follow a specific pattern:
        # Starting at index 0, the gap for the next allele increases by 1 each time.
        allele_pl_index = 0
        increment = 2
        for allele in alleles:
            current_pl = current_log_pl_scores[allele_pl_index]
            if allele == "A":
                log_gl_scores[pos, 0] = current_pl
            elif allele == "T":
                log_gl_scores[pos, 1] = current_pl
            elif allele == "G":
                log_gl_scores[pos, 2] = current_pl
            elif allele == "C":
                log_gl_scores[pos, 3] = current_pl
            elif allele == "<*>":
                # For allele "<*>", fill any missing score (-1) in the current row.
                for j in range(4):
                    if log_gl_scores[pos, j] == -1:
                        log_gl_scores[pos, j] = current_pl
            allele_pl_index += increment
            increment += 1

    # Return the negated and scaled genotype likelihoods along with INDEL dictionaries.
    return -log_gl_scores / 10.0, insertions, deletions


import numpy as np
from typing import Dict, List, Any

def call_likelihood(
    log_gl_scores: np.ndarray,
    node: Any,
    insertions: Dict[int, List[float]],
    deletions: Dict[int, List[float]],
    parent_likelihood: float = 0.0
) -> float:
    """
    Calculate the updated likelihood for a phylogenetic node based on SNPs and indels.

    This function updates the genotype likelihood (lh) of a haplogroup node by adjusting it
    according to SNP changes, backward SNP changes, and indel events (insertions and deletions).
    
    For each SNP record (of the form [old_allele, new_allele, position]), the likelihood is
    updated by subtracting the log-scaled likelihood of the old allele and adding that of the new allele.
    A similar adjustment is applied to backward SNPs. Additionally, for insertion and deletion events
    present in the node (if they exist in the provided dictionaries), the likelihood is adjusted by
    subtracting the first score and adding the second score.

    Parameters:
        log_gl_scores (np.ndarray): A 2D array of logarithmed genotype likelihoods, where each row corresponds to a genomic
                          position and each column corresponds to one of the nucleotides in the order A, T, G, C.
        node (Any): A node in the phylogenetic tree with the following attributes:
                    - snps: List of SNP records [old_allele, new_allele, position].
                    - backward_snps: List of backward SNP records [old_allele, new_allele, position].
                    - insertion: List of insertion positions.
                    - deletion: List of deletion positions.
        insertions (Dict[int, List[float]]): A dictionary mapping positions (int) to a list of two float values
                                               representing log-scaled likelihoods for insertion events.
        deletions (Dict[int, List[float]]): A dictionary mapping positions (int) to a list of two float values
                                              representing log-scaled likelihoods for deletion events.
        chrom (str): The chromosome identifier (currently unused, reserved for future use).
        lh (float, optional): The initial likelihood (e.g., from the parent haplogroup). Defaults to 0.0.

    Returns:
        float: The updated likelihood for the haplogroup node.
    """
    # Mapping from nucleotide to index in the gls array.
    nucleotide_to_index = {"A": 0, "T": 1, "G": 2, "C": 3}
    haplogroup_likelihood = parent_likelihood
    # Process forward SNPs: each SNP record is [old_allele, new_allele, position]
    for old_allele, new_allele, pos in node.snps:
        pos_index = pos - 1  # Convert 1-based position to 0-based index
        haplogroup_likelihood += log_gl_scores[pos_index][nucleotide_to_index[new_allele.upper()]] - log_gl_scores[pos_index][nucleotide_to_index[old_allele.upper()]]

    # Process backward SNPs similarly
    for old_allele, new_allele, pos in node.backward_snps:
        pos_index = pos - 1  # Convert 1-based position to 0-based index
        haplogroup_likelihood += log_gl_scores[pos_index, nucleotide_to_index[new_allele.upper()]] - log_gl_scores[pos_index, nucleotide_to_index[old_allele.upper()]]

    # Adjust likelihood based on insertion events present in the node
    for ins_pos in node.insertions:
        if ins_pos in insertions:
            # Subtract the first score and add the second score for insertions.
            haplogroup_likelihood += -insertions[ins_pos][0] + insertions[ins_pos][1]

    # Adjust likelihood based on deletion events present in the node
    for del_pos in node.deletions:
        if del_pos in deletions:
            haplogroup_likelihood += -deletions[del_pos][0] + deletions[del_pos][1]
    return haplogroup_likelihood


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


def calculate_haplogroups_likelihoods(
    node: Any,
    gls: np.ndarray,
    deletions: Dict[int, List[float]],
    insertions: Dict[int, List[float]],
    ref_lh: float,
) -> None:
    """
    Recursively update the genotype likelihood for each haplogroup node in the phylogenetic tree.

    For the root node (i.e. with no parent), the likelihood is set to the reference likelihood (ref_lh).
    For non-root nodes, the likelihood is computed by calling `call_likelihood` with the parent's likelihood.
    The function then recursively propagates the likelihood update to all descendant nodes.

    Parameters:
        node (Any): Current node in the haplogroup tree. Expected to have attributes:
                    - parent: Parent node (or None if root)
                    - children: List of child nodes
                    - lh: Current likelihood value
                    - snps, backward_snps, insertion, deletion: Mutation data for likelihood adjustment.
        gls (np.ndarray): 2D numpy array of SNP PL scores (rows: positions, columns: A, T, G, C).
        deletions (Dict[int, List[float]]): Dictionary mapping positions to deletion likelihood adjustments.
        insertions (Dict[int, List[float]]): Dictionary mapping positions to insertion likelihood adjustments.
        ref_lh (float): Reference likelihood to be used for the root node.
        chrom (str): Chromosome identifier (e.g. "chrM").
    """
    if node.parent is None:
        node.lh = ref_lh
    else:
        node.lh = call_likelihood(gls, node, insertions, deletions, node.parent.lh)
    
    for child in node.children:
        calculate_haplogroups_likelihoods(child, gls, deletions, insertions, ref_lh)



def calc_genomic_distance(node1: Node, node2: Node) -> int:
    """
    Calculate the genomic distance between two nodes in a phylogenetic tree.

    The genomic distance is defined as the total number of mutations (SNPs and indels)
    along the paths connecting the two nodes. This function uses a Walker to obtain the
    paths between node1 and node2 and aggregates the mutations (from SNPs, backward SNPs,
    insertions, and deletions) encountered along both the upward and downward traversal.

    Parameters:
        node1 (Node): The first node in the phylogenetic tree.
        node2 (Node): The second node in the phylogenetic tree.

    Returns:
        int: The total number of mutations (SNPs and indels) between the two nodes.
    """
    walker = Walker().walk(node1, node2)
    mutations: List[Any] = []
    
    # Process both segments of the walk (e.g., from node1 to LCA and from LCA to node2)
    for path in (walker[0], walker[2]):
        for node in path:
            mutations.extend(node.snps)
            mutations.extend(node.snps_back)
            mutations.extend(node.insertions)
            mutations.extend(node.deletions)
            
    return len(mutations)

def glhap(tree_file: str, chrom: str, bcf_path: str, verbose: bool = False) -> str:
    """
    Calculate haplogroup genotype likelihoods and return a formatted summary of top haplogroups.

    The function loads a phylogenetic tree from a file (using different loaders based on the
    chromosome), opens a BCF file, calculates a matrix of genotype likelihoods, and computes
    haplogroup likelihoods via a pruning algorithm. If verbose is True, detailed progress
    messages are printed.

    Parameters:
        tree_file (str): Path to the tree file.
        chrom (str): Chromosome identifier ("chrM" or "chrY").
        bcf_path (str): Path to the BCF (or VCF) file.
        verbose (bool): If True, prints detailed status messages.

    Returns:
        str: A formatted string containing the top 10 haplogroups (excluding placeholder nodes)
             and their PL scores.
    """
    # Load tree based on the chromosome.
    if chrom == "chrM":
        tree_root = build_mt_phylogenetic_tree(tree_file)
    elif chrom == "chrY":
        tree_root = import_yfull_tree(tree_file)
    else:
        msg = f"Unsupported chromosome: {chrom}"
        if verbose:
            print(msg)
        return msg

    if verbose:
        print("Tree is loaded.")

    # Open the BCF file.
    try:
        bcf_file = VariantFile(bcf_path)
    except OSError:
        msg = f"Unable to open BCF file: {bcf_path}"
        if verbose:
            print(msg)
        return msg

    if verbose:
        print("BCF file is opened.")

    # Compute the matrix of genotype likelihoods and associated indel data.
    gls, ins_mapping, del_mapping = compute_log_genotype_likelihood(bcf_file, chrom)
    if verbose:
        print(f"Genotype likelihood matrix calculated with shape {gls.shape}.")

    # Set the reference likelihood. (Currently hardcoded to 0.)
    ref_lh = 0
    if verbose:
        print("Reference likelihood is set to 0.")

    # Use empty lists for insertions and deletions (pruning uses these lists).
    insertions_list = []
    deletions_list = []
    calculate_haplogroups_likelihoods(tree_root, gls, deletions_list, insertions_list, ref_lh)
    if verbose:
        print("Haplogroup likelihoods are calculated using pruning.")

    # Collect nodes in pre-order traversal.
    nodes = list(PreOrderIter(tree_root))

    # Sort nodes by likelihood in descending order.
    sorted_nodes = sorted(nodes, key=lambda node: node.lh, reverse=True)

    # Build a formatted table with header and top 10 haplogroups (excluding nodes with name "-").
    header = f"{'Haplogroup':<25}|\t{'inverse of log-gl score':<15}"
    output_lines = [header]
    count = 0
    for node in sorted_nodes:
        if node.name != "-":
            output_lines.append(f"{node.name:<25}|\t{node.lh:<15.2f}")
            count += 1
            if count == 10:
                break

    if verbose:
        print("Top haplogroups summary generated.")
    print("\n".join(output_lines))
    likelihood_ratio = 2*(sorted_nodes[0].lh-sorted_nodes[1].lh)/np.log(2)
    print(f'Top 1 likelihood ratio test:{likelihood_ratio:<5.2f}')
    return "\n".join(output_lines)




# def get_this():
#     print('file\thaplogroup')
#     S_list = list()
#     for i in range(1,11):
#         bcf_in = VariantFile(f"ydna/in{i}.vcf.gz")
#         ref = FastaFile("reference/Homo_sapiens.GRCh38.dna.chromosome.Y.fa")
#         tree = make_tree_from_yFull_json("isogg.json")
#         gls, insertions, deletions = get_log_monozygous(bcf_in, "chrY")
#         pruning(tree, ref.fetch(chrom), gls, deletions, insertions, ref_lh=0, chrom='chrY')
#         S = list()
#         for j in PreOrderIter(tree):
#             S.append(j)
#         S.sort(key=lambda x: x.lh, reverse=True)
#         print(f'in{i}\t{S[0].name}')
#         S_list.append(S)
#     return S_list
    
    
    
# def get_confusion(vcf_fn, haplogroup_name):
#     # reference = FastaFile(reference_fn)
#     vcf = VariantFile(vcf_fn)
#     haplo_node = find(tree, lambda node : node.name == haplogroup_name)
#     node = haplo_node
#     snps = list()
#     snps_dict = dict()
#     for snp in node.snps:
#         snps_dict[snp[2]]=snp[1]
#     while node.parent is not None:
#         node = node.parent
#         for snp in node.snps:
#             snps_dict[snp[2]]=snp[1]
#     print(len(snps_dict))
#     m = 0
#     n = 0
#     vcf_fetch = vcf.fetch('chrY')
#     for rec in vcf_fetch:
#         if rec.alts is not None and len(rec.alts[0])==1:
#             # print(rec.alts)
#             if rec.pos in snps_dict:
#                 m += 1
#             else:
#                 n += 1
#         # else:
#             # if rec.pos in snps_dict:
#                 # n += 1
#     return m, n
                
# glhap('/Users/nikita/Desktop/HSE/genomic/glyhap/Phylotree/mtDNA.json', 'chrM', 'mtdna/in1.vcf.gz', verbose=True)


glhap('isogg.json', 'chrY', 'ydna/in1.cov_0.1.mpileup.vcf.gz', verbose=True)