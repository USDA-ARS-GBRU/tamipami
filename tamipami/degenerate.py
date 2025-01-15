from typing import List
import logging
from typing import Iterable
# -*- coding: utf-8 -*-
from itertools import product

import numpy as np
import textdistance
import treelib
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy

from tamipami.config import config

codes = config["codes"]

def expand_degenerate_strings(degenerate_strings: list[str]) -> list[str]:
    """
    Expand a list of degenerate sequences in IUPAC notation into all possible sequences.

    Args:
        degenerate_strings (list[str]): A list of IUPAC degenerate sequences.

    Returns:
        list[str]: A list of all possible sequences derived from the degenerate sequences.
    """

    def expand_string(degenerate_string):
        expanded = [codes[char] if char in codes else {char} for char in degenerate_string]
        return ["".join(chars) for chars in product(*expanded)]

    expanded_strings = set()
    for degenerate_string in degenerate_strings:
        expanded_strings.update(expand_string(degenerate_string))

    return list(expanded_strings)

def degenerate_representation(position_to_chars: list[set]) -> str:
    # Create a dictionary for fast lookup
    chars_to_code = {frozenset(options): code for code, options in codes.items()}
    
    # Create a list to hold potential representatives
    minimal_representations = []

    for chars in map(frozenset, position_to_chars):
        if len(chars) == 1:
            minimal_representations.append(next(iter(chars)))
        else:
            # Use the dictionary for fast lookup
            code = chars_to_code.get(frozenset(chars), chars)
            minimal_representations.append(code)
    return "".join(minimal_representations)


def unique_characters(expanded_strings: set | list[str]) -> list[set]:
    """Given a set of strings returns the unique characters at each position
    Args:
        expanded_string: A set of strings of equal length
    returns:
        A list with the unique characters at each position
    """
    length = len(next(iter(expanded_strings)))
    position_to_chars = [set() for _ in range(length)]
    # Collect unique characters for each position
    for string in expanded_strings:
        for index, char in enumerate(string):
            position_to_chars[index].add(char)
    return position_to_chars



def degenerate_representation(position_to_chars: list[set]) -> list:
    # Create a list to hold potential representatives
    minimal_representations = []

    for chars in position_to_chars:
        if len(chars) == 1:
            minimal_representations.append(chars.pop())
        else:
            # Find the best degenerate character representation
            found = False
            for code, options in codes.items():
                if options == chars:
                    minimal_representations.append(code)
                    found = True
                    break
            if not found:
                minimal_representations.append(chars)
    return "".join(minimal_representations)


def gen_all_combo(minimal_representations: list[set | list]) -> set[str]:
    """Takes a list of nucleotides at  each position and returns the full set of expanded strings
    Args:
        minimal_representations  a list of possible nucleotides at each position e.g. [{'A', 'C'}, {'G'}, {'T', 'C'}]
    Returns:
        a list of all sequences in expanded notation
    """
    all_possible_representations = set()
    for combination in product(
        *[
            list(chars) if isinstance(chars, set) else [chars]
            for chars in minimal_representations
        ]
    ):
        all_possible_representations.add("".join(combination))
    return all_possible_representations


def create_degenerate(subset: set[str]) -> str | None:
    """Create a degenerate sequence if exactly all members of the set can be represented by one degenerate sequence, else returns None
    Args:
        subset: takes a set of kmer strings
    Returns:
        a single degenerate sequence string or None if no string can be formed

    """
    position_to_chars = unique_characters(subset)
    if subset == gen_all_combo(position_to_chars):
        return degenerate_representation(position_to_chars)
    else:
        return None


### We  now have code that will return the degenerate code if all members of the group can be represented by one code.
### The next setep is to apply a subset selection algorithm that can  tst possible groups efficiently


def calculate_hamming_distance(strings: List[str]) -> np.array:
    """take a list of strings and return a condensed hamming distance matrix and the strings as labels
    Args:
        strings: a list of equal length sequence strings
    Returns:
        a condensed scipy distance matrix as a numpy array
    """
    if not all(len(s) == len(strings[0]) for s in strings):
        raise ValueError("All strings must be of equal length")
    # prepare 2 dimensional array M x N (M entries (3) with N dimensions (1))
    transformed_strings = np.array(strings).reshape(-1, 1)
    # calculate condensed distance matrix by wrapping the hamming distance function
    distance_matrix = pdist(
        transformed_strings, lambda x, y: textdistance.hamming(x[0], y[0])
    )
    # get square matrix
    return distance_matrix


def create_tree(labels: list[str], distance_matrix: np.array) -> treelib.Tree:
    """Convert the output of scipy.cluster.hierarchy.to_tree() into a treelib.Tree object.
    Args:
        distance_matrix: A condensed scipy distance matrix
        labels: A list of strings used as tags for the nodes in the same order as the nodes list.
    Returns:
        A treelib.Tree object representing the hierarchical structure.
    """

    Z = scipy.cluster.hierarchy.linkage(distance_matrix, method="ward")
    # Convert linkage matrix to a tree structure
    root = scipy.cluster.hierarchy.to_tree(Z, rd=False)

    def convert_scipy_tree_to_treelib(root, labels):
        """
        Convert a tree from scipy.cluster.hierarchy.to_tree to a treelib.Tree.

        Args:
            root: The root node of the scipy hierarchical tree.
            labels: A list of labels for the leaf nodes.

        Returns:
            A treelib.Tree object representing the hierarchical structure.
        """
        tree = treelib.Tree()
        node_dict = {}

        def add_children(scipy_node, parent_id=None):
            node_id = str(scipy_node.id)
            if scipy_node.is_leaf():
                tag = labels[scipy_node.id] if scipy_node.id < len(labels) else f"Leaf {scipy_node.id}"
            else:
                tag = f"Node {scipy_node.id}"

            tree.create_node(tag=tag, identifier=node_id, parent=parent_id)
            node_dict[scipy_node.id] = node_id

            if scipy_node.left is not None:
                add_children(scipy_node.left, node_id)
            if scipy_node.right is not None:
                add_children(scipy_node.right, node_id)

        add_children(root)
        return tree
    return convert_scipy_tree_to_treelib(root, labels)

def find_degenerates(tree: treelib.Tree) -> list[str]:
    """
    Traverse a tree using breadth-first search to identify and process leaf nodes.

    This function performs a breadth-first traversal of a given treelib.Tree
    object, processing each leaf node to determine if a degenerate sequence
    can be created from the leaf tags. If a degenerate sequence is found,
    it is added to the results list.

    Parameters:
    - tree: A treelib.Tree object representing the hierarchical structure to be traversed.

    Returns:
    - A list of degenerate sequence strings generated from the leaf nodes of the tree.
    """
    # Temporary list to store results from create_degenerate
    results = []
    # Expand the tree in breadth-first order to list tree node numbers as strings
    nodes = (node for node in tree.expand_tree(mode=treelib.Tree.WIDTH))
    node_skip_list = set()
    for node in nodes:
        if node not in node_skip_list:
            # Call create_degenerate with leaf tags
            seqs: set[str] = {item.tag for item in tree.leaves(node)}
            result = create_degenerate(seqs)
            if result is not None:  # If it returns a string, save it
                results.append(result)
                node_skip_list.update(tree.subtree(node).nodes.keys())
    return results




def seqs_to_degenerates(seqs: list[str]) -> list[str]:
    """Takes a list of equal length sequences and finds a minimal set of degenerate codes that represents them.
    Args:
        seqs: A list of sequences where each sequence is a string of equal length. The sequences should only contain valid nucleotide characters (e.g., A, T, C, G).
    Returns:
        a list of degenerate sequences representing the input list
    """
    if not seqs:
        return []

    try:
        logging.debug("Calculating Hamming distance matrix.")
        distance_matrix = calculate_hamming_distance(seqs)
        logging.debug("Creating tree from distance matrix.")
        tree = create_tree(seqs, distance_matrix)
        logging.debug("Finding degenerate sequences.")
        return find_degenerates(tree)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        return []


