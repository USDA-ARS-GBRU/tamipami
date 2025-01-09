# -*- coding: utf-8 -*-
from itertools import product

import numpy as np
import textdistance
import treelib
from scipy.spatial.distance import pdist
import scipy.cluster.hierarchy

from .config import config

codes = config["codes"]


def expand_degenerate_strings(degenerate_strings: list[str]) -> list[str]:
    """Takes a list of degenerate sequences in IUPAC notation and expands it to a list of the sequences that make up that degenerate list
    Args:
        degenerate_strings:  a list of IUPAC degenerate sequences
    Returns:
        a list of sequences that make up the degenerate sequences.
    """

    def expand_string(degenerate_string):
        expanded = []
        for char in degenerate_string:
            if char in codes:
                expanded.append(codes[char])
            else:
                expanded.append({char})
        return ["".join(chars) for chars in product(*expanded)]

    expanded_strings = set()
    for degenerate_string in degenerate_strings:
        expanded_strings.update(expand_string(degenerate_string))

    return list(expanded_strings)


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


def hammingdist(strings: list) -> np.array:
    """take a list of strings and return a condensed hamming distance matrix and the strings as labels
    Args:
        strings: a list of equal length sequence strings
    Returns:
        a condensed scipy distance matrix as a numpy array
    """
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
    root, nodes = scipy.cluster.hierarchy.to_tree(Z, rd=True)

    # Create a treelib tree
    tree = treelib.Tree()

    # Create a dictionary to store the nodes
    node_dict = {}

    # Add the nodes to the dictionary
    for i, node in enumerate(nodes):
        node_dict[node.id] = node

    # Add the root node to the tree with the tag "root"
    tree.create_node(tag="root", identifier=str(root.id))

    # Recursive function to add child nodes
    def add_children(parent_node, parent_id, labels):
        if parent_node.left is not None:
            # Add left child
            left_node = node_dict[parent_node.left.id]
            if left_node.left is None and left_node.right is None:
                # Leaf node, add label using its index in the nodes list
                tree.create_node(
                    tag=labels[left_node.id],
                    identifier=str(left_node.id),
                    parent=parent_id,
                )
            else:
                tree.create_node(tag="", identifier=str(left_node.id), parent=parent_id)
                add_children(left_node, str(left_node.id), labels)

        if parent_node.right is not None:
            # Add right child
            right_node = node_dict[parent_node.right.id]
            if right_node.left is None and right_node.right is None:
                # Leaf node, add label using its index in the nodes list
                tree.create_node(
                    tag=labels[right_node.id],
                    identifier=str(right_node.id),
                    parent=parent_id,
                )
            else:
                tree.create_node(
                    tag="", identifier=str(right_node.id), parent=parent_id
                )
                add_children(right_node, str(right_node.id), labels)

    # Start adding children from the root node
    add_children(root, str(root.id), labels)

    return tree


def find_degenerates(tree: treelib.Tree) -> list[str]:
    """
    Perform a breadth-first traversal of the treelib tree and process leaves.

    Parameters:
    - tree: A treelib.Tree object.

    Returns:
    - A list of strings returned by the create_degenerate function.
    """
    # Temporary list to store results from create_degenerate
    results = []
    # Expand the tree in breadth-first order to list tree node numbers as strings
    nodes = tree.expand_tree(mode=treelib.Tree.WIDTH)
    node_skip_list = []
    for node in nodes:
        if node not in node_skip_list:
            # Call create_degenerate with leaf tags
            seqs = [item.tag for item in tree.leaves(node)]
            result = create_degenerate(set(seqs))
            if result is not None:  # If it returns a string, save it
                results.append(result)
                node_skip_list.extend(list(tree.subtree(node).nodes))
    return results


def seqs_to_degenerates(seqs: list[str]) -> list[str]:
    """Takes a list of equal length sequences and finds a minimal set of degenerate codes  that represents them.
    Args:
        seqs: A list of sequences
    Returns:
        a list of degenerate sequences representing the input list
    """
    distance_matrix = hammingdist(seqs)
    tree = create_tree(seqs, distance_matrix)
    return find_degenerates(tree)
