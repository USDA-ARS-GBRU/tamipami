# -*- coding: utf-8 -*-
"""entropy :  a TamiPami module for working with entropy of sequence sets"""

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


def strings_to_char_array(strings: list[str]) -> np.ndarray:
    """
    Convert a list of equal-length strings into a 2D numpy character array.

    Args:
        strings (list of str): A list of strings, each of equal length.

    Returns:
        np.ndarray: A 2D numpy array where each column is a letter in the string
                    and each row is a string from the list, broken up by letter.
    """
    if not strings:
        raise ValueError("The input list must not be empty.")

    # Check that all strings are of equal length
    string_length = len(strings[0])
    if not all(len(s) == string_length for s in strings):
        raise ValueError("All strings must be of equal length.")

    # Convert the list of strings into a 2D numpy array
    char_array = np.array([list(s) for s in strings], dtype="str")

    return char_array


def calculate_character_fractions(char_array: np.ndarray) -> np.ndarray:
    """
    Calculate the fraction of each character ('A', 'G', 'C', 'T') in a specified range
    of rows and columns of a numpy character array.

    Args:
        char_array (np.ndarray): A 2D numpy character array.

    Returns:
        np.ndarray: A 2D numpy array with one row for each character ('A', 'G', 'C', 'T')
                    and one column for each selected column in the input array. Each element
                    represents the fraction of the column containing the respective character.
    """
    # Define the characters to count
    characters = np.array(["A", "G", "C", "T"])

    # Slice the array to get the specified range of rows and columns

    # Initialize an array to hold the fractions
    fractions = np.zeros((len(characters), char_array.shape[1]), dtype=float)

    # Calculate the fraction of each character in each column
    for i, char in enumerate(characters):
        # Count occurrences of the character in each column
        char_counts = np.sum(char_array == char, axis=0)
        # Calculate the fraction
        fractions[i] = char_counts / char_array.shape[0]

    return fractions


def entropy_of_frequency_array(freq_array):
    """
    Calculate the entropy for each column in a frequency array.

    Args:
        freq_array (np.ndarray): A 2D numpy array where each column represents
                                 the frequency distribution of a set of events.

    Returns:
        np.ndarray: A 1D numpy array containing the entropy of each column.
    """
    # Ensure the input is a numpy array
    freq_array = np.array(freq_array)

    # Normalize the frequency array to get probabilities
    prob_array = freq_array / freq_array.sum(axis=0, keepdims=True)

    # Calculate entropy for each column
    entropy_per_column = -np.nansum(
        np.where(prob_array > 0, prob_array * np.log2(prob_array), 0), axis=0
    )

    return entropy_per_column


def cumulative_entropy(strings: list[str]) -> np.ndarray:
    """Takes a list of ordered strings and returns the cumulative entropy  at each string position from right to left.

    Args:
        strings (list of str): A list of strings, each of equal length.

    Returns:
        np.ndarray: A 2D numpy array where each column is represents the position in the PAM/TAM
        and each row represents the cumulative entropy  to the set  of sequences from  that position until the end of the list
    """
    char_array = strings_to_char_array(strings)
    num_kmers = char_array.shape[0]
    cum_ent_list = []
    for i in range(num_kmers):
        sub_array = char_array[i:, :]
        sub_fracs = calculate_character_fractions(sub_array)
        cum_ent_list.append(entropy_of_frequency_array(sub_fracs))

    return np.array(cum_ent_list, dtype=np.float32)


def tot_ent(cumm_ent_array: np.ndarray, length: int, orientation: str) -> np.ndarray:
    """
    Calculate the total entropy for a given orientation and length from a cumulative entropy array.

    Args:
        cumm_ent_array (np.ndarray): A 2D array representing cumulative entropy values.
        length (int): The length of the sequence to consider for entropy calculation.
        orientation (str): The orientation of the sequence, either '3prime' or '5prime'.

    Returns:
        np.ndarray: An array containing the total entropy values for each sequence.

    Raises:
        ValueError: If the length is not within the valid range or if the orientation is invalid.
        TypeError: If the orientation is not a string.
    """
    cumm_ent_array = np.asarray(cumm_ent_array)
    if length < 0 or length > cumm_ent_array.shape[1]:
        raise ValueError("length must be less than the array width and greater than 0")
    if not isinstance(orientation, str):
        raise TypeError("orientation must be a string")
    if orientation == "3prime":
        sub_array = cumm_ent_array[:, :length]
    elif orientation == "5prime":
        sub_array = cumm_ent_array[:, length:]
    else:
        raise ValueError("Invalid orientation value. Expected '3prime' or '5prime'.")
    return np.sum(sub_array, axis=1)


def find_max_gradient_and_zscore(
    cumm_ent_array: np.ndarray, zscore: list, orientation: str
) -> pd.DataFrame:
    """
    Find the optimal cutpoint by combining cumulative entropy and z-score arrays,
    calculating gradients, and identifying maximum gradient values where z-score is positive.

    Args:
        cumm_ent_array (np.ndarray): A 2D numpy array of cumulative entropy values.
        zscore (np.ndarray): A 1D numpy array of z-score values.
        orientation (str): The orientation of the sequence, either '3prime' or '5prime'.

    Returns:
        pd.DataFrame: A DataFrame containing the maximum gradient and corresponding z-score
                      for each column in the cumulative entropy array.
    """
    # Ensure the input arrays are compatible
    if cumm_ent_array.shape[0] != len(zscore):
        raise ValueError("cumm_ent_array and zscore must have the same length")

    results = []
    # look at PAM candidates of 2 or more
    for i in range(1, cumm_ent_array.shape[1]):
        t_ent = tot_ent(
            cumm_ent_array=cumm_ent_array, length=i + 1, orientation=orientation
        )
        df = pd.DataFrame({"zscore": zscore, "t_ent": t_ent})
        df["gradient"] = np.gradient(
            df["zscore"],
            df["t_ent"],
        )
        filtered_df = df[df["zscore"] > 0]
        max_gradient_row = filtered_df.loc[filtered_df["gradient"].idxmax()]
        zscore_value = max_gradient_row["zscore"]
        gradient_value = max_gradient_row["gradient"]

        results.append(
            {
                "pam_length": i,
                "max_gradient": gradient_value,
                "zscore_at_max_gradient": zscore_value,
            }
        )
    return pd.DataFrame(results)


def analyze_clusters(df: pd.DataFrame) -> tuple[int, float]:
    """
    Perform k-means clustering on the DataFrame, identify the group with the largest
    max_gradient, and return the row with the largest pam_position in that group.

    Args:
        df (pd.DataFrame): A DataFrame with columns 'pam_position', 'max_gradient',
                           and 'zscore_at_max_gradient'.

    Returns:
        tuple: A tuple containing the pam_length zscore_at_max_gradient
               of the row with the largest pam_length in the group with the largest max_gradient.
    """
    # Perform k-means clustering with 2 clusters
    kmeans = KMeans(n_clusters=2, random_state=42)
    df["cluster"] = kmeans.fit_predict(df[["max_gradient", "zscore_at_max_gradient"]])
    print(df)
    # Identify the cluster with the largest max_gradient
    cluster_max_gradients = df.groupby("cluster")["max_gradient"].max()
    target_cluster = cluster_max_gradients.idxmax()

    # Filter the DataFrame to the target cluster
    target_group = df[df["cluster"] == target_cluster]

    # Find the row with the largest pam_position in the target cluster
    target_row = target_group.loc[target_group["pam_length"].idxmax()]
    # pam_length is zero-indexed so add one
    return (
        int(target_row["pam_length"] + 1),
        float(target_row["zscore_at_max_gradient"]),
    )


def main(
    kmers: list[str] | np.ndarray | pd.Series,
    zscore: list[str] | np.ndarray | pd.Series,
    orientation: str,
) -> tuple[int, float]:
    """
    Determine the optimal PAM/TAM length and breakpoint for sequence cutting based on cumulative entropy and z-scores.

    Args:
        kmers (list[str] | np.ndarray | pd.Series): A list-like object containing ordered kmers.
        zscore (list[str] | np.ndarray | pd.Series): A list-like object containing corresponding z-scores.
        orientation (str): The orientation of the sequence, either '3prime' or '5prime'.

    Returns:
        tuple: A tuple containing the best PAM/TAM length and the optimal breakpoint z-score.
    """
    if isinstance(kmers, (np.ndarray, pd.Series)):
        kmers = kmers.tolist()
    if isinstance(zscore, (np.ndarray, pd.Series)):
        zscore = zscore.tolist()
    cumm_ent_array = cumulative_entropy(kmers)
    mgz_df = find_max_gradient_and_zscore(
        cumm_ent_array=cumm_ent_array, zscore=zscore, orientation=orientation
    )
    return analyze_clusters(mgz_df)
