# -*- coding: utf-8 -*-
"""Utility functions to import and export shared data and graphics 
"""
import sys
import tempfile
import logging
import tempfile
import os

import pandas as pd
import altair as alt


## HDF5 io ###


def export_hdf(
    multikmerdict: dict, to_buffer: bool = True, filename: str = "tamipami.h5"
) -> bytes:
    """
    Convert a dictionary of DataFrames to an HDF5 file.

    This function takes a dictionary where keys represent k-mer lengths and values
    are pandas DataFrames, and exports them to an HDF5 file. Each DataFrame is stored
    under a key prefixed with 'k_' followed by the k-mer length. The function can
    either write directly to a file or return the file as a byte buffer.

    Parameters:
        multikmerdict (dict): A dictionary with k-mer lengths as keys and DataFrames as values.
        to_buffer (bool): If True, returns the HDF5 file as a byte buffer. Defaults to True.
        filename (str): The name of the HDF5 file to write to. Defaults to 'tamipami.h5'.

    Returns:
        bytes: The HDF5 file as a byte buffer if to_buffer is True, otherwise None.
    """
    if to_buffer:
        driver_core_backing_store = 0
    else:
        driver_core_backing_store = 1
    try:
        with pd.HDFStore(
            filename,
            mode="a",
            driver="H5FD_CORE",
            driver_core_backing_store=driver_core_backing_store,
            complevel=6,
        ) as out:
            for kmerlen, df in multikmerdict.items():
                key = "k_" + str(kmerlen)
                out.put(key, df, format="table", append=True)
            if to_buffer:
                return out._handle.get_file_image()
            else:
                out.close()
    except (OSError, KeyError, ValueError, TypeError, pd.errors.PossibleDataLossError) as e:
        logging.error(f"Could not export the HDF5 file: {e}")
        raise e


def import_hdf(from_buffer: bool = False, filename: str = None) -> dict:
    """
    Imports data from an HDF file or a temporary file created from stdin.

    Parameters:
        from_buffer (bool): If True, reads data from stdin and saves it to a temporary file.
        filename (str): The path to the HDF file to import. Required if 'from_buffer' is False.

    Returns:
        dict: A dictionary where keys are k-mer lengths and values are the corresponding dataframes.

    Raises:
        ValueError: If 'from_buffer' is False and no filename is provided.
        OSError, KeyError, pd.errors.HDFStoreError: If there is an error reading the HDF file or stream.
    """
    try:
        multikmerdict = {}
        inputdat: str
        if from_buffer:
            inputdat = store_stdin_binary_to_tempfile()
        else:
            if filename is None:
                raise ValueError(
                    "A valid filename must be provided when 'from_buffer' is False."
                )
            inputdat = filename
        with pd.HDFStore(inputdat, "r") as store:
            kmerkeys = store.keys()
        for kmerlen in kmerkeys:
            multikmerdict[int(kmerlen[3:])] = pd.read_hdf(inputdat, key=kmerlen)
        return multikmerdict
    except (OSError, KeyError, ValueError, pd.errors.PossibleDataLossError) as e:
        logging.error(f"Could not import the hdf file or file stream: {e}")
        raise e
    finally:
        if from_buffer and os.path.exists(inputdat):
            os.remove(inputdat)


def store_stdin_binary_to_tempfile():
    """
    Creates a temporary file to store binary data read from standard input (stdin).

    This function checks if stdin is not empty and reads binary data from it.
    The data is then written to a named temporary file, which is not deleted
    automatically. The function logs the creation of the temporary file and
    any errors encountered during the process. It returns the name of the
    temporary file if successful, or None if an error occurs.
    """
    # Check if stdin is not empty
    if not sys.stdin.isatty():
        try:
            # Create a named temporary file
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                temp_file_name = temp_file.name
                logging.info(f"Temporary file created: {temp_file_name}")
                # Read binary data from stdin
                binary_data = sys.stdin.buffer.read()
                # Write binary data to the temp file
                temp_file.write(binary_data)
        except (OSError, IOError) as e:
            logging.error(f"An error occurred: {e}")
            return None

        logging.info(f"Temporary file will be returned: {temp_file_name}")
        # Return the name of the temp file
        return temp_file_name


### Graphics ###


def histogram_plot(
    source: pd.DataFrame, maxbins: int, cutoff: float = None, filename: str = None
) -> alt.Chart:
    """
    Generates a histogram plot from a pandas DataFrame using Altair.

    Parameters:
        source (pd.DataFrame): The data source for the histogram.
        maxbins (int): The maximum number of bins for the histogram.
        cutoff (float, optional): A cutoff value to draw a vertical line on the plot.
        filename (str, optional): The file path to save the plot. If not provided, the plot is not saved.

    Returns:
        alt.Chart: The generated Altair chart object.

    Raises:
        ValueError: If `source` is not a DataFrame, `maxbins` is not a positive integer, or `cutoff` is not a number.
    """
    # Input validation
    if not isinstance(source, pd.DataFrame):
        raise ValueError("`source` must be a pandas DataFrame.")
    if not isinstance(maxbins, int) or maxbins <= 0:
        raise ValueError("`maxbins` must be a positive integer.")
    if cutoff is not None and not isinstance(cutoff, (int, float)):
        raise ValueError("`cutoff` must be a number if provided.")
    # Base chart with bars
    chart = (
        alt.Chart(source)
        .mark_bar(binSpacing=0)
        .encode(
            alt.X("zscore", bin=alt.BinParams(maxbins=maxbins)),
            y="count()",
        )
        .properties(width=700, height=500)
        .interactive()
    )

    def create_cutoff_line(cutoff):
        return (
            alt.Chart(pd.DataFrame({"cutoff": [cutoff]}))
            .mark_rule(color="red")
            .encode(x="cutoff:Q")
        )

    if cutoff is not None:
        cutoff_line = create_cutoff_line(cutoff)
        chart = chart + cutoff_line

    if filename:
        try:
            chart.save(filename)
            logging.info(f"Chart successfully saved to {filename}")
        except IOError as e:
            logging.error(f"Failed to save chart to {filename}: {e}")

    return chart
