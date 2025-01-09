# -*- coding: utf-8 -*-
"""Utility functions to import and export and shared graphics 
"""
import pandas as pd
import logging
import sys
import tempfile
import os

import altair as alt


## HDF5 io ###


def export_hdf(
    multikmerdict: dict, to_buffer: bool = True, filename: str = "tamipami.h5"
) -> bytes:
    """A function to convert the pamseqboject multikmerdict to an HDF5 file"""
    if to_buffer:
        dcb = 0
    else:
        dcb = 1
    with pd.HDFStore(
        filename,
        mode="a",
        driver="H5FD_CORE",
        driver_core_backing_store=dcb,
        complevel=6,
    ) as out:
        for kmerlen, df in multikmerdict.items():
            key = "k_" + str(kmerlen)
            out.put(key, df, format="table", append=True)
        if to_buffer:
            return out._handle.get_file_image()
        else:
            out.close()


def _save_stdin_to_temp_file():
    # Create a named temporary file
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_name = temp_file.name
        # Read binary data from stdin
        binary_data = sys.stdin.buffer.read()
        # Write binary data to the temp file
        temp_file.write(binary_data)

    # Return the name of the temp file
    return temp_file_name


def import_hdf(from_buffer: bool = False, filename: str = None) -> dict:
    try:
        multikmerdict = {}
        if from_buffer:
            inputdat = _save_stdin_to_temp_file()
        else:
            inputdat = filename
        with pd.HDFStore(inputdat, "r") as store:
            kmerkeys = store.keys()
        for kmerlen in kmerkeys:
            multikmerdict[int(kmerlen[3:])] = pd.read_hdf(inputdat, key=kmerlen)
        return multikmerdict
    except Exception as e:
        logging.error("Could not import the hdf file or file stream")
        raise e
    finally:
        if from_buffer:
            os.remove(inputdat)


### Graphics ###


def histogram_plot(source, maxbins, cutoff=None, filename: str = None) -> alt.Chart:
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
    if cutoff is not None:
        # Vertical line to indicate the cutoff
        cutoff_line = (
            alt.Chart(pd.DataFrame({"cutoff": [cutoff]}))
            .mark_rule(color="red")
            .encode(x="cutoff:Q")
        )
        chart = chart + cutoff_line
    if filename:
        chart.save(filename)
    return chart
