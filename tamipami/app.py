#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""app: a TamiPami module for a Streamlit app interface to TamiPami

This app can process either:
- Two FASTQ files (forward reads for control and experimental samples, single-end mode)
- Four FASTQ files (forward and reverse reads for both control and experimental samples, paired-end mode)

You must provide either both reverse files (cont2 and exp2) for paired-end, or omit both for single-end.
"""

import os
import uuid
import gzip
import shutil
from pathlib import Path
import json

import streamlit as st
import pandas as pd

from tamipami.config import config
from tamipami import pam
from tamipami import fastq
from tamipami import degenerate
from tamipami import tpio
from tamipami._version import __version__


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # This is your Project Root


def create_session_dir():
    """
    Creates a unique session directory for storing session-specific data.

    Generates a secure UUID to name the directory, creates the directory
    within the ROOT_DIR, and stores the path in Streamlit's session state.
    If directory creation fails, an error message is displayed.

    Raises:
        Exception: If the session directory cannot be created.
    """
    try:
        session_id = uuid.uuid4()
        datadir = Path(os.path.join(ROOT_DIR, str(session_id)))
        datadir.mkdir(parents=True, exist_ok=True)
        st.session_state["datadir"] = datadir
    except Exception as e:
        st.error(f"Failed to create session directory: {e}")


def delete_session_dir():
    """
    Delete the session directory stored in Streamlit's session state.

    This function checks if a 'datadir' key exists in the session state.
    If it does, it retrieves the directory path, deletes the directory
    if it exists, and then removes the 'datadir' key from the session state.
    """
    if "datadir" in st.session_state:
        session_dir = st.session_state["datadir"]
        if os.path.exists(session_dir):
            shutil.rmtree(session_dir)
        del st.session_state["datadir"]


# Define input parameters and widgets
apptitle = "TamiPami"
st.set_page_config(page_title=apptitle, page_icon=":material/genetics:")
st.logo(os.path.join(ROOT_DIR, "assets/ARS_UF_combined.png"), size="large")
st.image(os.path.join(ROOT_DIR, "assets/tami_postcard.jpeg"))
st.subheader("Identify the PAMs of new Cas enzymes or TAMs of TnpB endonucleases")

st.markdown(
    """
            ## Overview 
            When a new Cas or TnpB guided endonuclease is discovered or engineered, one of the first tasks is identifying its PAM or TAM recognition site. 
            This can be done by treating a pool of plasmid DNA containing a target site adjacent to a random region. The random regions containing 
            the PAM/TAM site are recognized and cut in the presence of the endonuclease and a guide RNA. These become depleted in the sequencing library. By comparing an uncut control library  to a cut 
            experimental library it is possible to identify the PAM/TAM site.  The method was introduced by [Walton et al. 2021]( https://doi.org/10.1038/s41596-020-00465-2).
            That work deposited the plasmid pools with [Addgene]( https://www.addgene.org/pooled-library/kleinstiver-ht-pamda/), making the lab protocol accessible.   
            
            This web application builds on the work by creating TamiPami, a Web application that simplifies the analysis of the sequencing data and adds rich interactive 
            visualizations for selecting the PAM/TAM site.

        """
)

with st.expander(label="Use Instructions", icon=":material/help_center:"):
    st.markdown(
        """
                1. **Upload your FASTQ files:**  
                   - For single-end: Upload only the forward reads for control and experimental samples (cont1 and exp1).  
                   - For paired-end: Upload both forward and reverse reads for control and experimental samples (cont1, cont2, exp1, exp2).  
                   - You must provide either both reverse files (cont2 and exp2) for paired-end, or omit both for single-end.
                2. Select the Addgene library used for the experiment. Or...
                3. If no library is selected, enter the target sequence and the orientation (5prime is PAM-Target, 3prime is Target-PAM (like spCas9)).
                4. Select the maximum length to analyze. You can compare all smaller lengths once you have analyzed the data
                5. Hit 'Submit'. This will process your files.
                6. After you have hit 'Submit' you can explore your data interactively. The key interface is the Zscore slider bar. Moving that will set the 
                cutoff value to separate kmers that cut from those that did not. This will update the histogram of sequence zscores, the table of reads  the degenerate sequences created and the sequence motif.
                7. Be sure to explore each length tab, to select the PAM/TAM site best supported by your data.
                8. Export the raw run data.

                If you need help or encounter errors please request support on the [Tamipami Github issues page](https://github.com/USDA-ARS-GBRU/tamipami/issues)   
    """
    )

args = {}
with st.sidebar:
    with st.expander(label="Need test data for TamiPami?", icon=":material/genetics:"):
        st.markdown(
            """
                    Download the zip file `Example.zip` with the button below. It contains 4 files: 
                    - `example_cont1.fastq.gz`
                    - `example_cont2.fastq.gz`
                    - `example_exp1.fastq.gz`
                    - `example_exp2.fastq.gz`

                    Each file contains 250,000 reads. the "cont" files contain the forward and reverse reads for the
                    undigested "RTW554" plasmid pool library. Those are uploaded in the control forward and control reverse
                    boxes below. The "exp" files contain the forward and reverse reads for the **RTW554** library digested with
                    SpCas9 nuclease. The example data is a sample of the first 250,000 read pairs from NCBI SRA **SRR34761011**
                    (cont) and **SRR34761004** (exp).
                    """
        )
        file_path = os.path.join(ROOT_DIR, "assets/Example.zip")
        with open(file_path, "rb") as f:
            file_bytes = f.read()
        st.download_button(
            label="Download example SpCas9 data",
            data=file_bytes,
            file_name="Example.zip",
            mime="application/zip",  # or "application/gzip" etc.,
            type="primary",
        )
    with st.form(key="newdata"):
        st.markdown("## Upload your FASTQ data")
        st.markdown(
            """
            - For **single-end**: Upload only the forward reads for control and experimental samples (cont1 and exp1).
            - For **paired-end**: Upload both forward and reverse reads for control and experimental samples (cont1, cont2, exp1, exp2).
            - You must provide either both reverse files (cont2 and exp2) for paired-end, or omit both for single-end.
            """
        )
        args["cont1"] = st.file_uploader(
            "Control library forward .fastq, .fq, .fastq.gz or .fq.gz file.",
            type=[".gz", ".fastq", ".fq"],
            key="cont1",
        )
        args["cont2"] = st.file_uploader(
            "Control library reverse .fastq, .fq, .fastq.gz or .fq.gz file (optional, required for paired-end).",
            type=[".gz", ".fastq", ".fq"],
            key="cont2",
        )
        args["exp1"] = st.file_uploader(
            "Experimental library forward .fastq, .fq, .fastq.gz or .fq.gz file.",
            type=[".gz", ".fastq", ".fq"],
            key="exp1",
        )
        args["exp2"] = st.file_uploader(
            "Experimental library reverse .fastq, .fq, .fastq.gz or .fq.gz file (optional, required for paired-end).",
            type=[".gz", ".fastq", ".fq"],
            key="exp2",
        )
        st.markdown("## Set your search settings")
        args["length"] = st.select_slider(
            "The maximum length of the PAM or TAM sequences",
            options=list(range(3, 7)),
            value=5,
            key="length",
        )
        args["library"] = st.selectbox(
            "The Addgene library pool. For custom pools use the --spacer and --orientation flags",
            ["RTW554", "RTW555", "RTW572", "RTW574"],
            index=None,
        )
        st.markdown("_OR_")
        args["spacer"] = st.text_input(
            "The spacer sequence for the guide RNA. Not needed if ---library is used",
            key="spacer",
        )
        args["orientation"] = st.selectbox(
            "The side of the spacer the PAM/TAM is on",
            ["3prime", "5prime"],
            index=None,
            key="orientation",
        )
        newrun = st.form_submit_button("Submit")


def write_input_file(datadir: str, stream, fname) -> None:
    """
    Writes a FASTQ file to the specified directory from a given stream.

    This function checks if the specified directory exists and is writable.
    If the directory does not exist, it creates it. The function then writes
    the content of the provided stream to a gzipped FASTQ file in the directory.
    The file is named using the provided filename with a '.fastq.gz' extension.

    Parameters:
        datadir (str): The directory path where the file will be saved.
        stream: The input stream containing the file data.
        fname: The base name for the output file.

    Raises:
        PermissionError: If the directory is not writable.
        IOError: If an I/O error occurs during file operations.
    """
    try:
        if not os.path.exists(datadir):
            os.makedirs(datadir)
        if not os.access(datadir, os.W_OK):
            raise PermissionError(f"Directory {datadir} is not writable.")

        filename = os.path.join(datadir, f"{fname}.fastq.gz")
        if stream.__getattribute__("type") == "application/x-gzip":
            with open(filename, "wb") as temp1:
                temp1.write(stream.getbuffer())
        else:
            with gzip.open(filename, "wb") as temp1:
                temp1.write(stream.getbuffer())
    except (IOError, PermissionError) as e:
        st.exception(e)
        raise (e)


def save_input_files(datadir: str, args: dict) -> None:
    """
    Save uploaded FASTQ files to the session directory.

    This function validates that either both reverse files (cont2 and exp2) are provided (paired-end),
    or both are omitted (single-end). It writes only the files that are present.

    Parameters:
        datadir (str): The directory path where the files will be saved.
        args (dict): Dictionary containing file upload objects.

    Raises:
        ValueError: If only one of cont2 or exp2 is provided.
        AssertionError: If required forward files are missing.
    """
    # Validation for paired-end vs single-end
    cont2_present = args.get("cont2") is not None
    exp2_present = args.get("exp2") is not None
    if (cont2_present and not exp2_present) or (exp2_present and not cont2_present):
        raise ValueError(
            "Invalid input: provide both reverse files (cont2 and exp2) for paired-end processing, or omit both for single-end."
        )
    # Always require forward files
    if args.get("cont1") is None or args.get("exp1") is None:
        raise AssertionError(
            "Both control and experimental forward files (cont1 and exp1) are required."
        )

    # Write files that are present
    for key in ["cont1", "cont2", "exp1", "exp2"]:
        stream = args.get(key)
        if stream is not None:
            write_input_file(datadir, stream, key)


def download_json(dfdict):
    tempdict = {}
    for length, df in dfdict.items():
        tempdict[length] = df.to_dict()
    json_data = json.dumps(tempdict)
    json_bytes = json_data.encode("utf-8")
    return json_bytes


def parse_lib(args: dict) -> tuple[str]:
    if args["library"]:
        spacer = config["spacer_dict"][args["library"]]["spacer"]
        orientation = config["spacer_dict"][args["library"]]["orientation"]
    else:
        spacer = args["spacer"]
        orientation = args["orientation"]
    return spacer, orientation


@st.cache_data
def process(args):
    """
    Process control and experimental FASTQ files to generate a pamSeqExp object and run summary.

    You must provide either:
    - Only cont1 and exp1 for single-end mode (forward reads only)
    - Or all of cont1, cont2, exp1, exp2 for paired-end mode (forward and reverse reads)

    Raises:
        ValueError: If only one of cont2 or exp2 is provided.
        AssertionError: If required forward files are missing.
        Exception: For unexpected errors during processing.
    """
    try:
        create_session_dir()
        datadir = st.session_state["datadir"]
        save_input_files(datadir, args)

        def safe_value(val):
            if hasattr(val, "name"):
                return val.name
            return str(val)

        with st.expander("Run Configuration"):
            param_df = pd.DataFrame(
                {
                    "Parameter": list(args.keys()),
                    "Value": [safe_value(v) for v in args.values()],
                }
            )
            st.dataframe(param_df, hide_index=True)
            st.write("Data directory: {}".format(datadir))
        spacer, orientation = parse_lib(args)
        run_summ = {"cont": {}, "exp": {}}

        # Determine paired-end or single-end
        cont2_path = (
            os.path.join(datadir, "cont2.fastq.gz") if args.get("cont2") else None
        )
        exp2_path = os.path.join(datadir, "exp2.fastq.gz") if args.get("exp2") else None

        cont_raw, run_summ["cont"]["tot"], run_summ["cont"]["targets"] = fastq.process(
            fastq=os.path.join(datadir, "cont1.fastq.gz"),
            fastq2=cont2_path,
            pamlen=args["length"],
            spacer=spacer,
            orientation=orientation,
        )
        exp_raw, run_summ["exp"]["tot"], run_summ["exp"]["targets"] = fastq.process(
            fastq=os.path.join(datadir, "exp1.fastq.gz"),
            fastq2=exp2_path,
            pamlen=args["length"],
            spacer=spacer,
            orientation=orientation,
        )
        pamexpobj = pam.pamSeqExp(ctl=cont_raw, exp=exp_raw, position=orientation)
        return pamexpobj, run_summ
    except Exception as e:
        st.error(
            f"There was an error processing the FASTQ files, Please verify your input files."
        )
        raise e
    finally:
        delete_session_dir()


def update_slider(n, sliderkey):
    breakpoint = st.session_state.pamexpobj.find_breakpoint(length=n, type="zscore")
    st.session_state[sliderkey] = breakpoint


def data_checks(length, sdcutoff):
    results = st.session_state.pamexpobj.check_N(length)
    mess = "Shot (Poisson) noise of the control library is {:.1%} and the experimental library is {:.1%}.".format(
        results[0], results[1]
    )
    recc = "Consider analyzing your data at a shorter length or sequencing more deeply"
    if any([x > sdcutoff for x in results]):
        return st.warning(mess + recc)
    else:
        return st.info(mess)


def read_count_check(run_summ, minfrac):
    expfrac = run_summ["exp"]["targets"] / run_summ["exp"]["targets"]
    contfrac = run_summ["cont"]["targets"] / run_summ["cont"]["targets"]
    if (expfrac or contfrac) < minfrac:
        return st.warning(
            "Only {:.1%} of merged experimental reads and {:.1%} \
                        of merged control reads contained a target. Please verify your Library, Spacer and Orientation settings".format(
                expfrac, contfrac
            )
        )


# Wrapper functions for caching


@st.cache_data
def compute_degenerate(filtered_kmers_tuple):
    return degenerate.seqs_to_degenerates(list(filtered_kmers_tuple))


def main(args):
    if newrun:
        try:
            pamexpobj, run_summ = process(args)
        except (ValueError, AssertionError) as e:
            st.error(str(e))
            return
        st.session_state["run_summ"] = run_summ
        st.session_state["pamexpobj"] = pamexpobj
        st.markdown("## Library information")
        spacer, orientation = parse_lib(args)
        if "library" in args:
            st.write("Plasmid Library: {}".format(args["library"]))
        st.write("Library spacer sequence: {}".format(spacer))
        st.write("PAM/TAM orientation: {}".format(orientation))
        if orientation == "5prime":
            st.image(os.path.join(ROOT_DIR, "assets/5prime.jpg"))
            st.write(
                "image from [Walton et al. 2021]( https://doi.org/10.1038/s41596-020-00465-2)"
            )
        elif orientation == "3prime":
            st.image(os.path.join(ROOT_DIR, "assets/3prime.jpg"))
            st.write(
                "image from [Walton et al. 2021]( https://doi.org/10.1038/s41596-020-00465-2)"
            )

    if "pamexpobj" in st.session_state:
        length_keys = list(st.session_state.pamexpobj.multikmerdict.keys())
        tablabels = [f"Length {x}" for x in length_keys]
        num_tabs = len(tablabels)

        # Maintain persistent tab selection in session state
        if "active_length_tab" not in st.session_state:
            st.session_state["active_length_tab"] = 0

        # The last active tab index is restored by using `active_length_tab`
        selected_idx = st.session_state.get("active_length_tab", 0)
        tabs = st.tabs(tablabels)

        for i, (tab, key) in enumerate(zip(tabs, length_keys)):
            with tab:
                # Assign current tab as active if visible
                st.session_state["active_length_tab"] = i
                df = st.session_state.pamexpobj.multikmerdict[key]
                slider_key = f"slider_{key}"

                # Initialize slider at autosplit value on first visit
                if slider_key not in st.session_state:
                    update_slider(key, slider_key)

                slider, slidebutton = st.columns(
                    spec=[0.8, 0.2], vertical_alignment="center"
                )
                with slider:
                    # Slider always uses key (session-state-controlled)
                    st.slider(
                        label="Select the Zscore cutoff:",
                        min_value=float(df["zscore"].min()),
                        max_value=float(df["zscore"].max()),
                        key=slider_key,
                    )
                with slidebutton:
                    st.button(
                        "Auto Split",
                        on_click=update_slider,
                        kwargs={"n": key, "sliderkey": slider_key},
                        key=f"autosplit_{key}",
                    )

                # --- Key FIX: Always get cutoff from session state ---
                cutoff = (
                    st.session_state[slider_key]
                    if slider_key in st.session_state
                    else None
                )

                althist = tpio.histogram_plot(
                    df,
                    maxbins=config["histogram_bins"],
                    cutoff=cutoff,  # Always pass the cutoff, even on tab change
                )

                hist, degenerates = st.columns(
                    spec=[0.8, 0.2], vertical_alignment="center"
                )
                with hist:
                    st.altair_chart(althist)
                filtered_df = df[df["zscore"] >= cutoff] if cutoff is not None else df
                dseqs = compute_degenerate(filtered_df["kmers"].tolist())
                with degenerates:
                    st.dataframe({"PAM/TAM site": dseqs}, hide_index=True)

                st.subheader(f" Review filtered data for length {key}:")
                styled_df = filtered_df.style.format(
                    {"pvalue": "{:.3e}", "p_adjust_BH": "{:.3e}"}
                )
                st.write(styled_df)

                st.subheader(
                    f" Review sequence motif for length {key} and selected cutoff:"
                )
                logo = st.session_state.pamexpobj.make_logo(
                    length=key, cutoff=cutoff, score_type="zscore", above=False
                )

                st.pyplot(logo)

        st.divider()
        st.markdown(
            """
            ### Data Export
            Data from the run can be exported as a json for further analysis.  This file contains the 
            unfiltered data tables for each length. It does not contain the selected cutoffs, plots 
            or degenerate sequences. Those can be downloaded by selecting the elements on the screen.
            """
        )
        st.download_button(
            label="📥 Download the full run data as a HDF5 file",
            data=tpio.export_hdf(
                st.session_state.pamexpobj.multikmerdict,
                to_buffer=True,
                filename="tamipami.h5",
            ),
            file_name="tamipami.h5",
            mime="application/x-hdf5",
            key=f"download",
        )


if __name__ == "__main__":
    main(args)
    with st.expander(label="CCO Licence Information", icon=":material/copyright:"):
        st.write(
            """
                As a work of the United States government, [USDA Agricultural Research Service](https://www.ars.usda.gov/), this project is in the public domain within the United States of America under 17 U.S.C. 105.
                Additionally, we waive copyright and related rights in the work worldwide through the CC0 1.0 Universal public domain dedication.
                ## CC0 1.0 Universal Summary
                This is a human-readable summary of the [Legal Code (read the full text)](https://creativecommons.org/publicdomain/zero/1.0/legalcode).
                ### No Copyright
                The person who associated a work with this deed has dedicated the work to the public domain by waiving all of their rights to the work worldwide under copyright law, including all related and neighboring rights, to the extent allowed by law.
                You can copy, modify, distribute, and perform the work, even for commercial purposes, all without asking permission.
                ### Other Information
                In no way are the patent or trademark rights of any person affected by CC0, nor are the rights that other persons may have in the work or in how the work is used, such as publicity or privacy rights.
                Unless expressly stated otherwise, the person who associated a work with this deed makes no warranties about the work, and disclaims liability for all uses of the work, to the fullest extent permitted by applicable law. When using or citing the work, you should not imply endorsement by the author or the affirmer.
                """
        )
    with st.expander(label="Citation", icon=":material/ink_pen:"):
        st.write(
            "We are preparing the TamiPami mansucript: Orosco, Carlos. Jain, Piyush. Rivers, Adam R. TamiPami: a TAM/PAM identification interface for CRISPR and Omega systems. In Prep."
        )
        st.write(
            "[Tamipami Github repository](https://github.com/USDA-ARS-GBRU/tamipami)"
        )
st.write("Tamipami version {}".format(__version__))
