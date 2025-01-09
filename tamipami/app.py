#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" app: a TamiPami module for a Streamlit app interface to TamiPami
"""

import os
import uuid
import gzip
import shutil
from pathlib import Path
import json

import streamlit as st
import altair as alt
import pandas as pd

from .config import config
from . import pam
from . import fastq
from . import degenerate
from . import tpio


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # This is your Project Root
# print("ROOT DIR")
# print(ROOT_DIR)
# create session-specific data dir


def create_session_dir():
    datadir = Path(os.path.join(ROOT_DIR, str(uuid.uuid4())))
    datadir.mkdir(parents=True, exist_ok=True)
    st.session_state["datadir"] = datadir


def delete_session_dir():
    if "datadir" in st.session_state:
        session_dir = st.session_state["datadir"]
    if os.path.exists(session_dir):
        shutil.rmtree(session_dir)
    del st.session_state["datadir"]


# Define input parameters and widgets

st.logo(os.path.join(ROOT_DIR, "assets/USDAARSIdentityRGB3.png"), size="large")
apptitle = "TamiPami"
st.set_page_config(page_title=apptitle, page_icon=":dna:")
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

with st.expander("Use Instructions"):
    st.markdown(
        """
                1. Load the forward and reverse FASTQ files for a single control and experimental library in the four input boxes on the sidebar.
                2. Select the Addgene library used for the experiment. Or...
                3. If no library is selected, enter the target sequence and the orientation (5prime is PAM-Target, 3prime is Target-PAM (like spCas9)).
                4. Select the maximum length to analyze. You can compare all smaller lengths once you have analyzed the data
                5. Hit 'Submit'. This will process your files.
                6. After you have hit 'Submit' you can explore your data interactively. The key interface is the Zscore slider bar. Moving that will set the 
                   cutoff value to separate kmers that cut from those that did not. This will update the histogram of sequence zscores, the table of reads  the degenerate sequences created and the sequence motif.
                7. Be sure to explore each length tab, to select the PAM/TAM site best supported by your data.
                8. Export the raw run data.  
    """
    )

args = {}
with st.sidebar:
    with st.form(key="newdata"):
        st.markdown("## Upload your FASTQ data")
        args["cont1"] = st.file_uploader(
            "Control library forward .fastq, .fq, .fastq.gz or .fq.gz file.",
            type=[".gz", ".fastq", ".fq"],
            key="cont1",
        )
        args["cont2"] = st.file_uploader(
            "Control library reverse .fastq, .fq, .fastq.gz or .fq.gz file.",
            type=[".gz", ".fastq", ".fq"],
            key="cont2",
        )
        args["exp1"] = st.file_uploader(
            "Experimental library forward .fastq, .fq, .fastq.gz or .fq.gz file.",
            type=[".gz", ".fastq", ".fq"],
            key="exp1",
        )
        args["exp2"] = st.file_uploader(
            "Experimental library reverse .fastq, .fq, .fastq.gz or .fq.gz file.",
            type=[".gz", ".fastq", ".fq"],
            key="exp2",
        )
        # args['log'] = st.sidebar.text_input('Log file', value=os.path.join(DATA_DIR,'tamipami.log'), key='log')
        st.markdown("## Set your search settings")
        args["length"] = st.select_slider(
            "The Maximum length of the PAM or TAM sequences",
            options=list(range(3, 9)),
            value=6,
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
    # with st.form(key='olddata'):
    #    st.markdown("## Or, visualize a previous run")
    #    args['olddata'] = st.file_uploader('a Tamipami output file in .json format', type=['.json'], key='olddata' )
    #    oldrun = st.form_submit_button('Submit' )


def write_input_file(datadir: str, stream, fname) -> None:
    try:
        filename = os.path.join(datadir, fname + ".fastq.gz")
        if stream.__getattribute__("type") == "application/x-gzip":
            with open(filename, "wb") as temp1:
                temp1.write(stream.getbuffer())
        else:
            with gzip.open(filename, "wb") as temp1:
                temp1.write(stream.getbuffer())
        temp1.close()
    except IOError as e:
        st.exception(e)
        raise (e)


def save_input_files(datadir: str, args: dict) -> None:
    argkeys = ["cont1", "cont2", "exp1", "exp2"]
    for akey, avalue in args.items():
        try:
            if akey in argkeys:
                assert (
                    avalue is not None
                ), "An input file is missing, please verify that the input files were uploaded correctly"
                write_input_file(datadir, avalue, akey)
        except AssertionError as e:
            st.exception(e)


# downloads
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
    try:
        create_session_dir()
        datadir = st.session_state["datadir"]
        save_input_files(datadir, args)
        with st.expander("Run Configuration"):
            st.dataframe(
                {"Parameter": args.keys(), "Value": args.values()}, hide_index=True
            )
            st.write("Data directory: {}".format(datadir))
        spacer, orientation = parse_lib(args)
        print("spacer: {}".format(spacer))
        print("orientation: {}".format(orientation))
        run_summ = {"cont": {}, "exp": {}}
        cont_raw, run_summ["cont"]["tot"], run_summ["cont"]["targets"] = fastq.process(
            fastq=os.path.join(datadir, "cont1.fastq.gz"),
            fastq2=os.path.join(datadir, "cont2.fastq.gz"),
            pamlen=args["length"],
            mergedfile=os.path.join(datadir, "cont_merged.fastq.gz"),
            spacer=spacer,
            orientation=orientation,
        )
        exp_raw, run_summ["exp"]["tot"], run_summ["exp"]["targets"] = fastq.process(
            fastq=os.path.join(datadir, "exp1.fastq.gz"),
            fastq2=os.path.join(datadir, "exp2.fastq.gz"),
            pamlen=args["length"],
            mergedfile=os.path.join(datadir, "exp_merged.fastq.gz"),
            spacer=spacer,
            orientation=orientation,
        )
        pamexpobj = pam.pamSeqExp(ctl=cont_raw, exp=exp_raw, position=orientation)
        return pamexpobj, run_summ
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
    print(results)
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


def main(args):
    if newrun:
        pamexpobj, run_summ = process(args)
        # if 'pamexpobj' not in st.session_state:
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
    # Check if pamdict is in session_state
    if "pamexpobj" in st.session_state:
        tablabels = [
            "Length " + str(x) for x in st.session_state.pamexpobj.multikmerdict.keys()
        ]
        st.markdown("### Select the PAM length to review")
        tabs = st.tabs(tablabels)
        for n, (key, df) in enumerate(st.session_state.pamexpobj.multikmerdict.items()):
            with tabs[n]:
                # Warning box
                data_checks(key, config["shot_cutoff"])
                read_count_check(
                    st.session_state["run_summ"], config["min_target_frac"]
                )
                st.markdown("### Select the cutoff point on the histogram")
                slider_key = f"slider_{key}"
                slider, slidebutton = st.columns(
                    spec=[0.8, 0.2], vertical_alignment="center"
                )

                with slider:
                    defaultval = 0.0
                    st.slider(
                        label="Select the Zscore cutoff:",
                        min_value=float(df["zscore"].min()),
                        max_value=float(df["zscore"].max()),
                        value=st.session_state.get(slider_key, defaultval),
                        key=slider_key,
                    )
                with slidebutton:
                    st.button(
                        "Auto Split",
                        on_click=update_slider,
                        kwargs={"n": key, "sliderkey": slider_key},
                        key=f"autosplit_{key}",
                    )

                # Filtered DataFrame based on slider value
                cutoff = st.session_state[slider_key]
                filtered_df = df[df["zscore"] >= cutoff]

                # create degernerate represenations
                althist = tpio.histogram_plot(
                    df, maxbins=config["histogram_bins"], cutoff=cutoff
                )
                dseqs = degenerate.seqs_to_degenerates(filtered_df["kmers"].tolist())
                # Plot histogram with vertical line at cutoff
                hist, degenerates = st.columns(spec=[0.8, 0.2])
                with hist:
                    st.altair_chart(althist)
                with degenerates:
                    st.dataframe({"PAM/TAM site": dseqs}, hide_index=True)

                st.subheader(f" Review filtered data for length {key}:")
                st.write(filtered_df)

                st.subheader(
                    f" Review sequence motif for length {key} and selected cutoff:"
                )
                logo = st.session_state.pamexpobj.make_logo(
                    length=key, cutoff=cutoff, type="zscore", above=False
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
            # data=download_json(st.session_state.pamexpobj.multikmerdict),
            # file_name="tamipami.json",
            # mime='application/json',
            key=f"download",
        )


if __name__ == "__main__":
    main(args)
    with st.expander("CCO Licence Information"):
        st.write(
            """
                As a work of the United States government, [USDA Agricultural Research Service](https://www.ars.usda.gov/), this project is in the public domain within the United States of America.
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
