

import os
import uuid
import gzip
import shutil
from pathlib import Path

import streamlit as st
import altair as alt
import pandas as pd

from config import config
import pam
import fastq 
import degenerate

st.logo("USDAARSIdentityRGB3.png", size= "large")
apptitle = 'TamiPami'
st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")
st.sidebar.markdown("## Select Parameters to identiy the PAM/TAM")

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))  # This is your Project Root

# create session-specific data dir

def create_session_dir():
    datadir = Path(os.path.join(ROOT_DIR, str(uuid.uuid4())))
    datadir.mkdir(parents=True, exist_ok=True)
    st.session_state['datadir'] = datadir


def delete_session_dir():
    if 'datadir' in st.session_state:
        session_dir = st.session_State['datadir']
    if os.path.exists(session_dir):
        shutil.rmtree(session_dir)
        st.info(f"Deleted session directory: {session_dir}")

#Define input parameters and widgets

st.title('TamiPami')
st.subheader("Identify PAM or TAM sites in new Cas enzymes or  TnpB Transposoons")

st.markdown('''
            ## Overview 
            When a new Cas or IS200 transposon is discovered or engineered one of the first tasks is identifying its PAM or TAM recognition site. 
             This can be done by treating a pool of plasmids containing a target site adjacent to a random region. The random regains containing 
            the PAM/TAM site are recognized and cut. These become depleted in the sequencing library.  By comparing an uncut library to a cut 
            library it is possible to identify the PAM site.  The method was introduced by [Walton et al. 2021]( https://doi.org/10.1038/s41596-020-00465-2).
             That work deposited the plasmid pools with [Addgene]( https://www.addgene.org/pooled-library/kleinstiver-ht-pamda/), making the lab protocol accessible.   
            
            This web application builds on the work by simplifying the analysis of the sequencing data and adding rich interactive visualizations for seleting the PAM/TAM site.

         ''')

with st.expander("Instructions"):
    st.markdown('''
                1. Load the forward and reverse FASTQ files for a single control and experimental library in the four input boxes on the sidebar
                2. Select the Addgene library used for the experiment OR
                3. If no library is selected, enter the target sequence and the orientation (5prime is PAM-Target, 3prime is Target-PAM (like spCas9))
                4. Select the maximum length to analize. 5or 6 is a good length. you can compare all smaller lengths once you have analized the data
                5. Hit Run. This will process your files.
                6. After you have hit run you can explore your data interactivly.  The key interface is the Zscore slider bar. Moving that will set the cutoff.
                   This will update the histogram of seuence abundance, the  table of reads  the degenerate sequences created and the sequence motif.
                7. Be sure to explore each length tab, to select the PAM/TAM site best supported by your data.
                8. Export the raw run data.  
    ''')

args = {}
args['cont1'] = st.sidebar.file_uploader('Control library forward .fastq, .fq, .fastq.gz or .fq.gz file.', type=['.gz', '.fastq', '.fq'], key='cont1' )
args['cont2'] = st.sidebar.file_uploader('Control library reverse .fastq, .fq, .fastq.gz or .fq.gz file.',  type=['.gz', '.fastq', '.fq'],  key='cont2')
args['exp1'] = st.sidebar.file_uploader('Experimental library forward .fastq, .fq, .fastq.gz or .fq.gz file.',type=['.gz', '.fastq', '.fq'], key='exp1')
args['exp2'] = st.sidebar.file_uploader('Experimental library reverse .fastq, .fq, .fastq.gz or .fq.gz file.',type=['.gz', '.fastq', '.fq'], key='exp2',)
args['library'] = st.sidebar.selectbox('The Addgene library pool. For custom pools use the --spacer and --orientation flags', ["RTW554", "RTW555", "RTW572", "RTW574"], index=None)
args['spacer'] = st.sidebar.text_input('The spacer sequence for the guide RNA. Not needed if ---library is used', key='spacer')
args['orientation'] = st.sidebar.selectbox('The side of the spacer the PAM/TAM is on', ["3prime", "5prime"],index=None, key='orientation')
# args['log'] = st.sidebar.text_input('Log file', value=os.path.join(DATA_DIR,'tamipami.log'), key='log')
args['length'] = st.sidebar.select_slider('The Maxumum length of the PAM or TAM sequences', options=list(range(3, 9)), value=6, key='length')


def write_input_file(datadir: str, stream, fname) -> None:
        try:
            filename = os.path.join(datadir, fname + ".fastq.gz")
            if stream.__getattribute__("type") == 'application/x-gzip':
                with open(filename,'wb') as temp1:
                    temp1.write(stream.getbuffer()) 
            else:
                with gzip.open(filename,'wb') as temp1:
                    temp1.write(stream.getbuffer())  
            temp1.close()
        except IOError as e :
            st.exception(e)
            raise(e) 

def save_input_files(datadir: str, args: dict) -> None:
    argkeys = ["cont1", "cont2", "exp1", "exp2"]
    for akey, avalue in args.items():
        try:
            if akey in argkeys:
                assert avalue is not None, "An input file is missing, please verify that the input files were uploaded correctly"
                write_input_file(datadir, avalue, akey)
        except AssertionError as e:
            st.exception(e)


#@st.cache_data
def histogram_plot(source, maxbins, cutoff=None):
    # Base chart with bars
    chart = alt.Chart(source).mark_bar(binSpacing=0).encode(
        alt.X("zscore", bin=alt.BinParams(maxbins=maxbins)),
        y='count()',
    ).properties(
        width=800,
        height=600
    ).interactive()

    if cutoff is not None:
        # Vertical line to indicate the cutoff
        cutoff_line = alt.Chart(pd.DataFrame({'cutoff': [cutoff]})).mark_rule(color='red').encode(
            x='cutoff:Q'
        )
        chart = chart + cutoff_line

    return chart


@st.cache_data
def process(datadir, args):
    with st.expander("Run Configuration"):
        st.dataframe({"Parameter" :  args.keys(), "Value":args.values()},hide_index=True)
        st.write("Data directory: {}".format(DATA_DIR))
    if args['library']:
        spacer = config['spacer_dict'][args['library']]['spacer']
        orientation = config['spacer_dict'][args['library']]['orientation']
    else:
        spacer = args['spacer']
        orientation = ['orientation']
    cont_raw = fastq.process(fastq=os.path.join(datadir, "cont1.fastq.gz"),
                             fastq2=os.path.join(datadir, "cont2.fastq.gz"),
                             pamlen=args['length'],
                             mergedfile=os.path.join(datadir,"cont_merged.fastq.gz"),
                             spacer=spacer,
                             orientation=orientation,
                             )
    exp_raw = fastq.process(fastq=os.path.join(datadir, "exp1.fastq.gz"),
                            fastq2=os.path.join(datadir, "exp2.fastq.gz"),
                            pamlen=args['length'],
                            mergedfile=os.path.join(datadir,"exp_merged.fastq.gz"),
                            spacer=spacer,
                            orientation=orientation,
                            )
    pamexpobj = pam.pamSeqExp(ctl=cont_raw, exp=exp_raw, position=orientation)
    return pamexpobj

def update_slider(n, sliderkey):
    breakpoint = st.session_state.pamexpobj.find_breakpoint(length=n, type='zscore')
    st.session_state[sliderkey] = breakpoint

def main(args):
    if st.sidebar.button("Run"):
        save_input_files(DATA_DIR, args)
        pamexpobj= process(DATA_DIR, args)
        #if 'pamexpobj' not in st.session_state:
        st.session_state['pamexpobj'] = pamexpobj

    # Check if pamdict is in session_state
    if 'pamexpobj' in st.session_state:
        tablabels = ["Length " + str(x) for x in st.session_state.pamexpobj.multikmerdict.keys()]
        tabs = st.tabs(tablabels)

        for n, (key, df) in enumerate(st.session_state.pamexpobj.multikmerdict.items()):

            # Warning box
            ctl_raw_sd, exp_raw_sd = st.session_state.pamexpobj.check_N(key)
            if ctl_raw_sd > 0.25:
                tabs[n].warning("Poisson noise of the control library is expected to be {:.1%} of the total variance. Consider using smaller length or sequencing more deeply".format(ctl_raw_sd))
            else:
                tabs[n].info("Poisson noise of the control library is expected to be {:.1%} of the total variance".format(ctl_raw_sd))
            if exp_raw_sd > 0.25:
                tabs[n].warning("Poisson noise of the experimental library is expected to be {:.1%} of the total variance. Consider using smaller length or sequencing more deeply".format(ctl_raw_sd))
            else:
                tabs[n].info("Poisson noise of the experimental library is expected to be {:.1%} of the total variance".format(ctl_raw_sd))
            # Update the slider value in session state
            slider_key = f"slider_{key}"
            #defaultval = st.session_state.pamexpobj.find_breakpoint(length=key, type='zscore')
            defaultval = 0.0
            tabs[n].slider(
                label="Select the Zscore cutoff:",
                min_value=float(df['zscore'].min()),
                max_value=float(df['zscore'].max()),
                value=st.session_state.get(slider_key, defaultval),
                key=slider_key
            )
            tabs[n].button("Auto Split", on_click=update_slider, kwargs={"n":key, "sliderkey":slider_key}, key=f"autosplit_{key}")
            # Filtered DataFrame based on slider value
            cutoff = st.session_state[slider_key]
            filtered_df = df[df['zscore'] >= cutoff]

            # create degernerate represenations
            dseqs = degenerate.seqs_to_degenerates(filtered_df['kmers'].tolist())
            # Plot histogram with vertical line at cutoff
            althist = histogram_plot(df, maxbins=100, cutoff=cutoff)
            tabs[n].col1.altair_chart(althist)

            tabs[n].subheader(f"Data for length {key}:")
            tabs[n].write(filtered_df)

            tabs[n].markdown("## Degenerate PAM/TAM sequences")
            tabs[n].dataframe({"PAM/TAM site":dseqs}, hide_index=True)
            tabs[n].markdown('## Sequence motif')
            logo = st.session_state.pamexpobj.make_logo(length=key, cutoff=cutoff, type='zscore',above=False)
            tabs[n].pyplot(logo)


if __name__ == "__main__":
    main(args)


