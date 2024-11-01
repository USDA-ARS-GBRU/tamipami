import streamlit as st
import numpy as np
import pandas as pd
import os

st.logo("USDAARSIdentityRGB3.png", size= "large")
apptitle = 'TamiPami'
st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")

st.sidebar.markdown("## Select Parameters to identiy the PAM/TAM")


DATA_DIR = os.path.dirname(os.path.abspath(__file__))  # This is your Project Root

   # Define input parameters and widgets


st.title('TamiPami: A web application to identify PAM or TAM sites in new Cas enzymes based on sequencing experiments')

cont1 = st.sidebar.file_uploader('A forward .fastq, .fq, .fastq.gz or .fq.gz file.', key='cont1' )
cont2 = st.sidebar.file_uploader('A reverse .fastq, .fq, .fastq.gz or .fq.gz file.', key='cont2')
exp1 = st.sidebar.file_uploader('A forward .fastq, .fq, .fastq.gz or .fq.gz file.', key='exp1')
exp2 = st.sidebar.file_uploader('A reverse .fastq, .fq, .fastq.gz or .fq.gz file.', key='exp2',)
library = st.sidebar.selectbox('The Addgene library pool. For custom pools use the --spacer and --orientation flags', ["RTW554", "RTW555", "RTW572", "RTW574"])
spacer = st.sidebar.text_input('The spacer sequence for the guide RNA. Not needed if ---library is used', key='spacer')
orientation = st.sidebar.selectbox('The side of the spacer the PAM/TAM is on', ["3prime", "5prime"], key='orientation')
log = st.sidebar.text_input('Log file', value='tamipami.log', key='log')
length = st.sidebar.select_slider('The MAxumum length of the PAM or TAM sequences', options=list(range(1, 11)), value=4, key='length')

if st.sidebar.button("Run"):
    st.write(f"""
    ### Configuration:
    - Forward file: {cont1}
    - Reverse file: {cont2}
    - Experiment forward file: {exp1}
    - Experiment reverse file: {exp2}
    - Library: {library}
    - Spacer: {spacer}
    - Orientation: {orientation}
    - Log file: {log}
    - Length: {length}
    """)

#@st.cached_data
def process_inputs(args):
    pass        