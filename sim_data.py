# -*- coding: utf-8 -*-
# #!/usr/bin/local python
"""
A script to generate synthetic data

"""

import numpy as np
from skbio.stats import composition
from scipy.stats import nbinom

import main

def create_sigvect(sig_pams, kmers):
    if sig_pams:
        sigvect = sig_pams.copy()
    else:
        sigvect = {}
    for kmer in kmers:
        if kmer not in sigvect.keys():
            sigvect[kmer] = 0.0
    return dict(sorted(sigvect.items()))


def rnbinom_python(n, mu_vector, size):    
    # Ensure mu_vector is a list or array
    if not isinstance(mu_vector, list):
        mu_vector = [mu_vector]
    # Initialize the output vector
    output_vector = []
    for mu in mu_vector:
        # Calculate the probability of success p for each mu
        p = size / (size + mu)
        # Generate n negative binomially distributed values for each mu
        values = nbinom.rvs(size, p, size=n )#, random_state=rs))
        output_vector.extend(values)    
    return [a.tolist() for a in output_vector]


def make_mu_vect(reads, kmers, log_reduction):
    nullmu = reads/len(kmers)
    return np.exp(list(log_reduction.values()))*nullmu


def make_data(length, sig_pams, reads, size):
    def process(length, sig_pams, reads, size, control):
        kmers = main.iterate_kmer(length)
        if control:
            sig_pams = None
        log_reduction = create_sigvect(sig_pams, kmers)
        mu_vector = make_mu_vect(reads, kmers, log_reduction)
        nbvals = rnbinom_python(n = 1, mu_vector=mu_vector.tolist(), size=size)
        raw_dict = dict(zip(kmers, nbvals))
        mc = composition.multi_replace(list(nbvals))
        dataclr = composition.clr(mc)
        return raw_dict, dataclr
    # Controls
    cont_vals, cont_clr = process(length, sig_pams, reads, size, control=True)
    # experimental
    exp_vals, exp_clr = process(length, sig_pams, reads, size, control=False)
    return cont_vals, cont_clr, exp_vals, exp_clr

#run

reads = 1000000
size = 17
length = 3
sig_pams = {"AGG": -5, "TGG":-5.1, "CGG": -5.2, "GGG":-5.5, 'ATA':-6,'ATT':-6.05,'ATC':-4.75}

cont_raw, cont_clr, exp_raw, exp_clr = make_data(length, sig_pams, reads, size)
df = main.make_df(cont_raw=cont_raw, cont_clr=cont_clr, exp_raw=exp_raw, exp_clr=exp_clr)
print(df)
main.make_logo(df=df, padjust=0.05, filename="sim_logo.pdf" )
df.to_csv("sim_dataframe.csv")