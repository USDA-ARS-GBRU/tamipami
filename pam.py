# -*- coding: utf-8 -*-
# #!/usr/bin/local python
"""pam:  a TamiPami module for working with pam spacer sequence 
"""

import math
import statistics 
import logging


import scipy
import pandas as pd
import matplotlib.pyplot as plt
from skbio.stats import composition
import ckmeans

import logomaker

class pamSeqExp():
    """ A class for saving and processing kmer cound data from PAM/TAM library sequencing experiments
    """
    def __init__(self, ctl: dict[str, int], exp: dict[str, int], position: str) -> None:
        """ 
            Args:
                ctl: a single count dictionary of kmers and their counts from the uncut control sequencing library
                exp: a single count dictionary of kmers and their counts fro mthe library cut with the endonuclease
                position: positional information on the orientation of the pam relative to the guide site: 3prime or 5prime
            Returns:
                None
        """
        self.ctl = ctl
        self.exp = exp
        self.position = position
        self.multikmerdict =  None
        self._create_multilevel()

    def _kmer_len(self, kmers:dict[str, int]) -> int:
        """Verify the length of all kmers in a dict
            Args:
                kmers:  a dict of kmers and their counts
            Returns:
                the length of the kmers:
            Raises:
                AssertionError if the linngth of all kmers is not the same.
        """
        kmerlen: int = None
        for kmer in kmers.keys():
            if not kmerlen:
                kmerlen = len(kmer)
            else:
                assert kmerlen == len(kmer)
        return kmerlen
    
    def _kmersummary(self, kmer_dict: dict[str,int]) -> dict[str,list]:
        """Creates a dict of kmers, counts and, the centered log ratio of the data as three lists
            Args:
                kmer_dict:  kmers:  a dict of kmers and their counts
            Retunrs: 
                dict of kmers, counts and, the centered log ratio of the data as three lists
        """
        kmers = kmer_dict.keys()
        counts = list(kmer_dict.values())
        assert sum(counts) > 0, "No targets were identified in your data, Please check your Library, spacer ,and orientation settings"
        refmc = composition.multi_replace(counts)
        clr = composition.clr(refmc).tolist()
        return {"kmers": kmers, "counts": counts, "clr": clr}
    
    def _combine_single_pair(self, exper: dict[str:list], ctl: dict[str,list]) -> pd.DataFrame:
        """Takes the control and experimental count and clr data and substracts the clrs, 
           and estimates zscores and signifigance, returning results as a dataframe.

        Args:
            exper: the experimental dict
            ctl: the control dict
        Returns: 
            a python data frame with kmers, ctl_raw, exp_raw, ctl_clr, exp_clr, diff, zscore, pvalue and BH adjusted p-value

        """
        assert exper['kmers'] == ctl['kmers'], "Kmers from the experimental and control data do not match."

        data = {"kmers": ctl['kmers'],
                "ctl_raw": ctl['counts'],
                "exp_raw": exper['counts'],
                 "ctl_clr": ctl['clr'],
                 "exp_clr": exper['clr']}
        df = pd.DataFrame.from_dict(data)
        df["diff"] = df["ctl_clr"] - df["exp_clr"]
        df["zscore"] = (df['diff'] - df['diff'].mean())/ df['diff'].std()
        df["pvalue"] = scipy.stats.norm.pdf(df['zscore'])
        df["p_adjust_BH"] = scipy.stats.false_discovery_control(df["pvalue"])
        df = df.sort_values(by=["kmers"])
        return df

    def _group_and_sum_kmers(self, df: pd.DataFrame, N: int, position: str ='3prime') -> tuple[dict[str,int],dict[str,int]]:
        """This function takes a df with longger kmers and sums up the counts at each shorter kmer. This is done so that 
            different lengths can be examines but the  reads do not need to be read multiple times.

            Args:
                df: The summary dataframe.
                N: the length of kmer to summarize too (must be smaller than the length of the input data )
                position: 5prime if the orientation is 5'-PAM-guide-3' or 3prime if the orientation is 5'-guide-PAM-3'
            Returns:
                A tuple of dicts 
        """
        if 'kmers' not in df.columns or 'ctl_raw' not in df.columns or 'exp_raw' not in df.columns:
            raise ValueError("DataFrame must contain 'kmers', 'ctl_raw', and 'exp_raw' columns")
        if not all(isinstance(x, int) for x in df['ctl_raw']) or not all(isinstance(x, int) for x in df['exp_raw']):
            raise ValueError("'ctl_raw' and 'exp_raw' columns must contain positive integers")
        if N <= 0:
            raise ValueError("N must be a positive integer")
        if position not in ['3prime', '5prime']:
            raise ValueError("position must be '3prime' or '5prime'")
        if position == '3prime':
            df['grouped_kmer'] = df['kmers'].str[:N]
        elif position == '5prime':
            df['grouped_kmer'] = df['kmers'].str[-N:]    
        grouped_df = df.groupby('grouped_kmer').agg({'ctl_raw': 'sum', 'exp_raw': 'sum'}).reset_index()
        grouped_df.rename(columns={"grouped_kmer": "kmers"}, inplace=True)
        return grouped_df.set_index('kmers')['ctl_raw'].to_dict(), grouped_df.set_index('kmers')['exp_raw'].to_dict()

    def _create_multilevel(self, mink: int = 3) -> None:
        multikmerdict = {}
        kl = self._kmer_len(self.ctl)
        ctlsum = self._kmersummary(self.ctl)
        expsum = self._kmersummary(self.exp)
        multikmerdict[kl] = self._combine_single_pair(exper=expsum, ctl=ctlsum)
        while kl > mink:
            short_ctl, short_exp = self._group_and_sum_kmers(df=multikmerdict[kl], N=kl - 1, position=self.position)
            multikmerdict[kl -1 ] = self._combine_single_pair(self._kmersummary(short_exp), self._kmersummary(short_ctl))
            kl = kl - 1
        self.multikmerdict = multikmerdict

    def find_breakpoint(self, length, type='zscore'):
        assert type in ['zscore','diff'], "Parameter 'type' should be 'diff' or 'zscore'"
        data = self.multikmerdict[int(length)][type].to_numpy()
        return ckmeans.breaks(data, 2)[0]
    
    def check_N(self, n: int) -> tuple[float]:
        """ Estimate the fraction of SD expected to come from shot noise
        The fraction of the SD expected to be attributable to Poisson or shot noise.
        A large value indicates that more reads or a shorter length value are required.  
        Args:
        vect (list): A list of kmer counts  
        Returns:
        a tuple of two floats showing  he fraction of SD attributable Poisson Noise in the control and experimental data sets
        """
        def check(vect):
            if vect.empty:
                return 0.0
            return math.sqrt(statistics.mean(vect)) / statistics.stdev(vect)
        ctl_raw_sd = check(self.multikmerdict[n]["ctl_raw"])
        exp_raw_sd = check(self.multikmerdict[n]["exp_raw"])
        return ctl_raw_sd, exp_raw_sd
    
    def make_logo(self, length: int, cutoff: float, type: str='zscore', above: bool=True) -> None:
        """
        Takes a pandas dataframe and saves a Sequence motif logo of the signifigant

        Args:
            length: the length of the kmers
            cutoff: the threshold for selected column to include in the motif
            filename:  a path for the pdf, jpg or png of figure
        Returns:
            None
        """
        try:
            df = self.multikmerdict[length]
            if above:
                df_filtered = df[df[type] < cutoff]
            else:
                df_filtered = df[df[type] >= cutoff]
            prob_df = logomaker.alignment_to_matrix(sequences=df_filtered['kmers'], to_type = 'probability', pseudocount = 0)
            logo_fig = logomaker.Logo(prob_df, color_scheme='classic')
            #plt.savefig(filename)
            return plt
        except Exception as e:
            logging.warning( "could not generate sequence motif graphic")
            raise e