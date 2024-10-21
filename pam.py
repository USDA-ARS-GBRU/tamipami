"""A class for working with pam spacer sequence counts using Padas dataframes
"""

import pandas as pd
from skbio.stats import composition
import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
from ckmeans_1d_dp import ckmeans

class pamSeqExp():
    def __init__(self, ctl: dict[str, int], exp: dict[str, int], position: str) -> None:
         self.ctl = ctl
         self.exp = exp
         self.position = position
         self.multikmerdict = None

    def _kmer_len(self, kmers:dict[str, int]) -> int:
        kmerlen: int = None
        for kmer in kmers.keys():
            if not kmerlen:
                kmerlen = len(kmer)
            else:
                assert kmerlen == len(kmer)
        return kmerlen
    
    def _kmersummary(self, kmer_dict: dict[str,int]) -> dict[str,list]:
        kmers = kmer_dict.keys()
        counts = list(kmer_dict.values())
        refmc = composition.multi_replace(counts)
        clr = composition.clr(refmc).tolist()
        return {"kmers": kmers, "counts": counts, "clr": clr}
    
    def _combine_single_pair(self, exper: dict[str:list], ctl: dict[str,list]) -> pd.DataFrame:
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

    def _create_multilevel(self, ) -> dict[int,pd.DataFrame]:
        multikmerdict = {}
        kl = self._kmer_len(self.ctl)
        ctlsum = self._kmersummary(self.ctl)
        expsum = self._kmersummary(self.exp)
        multikmerdict[kl] = self._combine_single_pair(exper=expsum, ctl=ctlsum)
        print(multikmerdict[kl])
        while kl > 1:
            short_ctl, short_exp = self._group_and_sum_kmers(df=multikmerdict[kl], N=kl - 1, position=self.position)
            multikmerdict[kl -1 ] = self._combine_single_pair(self._kmersummary(short_ctl), self._kmersummary(short_exp))
            kl = kl - 1
        self.multikmerdict = multikmerdict

    def plot_kmer_zscores(self, attribute='zscore', save_path=None):
        data_dict = self.multikmerdict
        num_plots = len(data_dict)
        grid_size = int(num_plots**0.5 + 0.5)  # Calculate the square grid size
        fig, axes = plt.subplots(grid_size, grid_size, figsize=(15, 15))
        axes = axes.flatten()  # Flatten the grid of axes
        assert attribute in ["zscore", "pvalue", "p_adjust_BH", "diff", "ctl_raw","ctl_clr","exp_raw","exp_clr"],"values of 'attribute' must be in ['zscore', 'pvalue', 'p_adjust_BH', 'diff', 'ctl_raw','ctl_clr','exp_raw','exp_clr']"
        for i, (kmer_length, df) in enumerate(data_dict.items()):
            sns.histplot(df[attribute], kde=True, ax=axes[i])
            axes[i].set_title(f'Kmer Length: {kmer_length}')        
        # Remove any empty subplots if num_plots is not a perfect square
        for j in range(i + 1, grid_size * grid_size):
            fig.delaxes(axes[j])       
        plt.tight_layout()      
        if save_path:
            plt.savefig(save_path, format='pdf')      
        plt.show()
    
    def find_breakpoint(self, length, type='zscore'):
        assert type in ['zscore','diff'], "Parameter 'type' should be 'diff' or 'zscore'"
        return ckmeans.breaks(self.multikmerdict[length][type], 2)[0]
    
    def select_pam(self, zscore: float = -3.0 ):
        longestkmer = max(list(self.combined_multikmerdict.keys()))
        if self.position='3prime'
            sigkmers = []
            for kmer, data in combined_multikmerdict[longestkmer].items():
                if data['zscore'] < zscore:
                    sigkmers.append(kmer)
    
    
            
a = pamSeqExp(ctl=cont_raw, exp=exp_raw,position='3prime')
a._create_multilevel()
a.plot_kmer_zscores()