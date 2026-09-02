# -*- coding: utf-8 -*-
"""pam:  a TamiPami module for working with PAM/TAM spacer sequences"""

import logging
from typing import Dict, List
from typing import Optional
import logging

import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
from skbio.stats import composition
import ckmeans
from sklearn.mixture import GaussianMixture
        

import logomaker


class pamSeqExp:
    """A class for managing and analyzing kmer count data from PAM/TAM sequencing experiments.

    This class provides methods to process kmer data from control and experimental sequencing libraries,
    calculate kmer lengths, summarize kmer data, combine control and experimental data, group and sum kmers,
    and create multilevel kmer dictionaries. It also includes functionality to find breakpoints, check noise levels,
    and generate sequence motif logos.

    Attributes:
        ctl: Dictionary of kmers and their counts from the control library.
        exp: Dictionary of kmers and their counts from the experimental library.
        position: Positional information of the PAM relative to the guide site.
        multikmerdict: Dictionary containing multilevel kmer data.

    Methods:
        _kmer_len: Verifies the length of all kmers in a dictionary.
        _kmersummary: Summarizes kmer data into counts and centered log ratio.
        _combine_single_pair: Combines control and experimental data, calculates z-scores and significance.
        _group_and_sum_kmers: Groups and sums kmer counts for different lengths.
        _create_multilevel: Creates a multilevel dictionary of kmer data.
        find_breakpoint: Finds breakpoints in kmer data based on z-score or difference.
        check_N: Estimates the fraction of standard deviation due to shot noise.
        make_logo: Generates a sequence motif logo based on kmer significance.
    """

    def __init__(
        self,
        ctl: dict[str, int] = None,
        exp: dict[str, int] = None,
        position: str = None,
        multikmerdict: dict[pd.DataFrame] = None,
    ) -> None:
        """
        Initializes the pamSeqExp object with kmer count data and positional information.

        Args:
            ctl: Dictionary of kmers and their counts from the uncut control sequencing library.
            exp: Dictionary of kmers and their counts from the library cut with the endonuclease.
            position: Positional information on the orientation of the PAM relative to the guide site: '3prime' or '5prime'.
            multikmerdict: Precomputed dictionary containing multilevel kmer data.

        Raises:
            ValueError: If neither or both sets of parameters (`ctl`, `exp`, `position` and `multikmerdict`) are provided.

        Returns:
            None
        """
        if not (all([ctl, exp, position]) ^ (multikmerdict is not None)):
            raise ValueError(
                "Please provide either `ctl`, `exp`, and `position` or a `multikmerdict` object."
            )
        if not multikmerdict:
            self.ctl = ctl
            self.exp = exp
            self.position = position
            self.multikmerdict = None
            self._create_multilevel()
        else:
            self.multikmerdict = multikmerdict

    def _kmer_len(self, kmers: dict[str, int]) -> int:
        """Verify the uniformity of kmer lengths in a dictionary.

        Args:
            kmers (dict[str, int]): A dictionary of kmers and their counts.

        Returns:
            int: The length of the kmers.

        Raises:
            ValueError: If the length of all kmers is not the same or if the dictionary is empty.
        """
        kmerlen: Optional[int] = None
        for kmer in kmers.keys():
            if not kmerlen:
                kmerlen = len(kmer)
            else:
                if kmerlen != len(kmer):
                    raise ValueError("All kmers must have the same length.")
        else:
            if kmerlen is None:
                raise ValueError(
                    "The `kmers` dictionary is empty, unable to determine kmer length."
                )
        return kmerlen

    def _kmersummary(self, kmer_dict: dict[str, int]) -> Dict[str, List]:
        """
        Summarizes kmer data into a dictionary containing kmers, their counts,
        and the centered log ratio (CLR) of the data.

        Args:
            kmer_dict (dict[str, int]): A dictionary of kmers and their counts.

        Returns:
            Dict[str, List]: A dictionary with keys 'kmers', 'counts', and 'clr',
            each containing a list of kmers, their counts, and the CLR of the data, respectively.

        Raises:
            ValueError: If the sum of counts is zero or negative, or if any count is negative.
        """

        kmers = list(kmer_dict.keys())
        counts = list(kmer_dict.values())

        if sum(counts) <= 0:
            logging.error(
                "Sum of counts is zero or negative. Check Library, spacer, and orientation settings."
            )
            raise ValueError(
                "No targets were identified in your data, Please check your Library, spacer ,and orientation settings"
            )

        if any(c < 0 for c in counts):
            logging.error("Counts list contains negative integers.")
            raise ValueError(
                "Counts must be a non-empty list of non-negative integers."
            )

        try:
            refmc = composition.multi_replace(counts)
            clr = composition.clr(refmc).tolist()
        except Exception as e:
            logging.error(f"Error in processing kmer data: {e}")
            raise ValueError(f"Error in processing kmer data: {e}")

        return {"kmers": kmers, "counts": counts, "clr": clr}


    def _combine_single_pair(
        self, exper: dict[str, Optional[list]], ctl: dict[str, Optional[list]]
    ) -> pd.DataFrame:
        import numpy as np
        import pandas as pd
        import scipy.stats
        from sklearn.mixture import GaussianMixture

        if exper["kmers"] != ctl["kmers"]:
            raise ValueError("Kmers from the experimental and control data do not match.")

        df = pd.DataFrame({
            "kmers": ctl["kmers"],
            "ctl_raw": ctl["counts"],
            "exp_raw": exper["counts"],
            "ctl_clr": ctl["clr"],
            "exp_clr": exper["clr"],
        })

        # Paired subtraction: positive value means depletion in experimental sample
        df["diff"] = df["ctl_clr"] - df["exp_clr"]

        # --- TWO-COMPONENT DECONVOLUTION ---
        diff_reshaped = df["diff"].values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, covariance_type="diag", random_state=42)
        gmm.fit(diff_reshaped)

        bg_component_idx = np.argmax(gmm.weights_)
        bg_peak_center = gmm.means_.flatten()[bg_component_idx]
        
        # FIX 1: Safely squeeze variances to handle 1D layout vector cleanly
        variances = gmm.covariances_.squeeze()
        
        # FIX 4: Protect against zero-variance component collapse with a floor epsilon
        sigma_diff_noise = max(np.sqrt(variances[bg_component_idx]), 1e-4)

        # --- Z-SCORE & STATISTICAL SIGNIFICANCE ---
        df["zscore"] = (df["diff"] - bg_peak_center) / sigma_diff_noise
        
        # FIX 2: One-tailed (right-tailed) survival function isolates true depletion cuts
        df["pvalue"] = scipy.stats.norm.sf(df["zscore"])
        
        # FIX 3: Prevent asymmetrical tail inflation in FDR multiple-testing controls.
        # We pass a two-tailed p-value matrix to the BH calculator to avoid ranking penalties 
        # from enriched/noise clusters, then apply the cut constraint natively.
        two_tailed_p = 2 * (1 - scipy.stats.norm.cdf(np.abs(df["zscore"])))
        df["p_adjust_BH"] = scipy.stats.false_discovery_control(two_tailed_p)
        
        # An item is a valid cut if it clears the BH threshold AND moves in the depletion direction
        df["is_statistically_cut"] = (df["p_adjust_BH"] <= 0.05) & (df["zscore"] > 0)

        # --- TRUE BIOLOGICAL INFERENCE (EFFECT SIZE) ---
        # FIX 3: Compositional correction applied directly in log-space
        corrected_log_ratio = df["diff"] - bg_peak_center
        true_survival_ratio = np.exp(-corrected_log_ratio)
        
        raw_depletion = (1.0 - true_survival_ratio) * 100.0

        df["true_percent_depleted"] = raw_depletion.clip(0, 100)
        df["shrunk_percent_depleted"] = np.where(
            df["is_statistically_cut"], raw_depletion, 0.0
        ).clip(0, 100)

        return df.sort_values(by="kmers")






    def _group_and_sum_kmers(
        self, df: pd.DataFrame, N: int, position: str = "3prime"
    ) -> tuple[dict[str, int], dict[str, int]]:
        """
        Groups and sums kmer counts from a DataFrame for a specified kmer length.

        This method processes a DataFrame containing kmer data, grouping and summing
        the counts of kmers based on a specified shorter kmer length. It allows for
        analysis of different kmer lengths without re-reading the data.

        Args:
            df (pd.DataFrame): The DataFrame containing kmer data with columns 'kmers',
                'ctl_raw', and 'exp_raw'.
            N (int): The length of the kmer to summarize to, which must be smaller than
                the length of the input kmers.
            position (str, optional): Specifies the orientation of the kmer grouping,
                either '3prime' for 5'-guide-PAM-3' or '5prime' for 5'-PAM-guide-3'.
                Defaults to '3prime'.

        Returns:
            tuple[dict[str, int], dict[str, int]]: A tuple containing two dictionaries
            with grouped kmers as keys and their summed counts as values for control
            and experimental data, respectively.

        Raises:
            ValueError: If the DataFrame does not contain the required columns, if
                'ctl_raw' or 'exp_raw' columns contain non-integer values, if N is not
                a positive integer, or if position is not '3prime' or '5prime'.
        """
        if (
            "kmers" not in df.columns
            or "ctl_raw" not in df.columns
            or "exp_raw" not in df.columns
        ):
            raise ValueError(
                "DataFrame must contain 'kmers', 'ctl_raw', and 'exp_raw' columns"
            )
        if not all(isinstance(x, int) for x in df["ctl_raw"]) or not all(
            isinstance(x, int) for x in df["exp_raw"]
        ):
            raise ValueError(
                "'ctl_raw' and 'exp_raw' columns must contain positive integers"
            )
        if N <= 0:
            raise ValueError("N must be a positive integer")
        if position not in ["3prime", "5prime"]:
            raise ValueError("position must be '3prime' or '5prime'")
        if position == "3prime":
            df["grouped_kmer"] = df["kmers"].str[:N]
        elif position == "5prime":
            df["grouped_kmer"] = df["kmers"].str[-N:]
        grouped_df = (
            df.groupby("grouped_kmer")
            .agg({"ctl_raw": "sum", "exp_raw": "sum"})
            .reset_index()
        )
        grouped_df.rename(columns={"grouped_kmer": "kmers"}, inplace=True)
        return (
            grouped_df.set_index("kmers")["ctl_raw"].to_dict(),
            grouped_df.set_index("kmers")["exp_raw"].to_dict(),
        )

    def _create_multilevel(self, mink: int = 3) -> None:
        multikmerdict = {}
        kl = self._kmer_len(self.ctl)
        ctlsum = self._kmersummary(self.ctl)
        expsum = self._kmersummary(self.exp)
        multikmerdict[kl] = self._combine_single_pair(exper=expsum, ctl=ctlsum)
        while kl > mink:
            short_ctl, short_exp = self._group_and_sum_kmers(
                df=multikmerdict[kl], N=kl - 1, position=self.position
            )
            multikmerdict[kl - 1] = self._combine_single_pair(
                self._kmersummary(short_exp), self._kmersummary(short_ctl)
            )
            kl = kl - 1
        self.multikmerdict = multikmerdict

    def find_breakpoint(self, length: int, type: str = "zscore") -> float:
        """
        Finds the breakpoint in kmer data for a specified length and type.

        This method retrieves kmer data from the multikmer dictionary for a given
        length and calculates the breakpoint using the ckmeans algorithm. The type
        of data used for breakpoint calculation can be either 'zscore' or 'diff'.

        For information on the algorithm see: Haizhou Wang and Mingzhou Song. 2011.
        Ckmeans.1d.dp: Optimal k-means Clustering in One Dimension by Dynamic Programming.
        The R Journal. 3:2, pages 29-33. doi:10.32614/RJ-2011-015.

        Args:
            length (int): The length of the kmer data to analyze.
            type (str, optional): The type of data to use for breakpoint calculation,
                either 'zscore' or 'diff'. Defaults to 'zscore'.

        Returns:
            float: The calculated breakpoint value.

        Raises:
            ValueError: If `type` is not 'zscore' or 'diff', or if `multikmerdict`
                is not initialized.
            KeyError: If the specified length is not found in `multikmerdict`.
        """
        if type not in ["zscore", "diff"]:
            raise ValueError("Parameter `type` should be 'diff' or 'zscore'")
        if self.multikmerdict is None:
            raise ValueError("`multikmerdict` is not initialized.")
        if length not in self.multikmerdict:
            raise KeyError(f"Length {length} not found in multikmerdict.")
        data = self.multikmerdict[int(length)][type].to_numpy()
        return ckmeans.breaks(data, 2)[0]

    def check_N(self, n: int) -> tuple[float]:
        """
        Estimates the fraction of standard deviation due to shot noise for a given kmer length.

        This method calculates the ratio of the square root of the mean to the standard deviation
        for control and experimental raw data vectors associated with a specified kmer length.
        It raises an error if the multikmer dictionary is not initialized or if the specified
        kmer length is not found.

        Args:
            n (int): The kmer length for which to perform the noise check.

        Returns:
            tuple[float]: A tuple containing the calculated noise estimates for control and
            experimental data.

        Raises:
            ValueError: If `multikmerdict` is not initialized.
            KeyError: If the specified kmer length is not found in `multikmerdict`.
        """
        if self.multikmerdict is None:
            raise ValueError("`multikmerdict` is not initialized.")

        def check(vect: List[float]):
            if vect.empty:
                return 0.0
            vect_array = np.array(vect)
            stdev = np.std(vect_array)
            if stdev == 0:
                return 0.0  # or another appropriate default value
            return np.sqrt(np.mean(vect_array)) / stdev

        try:
            ctl_raw_sd = check(self.multikmerdict[n]["ctl_raw"])
            exp_raw_sd = check(self.multikmerdict[n]["exp_raw"])
        except KeyError as e:
            logging.error(
                f"KeyError: {e} - The key {n} does not exist in multikmerdict."
            )
            raise
        return ctl_raw_sd, exp_raw_sd

    def make_logo(
        self,
        length: int,
        cutoff: float,
        score_type: str = "zscore",
        above: bool = True,
        filename: str = None,
    ) -> None:
        """
        Generates and saves a sequence motif logo based on kmer significance.

        This method filters kmer data from the multikmer dictionary based on the
        specified length, cutoff, and score type. It then creates a sequence motif
        logo using the Logomaker library and saves it to a file if a filename is
        provided.

        Args:
            length (int): The length of the kmers to analyze.
            cutoff (float): The threshold for filtering kmers based on the score type.
            score_type (str, optional): The type of score to use for filtering,
                defaults to "zscore".
            above (bool, optional): If True, filters kmers with scores below the cutoff;
                if False, filters kmers with scores above the cutoff. Defaults to True.
            filename (str, optional): The path to save the logo image file. Must end
                with .pdf, .jpg, or .png.

        Raises:
            KeyError: If the specified length is not found in multikmerdict.
            ValueError: If an invalid value is encountered during logo generation.
            Exception: For any unexpected errors during sequence motif graphic generation.

        Returns:
            None
        """
        try:
            if length not in self.multikmerdict:
                raise KeyError(f"Length {length} not found in multikmerdict.")
            df = self.multikmerdict[length]
            if above:
                df_filtered = df[df[score_type] < cutoff]
            else:
                df_filtered = df[df[score_type] >= cutoff]
            prob_df = logomaker.alignment_to_matrix(
                sequences=df_filtered["kmers"], to_type="probability", pseudocount=0
            )
            logo_fig = logomaker.Logo(prob_df, color_scheme="colorblind_safe")
            if filename:
                valid_extensions = (".pdf", ".jpg", ".png")
                if filename.endswith(valid_extensions):
                    plt.savefig(filename)
                else:
                    logging.warning(
                        "Invalid file extension. Please use one of the following: .pdf, .jpg, .png"
                    )
            else:
                return plt
        except KeyError as e:
            logging.error(
                "KeyError: Could not find the specified length in multikmerdict."
            )
            raise e
        except ValueError as e:
            logging.error(
                "ValueError: Invalid value encountered during logo generation."
            )
            raise e
        except Exception as e:
            logging.error(
                "Unexpected error occurred during sequence motif graphic generation."
            )
            raise e
