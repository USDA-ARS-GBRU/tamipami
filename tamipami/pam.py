# -*- coding: utf-8 -*-
"""pam:  a TamiPami module for working with PAM/TAM spacer sequences"""


import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import ckmeans
import logomaker
from skbio.stats import composition
from scipy.stats import gaussian_kde

logger = logging.getLogger(__name__)


class pamSeqExp:
    """A class for managing and analyzing kmer count data from PAM/TAM sequencing experiments.

    This class provides methods to process kmer data from control and experimental sequencing libraries,
    calculate kmer lengths, summarize kmer data, combine control and experimental data, group and sum kmers,
    and create multilevel kmer dictionaries. It also includes functionality to find breakpoints, check noise levels,
    and generate sequence motif logos.
    """

    def __init__(
        self,
        ctl: dict[str, int] = None,
        exp: dict[str, int] = None,
        position: str = None,
        multikmerdict: Dict[int, pd.DataFrame] = None,
    ) -> None:
        """
        Initializes the pamSeqExp object with kmer count data and positional information.

        Args:
            ctl: Dictionary of kmers and their counts from the uncut control sequencing library.
            exp: Dictionary of kmers and their counts from the library cut with the endonuclease.
            position: Orientation of the PAM relative to the guide site: '3prime' or '5prime'.
            multikmerdict: Optional precomputed dict[int -> DataFrame] containing multilevel kmer data.
        """
        self.ctl = ctl
        self.exp = exp
        self.position = position

        # Do NOT build in __init__; Streamlit hot-reloader can cause timing issues.
        # Use lazy initialization instead.
        if multikmerdict is not None and not isinstance(multikmerdict, dict):
            raise TypeError("multikmerdict must be a dict[int -> pd.DataFrame].")
        self._multikmerdict: Optional[Dict[int, pd.DataFrame]] = multikmerdict

        # Optional knobs (can be set by caller)
        self.alpha0: float = getattr(self, "alpha0", 0.5)              # Jeffreys prior
        self.dirichlet_draws: int = getattr(self, "dirichlet_draws", 8000)
        # Fixed default seed so identical inputs give identical results; set to None for a random seed
        self.random_seed: Optional[int] = getattr(self, "random_seed", 0)

    # ---------- Lazy initialization helpers ----------
    def ensure_multikmerdict(self, mink: int = 3) -> Dict[int, pd.DataFrame]:
        """Build multikmerdict if it hasn't been built yet."""
        if self._multikmerdict is None:
            self._create_multilevel(mink=mink)
        return self._multikmerdict

    @property
    def multikmerdict(self) -> Dict[int, pd.DataFrame]:
        """Property that lazily builds multikmerdict on first access."""
        return self.ensure_multikmerdict()

    @multikmerdict.setter
    def multikmerdict(self, value: Dict[int, pd.DataFrame]) -> None:
        self._multikmerdict = value

    # ---------- Core utilities ----------
    def _kmer_len(self, kmers: dict[str, int]) -> int:
        """Verify the uniformity of kmer lengths in a dictionary."""
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
        """
        kmers = list(kmer_dict.keys())
        counts = list(kmer_dict.values())

        if sum(counts) <= 0:
            logger.error(
                "Sum of counts is zero or negative. Check Library, spacer, and orientation settings."
            )
            raise ValueError(
                "No targets were identified in your data. Please check your Library, spacer, and orientation settings."
            )

        if any(c < 0 for c in counts):
            logger.error("Counts list contains negative integers.")
            raise ValueError(
                "Counts must be a non-empty list of non-negative integers."
            )

        try:
            refmc = composition.multi_replace(counts)
            clr = composition.clr(refmc).tolist()
        except Exception as e:
            logger.error(f"Error in processing kmer data: {e}")
            raise ValueError(f"Error in processing kmer data: {e}")

        return {"kmers": kmers, "counts": counts, "clr": clr}

    def _combine_single_pair(
        self, exper: dict[str, List], ctl: dict[str, List]
    ) -> pd.DataFrame:
        """
        Empirical-null compositional depletion test with finite-depth uncertainty
        propagation using a **paired baseline Dirichlet** sampler + **precision-weighted
        left-fold sigma** for the empirical null.

        Inputs MUST be outputs of _kmersummary:
            { "kmers": [...], "counts": [...], "clr": [...] }

        Steps:
        1) Paired Dirichlet: sample a shared baseline b from pooled counts,
            then sample control and experimental compositions conditioned on b.
        2) CLR per posterior draw; Diff = clr(control) - clr(experimental)
        3) Empirical null: HSM(diff_mean), sigma via left-fold reflection
            with inverse-variance weights (precision-weighted)
        4) Per-component tests: z_j = (Δ̂_j - μ_L) / sqrt(σ_L^2 + SE_j^2)
        5) BH-FDR
        6) 95% credible intervals added (ci_lo, ci_hi)
        """
        # ----- Validate alignment -----
        kmers_exp = np.asarray(exper["kmers"])
        kmers_ctl = np.asarray(ctl["kmers"])
        if not np.array_equal(kmers_exp, kmers_ctl):
            raise ValueError("Experimental and control kmers do not match.")

        kmers = kmers_ctl
        x_ctl = np.asarray(ctl["counts"], dtype=float)
        x_exp = np.asarray(exper["counts"], dtype=float)
        ctl_clr = np.asarray(ctl["clr"], dtype=float)
        exp_clr = np.asarray(exper["clr"], dtype=float)
        D = kmers.size

        # ----- Paired baseline Dirichlet sampler -----
        # Shared baseline prior from pooled counts (empirical Bayes)
        alpha0 = getattr(self, "alpha0", 0.5)
        M      = getattr(self, "dirichlet_draws", 8000)
        seed   = getattr(self, "random_seed", 0)
        rng    = np.random.default_rng(seed)

        # Pooled parameters (shared library bias)
        alpha_pool = x_ctl + x_exp + alpha0

        # Observed CLR difference -> componentwise ratio r (mean shift for experimental)
        # This uses your already-computed CLRs to define the experimental mean perturbation.
        diff_obs = exp_clr - ctl_clr
        r = np.exp(diff_obs)  # element-wise; experimental mean ~ b * r (renormalized)

        # Dirichlet concentration anchored by depths (allow optional scaling)
        k_ctl = float(np.sum(x_ctl))
        k_exp = float(np.sum(x_exp))
        kappa_scale = float(getattr(self, "kappa_scale", 1.0))  # e.g., 0.8–1.0 if you want mild overdispersion
        k_ctl_eff = max(1.0, kappa_scale * k_ctl)
        k_exp_eff = max(1.0, kappa_scale * k_exp)

        # (Optional) adaptive cap to avoid huge M*D in CI
        mc_cap = int(getattr(self, "mc_work_cap", 2_000_000))  # cap of total samples M*D
        M_eff  = min(M, max(1000, mc_cap // max(D, 1)))

        # Floor for Dirichlet draws: very small concentrations can underflow to exactly 0,
        # and log(0) would turn the whole CLR vector into NaN
        p_floor = np.finfo(float).tiny

        # Paired sampling loop
        deltas = np.empty((M_eff, D), dtype=float)
        for m in range(M_eff):
            # 1) Shared baseline draw
            b = rng.dirichlet(alpha_pool)

            # 2) Control draw around baseline
            p_c = rng.dirichlet(k_ctl_eff * b)

            # 3) Experimental mean: baseline * ratio, renormalized; then draw
            mean_e = b * r
            mean_e = mean_e / mean_e.sum()
            p_e = rng.dirichlet(k_exp_eff * mean_e)

            # 4) CLR and paired difference
            p_c = np.maximum(p_c, p_floor)
            p_e = np.maximum(p_e, p_floor)
            log_pc = np.log(p_c); clr_c = log_pc - log_pc.mean()
            log_pe = np.log(p_e); clr_e = log_pe - log_pe.mean()
            deltas[m] = clr_c - clr_e

        # Posterior summaries of the paired difference
        diff_mean = deltas.mean(axis=0)
        diff_se   = deltas.std(axis=0, ddof=1)
        ci_lo, ci_hi = np.quantile(deltas, [0.025, 0.975], axis=0)

        # ----- Empirical null: HSM on diff_mean -----
        def _half_sample_mode(x: np.ndarray) -> float:
            xs = np.sort(x)
            while xs.size > 3:
                n = xs.size
                w = int(np.ceil(n / 2))
                spans = xs[w:] - xs[:-w]
                i = int(np.argmin(spans))
                xs = xs[i : i + w]
            return float(np.median(xs))

        mu_left = _half_sample_mode(diff_mean)

        # ----- Precision-weighted left-fold sigma -----
        # Left offsets relative to mu_left; weights from posterior precision 1/se^2
        left_mask    = (diff_mean <= mu_left)
        left_offsets = diff_mean[left_mask] - mu_left
        diff_se_left = diff_se[left_mask]

        if left_offsets.size >= 5 and np.any(left_offsets != 0):
            eps   = float(getattr(self, "null_eps", 1e-12))            # stability
            w     = 1.0 / (diff_se_left**2 + eps)                      # inverse-variance weights
            qclip = float(getattr(self, "null_weight_clip_q", 0.95))   # clip extreme weights
            w_cap = np.quantile(w, qclip)
            w     = np.minimum(w, w_cap)

            # Reflect for symmetry; compute weighted SD
            reflected = np.concatenate([left_offsets, -left_offsets])
            w_sym     = np.concatenate([w, w])
            mean_w    = np.average(reflected, weights=w_sym)
            sigma_left = float(np.sqrt(np.average((reflected - mean_w)**2, weights=w_sym)))
        else:
            # Robust fallback: MAD -> sigma on diff_mean
            mad = float(np.median(np.abs(diff_mean - mu_left)))
            sigma_left = 1.4826 * mad

        sigma_left = max(sigma_left, 1e-8)

        # ----- Component-wise tests (right-tailed) -----
        var = sigma_left**2 + diff_se**2
        z   = (diff_mean - mu_left) / np.sqrt(var)
        p   = scipy.stats.norm.sf(z)

        # BH-FDR
        def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
            m = pvals.size
            order = np.argsort(pvals)
            ranked = pvals[order]
            bh = ranked * m / np.arange(1, m + 1)
            q_ranked = np.minimum.accumulate(bh[::-1])[::-1]
            q_ranked = np.minimum(q_ranked, 1.0)
            q = np.empty_like(q_ranked)
            q[order] = q_ranked
            return q

        q = _bh_adjust(p)
        is_cut = (q <= 0.05) & (z > 0)

        # ----- Effect sizes (centered at empirical null) -----
        corrected       = diff_mean - mu_left
        true_survival   = np.exp(-corrected)
        raw_depletion   = (1.0 - true_survival) * 100.0

        # ----- Assemble output DataFrame with final schema names -----
        df = pd.DataFrame({
            "kmers": kmers,
            "ctl_raw": x_ctl,
            "exp_raw": x_exp,

            # Observed CLR (from _kmersummary)
            "ctl_clr": ctl_clr,
            "exp_clr": exp_clr,

            # Effect & uncertainty from paired sampler
            "diff":   diff_mean,
            "diff_se": diff_se,
            "ci_lo":  ci_lo,
            "ci_hi":  ci_hi,

            # Empirical null params (scalars repeated for convenience)
            "mu_left":   mu_left,
            "sigma_left": sigma_left,

            # Testing & FDR
            "zscore": z,
            "pvalue": p,
            "p_adjust_BH": q,
            "is_statistically_cut": is_cut,

            # Effect size (% depleted), centered at empirical null
            "true_percent_depleted": np.clip(raw_depletion, 0, 100),
        })

        df["shrunk_percent_depleted"] = np.where(
            df["is_statistically_cut"], df["true_percent_depleted"], 0.0
        )

        return df.sort_values("kmers")

    # ---------- Group & sum kmers ----------


    def _group_and_sum_kmers(
        self, df: pd.DataFrame, N: int, position: str = "3prime"
    ) -> Tuple[dict[str, int], dict[str, int]]:
        """
        Groups and sums kmer counts from a DataFrame for a specified kmer length.
        Ensures output counts are integers (as required by downstream code).
        """
        required = {"kmers", "ctl_raw", "exp_raw"}
        if not required.issubset(df.columns):
            raise ValueError("DataFrame must contain 'kmers', 'ctl_raw', and 'exp_raw' columns")

        if N <= 0:
            raise ValueError("N must be a positive integer")
        if position not in ["3prime", "5prime"]:
            raise ValueError("position must be '3prime' or '5prime'")

        # Create grouping key
        df = df.copy()
        if position == "3prime":
            df["grouped_kmer"] = df["kmers"].str[:N]
        else:  # "5prime"
            df["grouped_kmer"] = df["kmers"].str[-N:]

        # Group & sum
        grouped_df = (
            df.groupby("grouped_kmer", as_index=False)
            .agg({"ctl_raw": "sum", "exp_raw": "sum"})
        )
        grouped_df.rename(columns={"grouped_kmer": "kmers"}, inplace=True)

        # ---- Robust dtype handling for integer counts ----
        def _ensure_integer_counts(series: pd.Series, name: str) -> pd.Series:
            # If already integer dtype, return as-is
            if pd.api.types.is_integer_dtype(series.dtype):
                return series

            # Convert to numeric (handles object/str), raising on bad values
            series_num = pd.to_numeric(series, errors="raise")

            # Check "integer-like": all values close to nearest integer
            rounded = np.round(series_num)
            if np.allclose(series_num, rounded, rtol=0, atol=1e-8):
                # Safe to cast
                return rounded.astype(int)

            # If not integer-like, fail loudly
            raise ValueError(
                f"Column '{name}' contains non-integer values after grouping. "
                f"Got dtype={series.dtype}, example values={series_num.iloc[:5].tolist()}"
            )

        grouped_df["ctl_raw"] = _ensure_integer_counts(grouped_df["ctl_raw"], "ctl_raw")
        grouped_df["exp_raw"] = _ensure_integer_counts(grouped_df["exp_raw"], "exp_raw")

        # Return dicts
        ctl_dict = grouped_df.set_index("kmers")["ctl_raw"].to_dict()
        exp_dict = grouped_df.set_index("kmers")["exp_raw"].to_dict()
        return ctl_dict, exp_dict


    # ---------- Create multilevel dict ----------
    def _create_multilevel(self, mink: int = 3) -> None:
        """Build multikmerdict: dict[int -> DataFrame] across k-mer lengths."""
        if self.ctl is None or self.exp is None or self.position is None:
            raise ValueError("ctl, exp, and position must be provided to build multilevel data.")

        multikmerdict: Dict[int, pd.DataFrame] = {}

        # Longest k-mer level
        kl = self._kmer_len(self.ctl)
        ctlsum = self._kmersummary(self.ctl)
        expsum = self._kmersummary(self.exp)
        multikmerdict[kl] = self._combine_single_pair(exper=expsum, ctl=ctlsum)

        # Shorter levels down to mink
        while kl > mink:
            short_ctl, short_exp = self._group_and_sum_kmers(
                df=multikmerdict[kl], N=kl - 1, position=self.position
            )
            ctlsum_short = self._kmersummary(short_ctl)
            expsum_short = self._kmersummary(short_exp)
            multikmerdict[kl - 1] = self._combine_single_pair(
                exper=expsum_short, ctl=ctlsum_short
            )
            kl -= 1

        # Store in private attribute (do NOT build in __init__)
        self._multikmerdict = multikmerdict

    # ---------- Breakpoint / diagnostics / logo ----------


    def find_breakpoint(self, length: int, type: str = "zscore") -> float:
        """
        Auto-split cutoff for a given k-mer length and score type ("diff" or "zscore").

        Strategy:
        1) Clean data (finite values).
        2) Get two clusters via ckmeans (k=2).
        3) Find KDE valley between cluster means (minimum density).
        4) Guardrails: for zscore, enforce non-negative cutoff; constrain to sensible percentiles.
        5) Fallbacks: μ_left + 2σ_left or 80th percentile when bimodality is weak.

        Returns:
        float cutoff value suitable for isolating the right tail (positives).
        """
        if type not in ["zscore", "diff"]:
            raise ValueError("Parameter `type` should be 'diff' or 'zscore'")

        df = self.multikmerdict.get(int(length))
        if df is None:
            raise KeyError(f"Length {length} not found in multikmerdict.")

        x = df[type].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        if x.size == 0:
            # No data; conservative default
            return 0.0 if type == "zscore" else float(self.multikmerdict[length]["mu_left"].iloc[0])

        # If the distribution is too small or near-constant, return a robust default
        if x.size < 16 or np.allclose(np.std(x), 0.0, atol=1e-12):
            if type == "zscore":
                return 0.0
            # For diff, a sensible split is μ_left + 2σ_left (repeated scalar per row)
            mu_left = float(self.multikmerdict[length]["mu_left"].iloc[0])
            sigma_left = float(self.multikmerdict[length]["sigma_left"].iloc[0])
            return mu_left + 2.0 * sigma_left

        # --- Cluster with ckmeans (k=2) ---
        # Note: ckmeans expects 1D array; breaks returns boundary between clusters
        try:
            boundary = float(ckmeans.breaks(x, 2)[0])
        except Exception:
            # Fallback if ckmeans fails
            boundary = float(np.median(x))

        # Identify cluster masks based on boundary
        left_mask  = x <= boundary
        right_mask = x > boundary

        # Guard against degenerate clustering
        if left_mask.sum() < 4 or right_mask.sum() < 4:
            # Not reliably bimodal; use percentile-based split
            # For zscore, prefer a non-negative cutoff
            if type == "zscore":
                return max(0.0, float(np.percentile(x, 80)))
            else:
                mu_left = float(self.multikmerdict[length]["mu_left"].iloc[0])
                sigma_left = float(self.multikmerdict[length]["sigma_left"].iloc[0])
                # Slightly more conservative than 80th percentile
                return max(mu_left + 2.0 * sigma_left, float(np.percentile(x, 75)))

        m1 = float(np.mean(x[left_mask]))
        m2 = float(np.mean(x[right_mask]))

        # Ensure m1 <= m2 (swap if necessary)
        if m1 > m2:
            m1, m2 = m2, m1

        # --- KDE valley between cluster means ---
        try:
            kde = gaussian_kde(x)
            # Build grid between the two cluster means; extend slightly to be safe
            pad = 0.05 * (m2 - m1 + 1e-12)
            grid = np.linspace(m1 - pad, m2 + pad, 512)
            dens = kde(grid)
            # Find the valley (minimum density) between the peaks
            idx_min = int(np.argmin(dens))
            cutoff = float(grid[idx_min])
        except Exception:
            # Fallback to boundary when KDE fails
            cutoff = boundary

        # --- Guardrails ---
        # Constrain to sensible percentiles to avoid pathological minima
        p60, p99 = np.percentile(x, [60, 99])
        cutoff = float(np.clip(cutoff, p60, p99))

        # For zscore, enforce non-negative cutoff (we only want the right tail)
        if type == "zscore":
            cutoff = max(0.0, cutoff)

        return cutoff

    def make_logo(
        self,
        length: int,
        cutoff: float,
        score_type: str = "zscore",
        above: bool = True,
        filename: str = None,
    ) -> plt.Figure | None:
        """
        Generates and (optionally) saves a sequence motif logo based on kmer significance.

        Returns:
            matplotlib.figure.Figure if successful; None only if an unexpected error occurs.
        """
        try:
            if length not in self.multikmerdict:
                raise KeyError(f"Length {length} not found in multikmerdict.")
            df = self.multikmerdict[length]

            if score_type not in df.columns:
                raise ValueError(f"score_type '{score_type}' not found in DataFrame columns.")

            # Note: keep your original semantics — above=True selects < cutoff
            if above:
                df_filtered = df[df[score_type] < cutoff]
            else:
                df_filtered = df[df[score_type] >= cutoff]

            # Build probability matrix and logo
            prob_df = logomaker.alignment_to_matrix(
                sequences=df_filtered["kmers"], to_type="probability", pseudocount=0
            )
            logo = logomaker.Logo(prob_df, color_scheme="colorblind_safe")

            # Retrieve the matplotlib Figure from the Logo's axes
            fig = logo.ax.figure

            # Save if requested
            if filename:
                valid_extensions = (".pdf", ".svg", ".eps", ".ps", ".png", ".jpg", ".jpeg", ".tiff")
                if filename.endswith(valid_extensions):
                    fig.savefig(filename, bbox_inches="tight")
                else:
                    logger.warning(
                        "Invalid file extension. Please use one of: %s", ", ".join(valid_extensions)
                    )

            # Always return the figure so the caller (Streamlit/app) can display it
            return fig

        except KeyError as e:
            logger.error("KeyError: Could not find the specified length in multikmerdict.")
            raise
        except ValueError as e:
            logger.error("ValueError: Invalid value encountered during logo generation.")
            raise
        except Exception as e:
            logger.error(
                "Unexpected error occurred during sequence motif graphic generation."
            )
            raise