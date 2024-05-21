from typing import IO, Union

import pandas as pd
import numpy as np
from scipy import interpolate

from .annotator_base import check_columns


class KinaseLibraryAnnotator:
    """Annotate pandas dataframe with highest scoring kinases from the kinase library.

    Johnson et al. 2023, https://doi.org/10.1038/s41586-022-05575-3

    Requires "Site sequence context" column in the dataframe to be present.
    The "Site sequence context" column can be generated with PeptidePositionAnnotator().

    Example:
        ::

            annotator = KinaseLibraryAnnotator(<path_to_motifs_file>, <path_to_quantiles_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(
        self,
        motifs_file: Union[str, IO],
        quantiles_file: Union[str, IO],
        score_cutoff: float = 3,
    ):
        """
        Initialize the input files and options for MotifAnnotator.

        Args:
            annotation_file: tab separated file with motifs and their identifiers

        """
        self.motifs_file = motifs_file
        self.quantiles_file = quantiles_file
        self.score_cutoff = score_cutoff
        self.odds_dict = None
        self.quantiles = None

    def load_annotations(self) -> None:
        """Reads in tab separated file with motif and quantile annotations."""
        odds_df = pd.read_csv(
            self.motifs_file, sep="\t", index_col=["Kinase", "Position", "AA"]
        )
        self.odds_dict = odds_df["Odds Ratio"].to_dict()

        quantile_matrix_df = pd.read_csv(
            self.quantiles_file, sep="\t", index_col="Score"
        ).T
        self.quantiles = {}
        for kinase, q in quantile_matrix_df.iterrows():
            self.quantiles[kinase] = _build_ecdf(scores=q.index, quantiles=q.values)

    @check_columns(["Site sequence context"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds column with motifs the site sequence context matches with.

        Adds the following annotation columns to dataframe\:
        
        - Motif Kinases = semicolon separated list of kinases that match with the site sequence contexts
        - Motif Scores = semicolon separated list of scores corresponding to Motif Kinases
        - Motif Percentiles = semicolon separated list of percentiles corresponding to Motif Kinases
        - Motif Totals = semicolon separated list of score*percentile corresponding to Motif Kinases

        Args:
            df: pandas dataframe with "Site sequence context" column
            inplace: Whether to modify the DataFrame rather than creating a new one.

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        site_sequence_plus_minus_5 = (
            annotated_df["Site sequence context"]
            .str.slice(start=10, stop=15)
            .str.upper()
            + annotated_df["Site sequence context"].str.slice(start=15, stop=16)
            + annotated_df["Site sequence context"]
            .str.slice(start=16, stop=21)
            .str.upper()
        )

        annotated_df[
            ["Motif Kinases", "Motif Scores", "Motif Percentiles", "Motif Totals"]
        ] = site_sequence_plus_minus_5.apply(
            lambda x: _find_upstream_kinase(
                x, self.quantiles, self.odds_dict, self.score_cutoff
            )
        )

        return annotated_df


def _score(seq, kinase, odds_dict, motif_size=5, sig_digits=4):
    """
    Score the motif based on the AA positional ODDS matrix given a sequence and a kinase
    """
    score = []
    for i, aa in enumerate(seq):
        pos = i - motif_size
        score.append(odds_dict.get((kinase, pos, aa), np.nan))
    return round(np.log2(np.nanprod(score)), sig_digits)


def _build_ecdf(scores, quantiles):
    """
    Build an empirical cumulative distribution function (ecdf) based on precomputed scores and quantiles

    Parameter
    ---------
    scores : array like
        input scores
    quantiles : array like
        input qunatiles

    Return
    ------
    ecdf function that maps x->cumprop
    """
    # Return the ecdf as function
    f = interpolate.interp1d(scores, quantiles, kind="linear")
    return f


def _find_upstream_kinase(seq, quantiles, odds_dict, score_cutoff=2):
    """
    Score all kinases against in Q and P to the input sequence
    """
    # Fast return for pY
    if len(seq) == 0:
        return pd.Series(["", "", "", ""])

    # score all kinases
    result = {}
    for kinase in quantiles.keys():
        s = _score(seq, kinase, odds_dict, motif_size=5, sig_digits=3)
        q = round(float(quantiles[kinase](s)), 3)
        t = round(s * q, 3)
        result[kinase] = (s, q, t)

    # prepare the output by sorting from high to low and report the top 5 hits
    top_5 = sorted(result.items(), key=lambda item: item[1][2], reverse=True)[:5]
    top_5 = [(k, tpl) for k, tpl in top_5 if tpl[2] > score_cutoff]
    if len(top_5) == 0:
        return pd.Series(["", "", "", ""])

    kinases, scores, percentiles, totals = zip(
        *[(k, s, q, t) for k, (s, q, t) in top_5]
    )

    kinases = ";".join(kinases)
    scores = ";".join(map(str, scores))
    percentiles = ";".join(map(str, percentiles))
    totals = ";".join(map(str, totals))

    return pd.Series([kinases, scores, percentiles, totals])
