from typing import IO, Union

import pandas as pd
import numpy as np

from .annotator_base import check_columns

ALLOWED_AA_CHARACTERS = {'A', 'C', 'D', 'E', 'F',
                         'G', 'H', 'I', 'K', 'L',
                         'M', 'N', 'P', 'Q', 'R',
                         'S', 'T', 'V', 'W', 'Y',
                         's', 't', 'y', '_'}


def validate_sequence(seq):
    invalid_chars = set(seq) - ALLOWED_AA_CHARACTERS
    if invalid_chars:
        raise ValueError(f"Sequence '{seq}' contains invalid characters: {''.join(invalid_chars)}")


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
        top_n: int = 5,
        score_cutoff: float = 3,
        split_sequences: bool = False,
        threshold_type="total",
        sort_type="total",
    ):
        """
        Initialize the input files and options for MotifAnnotator.

        Args:
            annotation_file: tab separated file with motifs and their identifiers

        """
        self.motifs_file = motifs_file
        self.quantiles_file = quantiles_file
        self.top_n = top_n
        self.score_cutoff = score_cutoff
        self.threshold_type = threshold_type
        self.split_sequences = split_sequences
        self.sort_type = sort_type
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
            self.quantiles[kinase] = (q.index, q.values)

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

        if self.split_sequences:
            annotated_df['Site sequence context'] = annotated_df['Site sequence context'].apply(lambda s: s.split(';'))
            annotated_df = annotated_df.explode('Site sequence context')

        # Throw an error if any sequence contains illegal characters
        annotated_df["Site sequence context"].apply(validate_sequence)

        def adjust_context_length(sequence, desired_length=11):
            # Edge case: Empty sequences are converted into '_' - they will not receive any kinases or scores
            if sequence == '':
                sequence = '_'

            if len(sequence) > desired_length:
                excess = (len(sequence) - desired_length) // 2
                res = sequence[excess: excess + desired_length]
            elif len(sequence) < desired_length:
                padding_length = (desired_length - len(sequence)) // 2
                res = '_' * padding_length + sequence + '_' * padding_length
            else:
                res = sequence
            return res[:desired_length//2].upper() + res[desired_length//2].lower() + res[(desired_length//2)+1:].upper()

        site_sequence_plus_minus_5 = annotated_df["Site sequence context"].apply(adjust_context_length)

        annotated_df[
            ["Motif Kinases", "Motif Scores", "Motif Percentiles", "Motif Totals"]
        ] = site_sequence_plus_minus_5.apply(
            lambda x: _find_upstream_kinase(
                x,
                self.quantiles,
                self.odds_dict,
                top_n=self.top_n,
                threshold=self.score_cutoff,
                threshold_type=self.threshold_type,
                sort_type=self.sort_type,
            )
        )

        return annotated_df


def _find_upstream_kinase(
    seq, Q, P, top_n=5, threshold=-np.inf, threshold_type="total", sort_type="total"
):
    """
    Score all kinases against input sequence based on Q-Matrix and P-Matrix.
    Percentile is the standard metic according to Johnson et al. and has the best perfromace in my hands as well.

    Input
    -----
    seq : str
        sequence to score
    Q : dict
        Quantile matrix Q
    P : dict
        Probabilty matrix P
    top_n : int
        considers the top n ranked kinases only, default 15.
    threshold: float
        filters for total value bigger than threshold, default -np.inf.
    threshold_type: float
        specifies which metric should be used for filtering (0=score, 1=percentile, 2=score*percentile), default is percentile [1]
    sort_type: float
        specifies which metric should be used for ranking (0=score, 1=percentile, 2=score*percentile), default is percentile [1]

    Returns
    -------
    (kinases, scores, percentiles, totals)
    each is a string with semicolon sorted values
    """
    if len(seq) == 0 or seq == len(seq) * '_':
        return pd.Series(["", "", "", ""])

    # Map the different parameter options
    str_to_int_map = {
        "score": 0,
        "percentile": 1,
        "total": 2,
    }
    if threshold_type not in str_to_int_map:
        raise ValueError("threshold_type")
    if sort_type not in str_to_int_map:
        raise ValueError("sort_type")
    threshold_type = str_to_int_map[threshold_type]
    sort_type = str_to_int_map[sort_type]

    scores = []
    kinases = []
    for kinase in Q.keys():
        s = _score(seq, kinase, P, motif_size=5)
        if s <= 0:
            continue
        scores.append(s)
        kinases.append(kinase)
    scores = np.log2(np.array(scores))

    quantiles = []
    for kinase, s in zip(kinases, scores):
        quantiles.append(_quantile(s, Q[kinase]))
    quantiles = np.array(quantiles)
    totals = scores * quantiles

    result = {k: (s, q, t) for k, s, q, t in zip(kinases, scores, quantiles, totals)}
    # Sort by sort_type, then filter threshold_type > threshold, and then take the topN of the list
    out = sorted(result.items(), key=lambda item: item[1][sort_type], reverse=True)
    out = [(k, *tpl) for k, tpl in out if tpl[threshold_type] > threshold]
    out = out[:top_n]

    if len(out) == 0:
        return pd.Series(["", "", "", ""])

    kinases, scores, quantiles, totals = zip(*out)
    scores = np.round(scores, 3)
    quantiles = np.round(quantiles, 3)
    totals = np.round(totals, 3)

    # Transform to ; separated lists
    kinases = ";".join(kinases)  # kinase
    scores = ";".join(map(str, scores))  # [0] score
    percentiles = ";".join(map(str, quantiles))  # [1] percentile
    totals = ";".join(map(str, totals))  # [2] total = s*q
    return pd.Series([kinases, scores, percentiles, totals])


def _score(seq, kinase, P, motif_size=5):
    """
    Score the motife based on the AA positional ODDS matrix given a sequence and a kinase
    """
    assert len(seq) == (2 * motif_size + 1), \
        f'Sequence expected to have length {(2 * motif_size + 1)} but had length {len(seq)} instead'
    score = 1.0
    for i, aa in enumerate(seq):
        pos = i - motif_size
        score *= P.get((kinase, pos, aa), 1.0)
    return score


def _quantile(s, Q_kinase):
    scores, quantiles = Q_kinase
    index = np.searchsorted(scores, s)
    if index + 1 >= len(scores):
        quantile = quantiles[-1]
    else:
        y1 = quantiles[index]
        y2 = quantiles[index + 1]
        x1 = scores[index]
        x2 = scores[index + 1]
        quantile = y1 + (s - x1) * (y2 - y1) / (x2 - x1)
    return float(quantile)
