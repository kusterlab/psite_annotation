# adapted from phospho_delocalization.py (Florian P. Bayer - 2025)

import re
import itertools

import pandas as pd
import numpy as np
from scipy.cluster import hierarchy

from .annotator_base import check_columns

PHOSPHORYLATION_PATTERN = re.compile(r"\(ph\)")


class ModifiedSequenceGroupAnnotator:
    """Annotate pandas dataframe with modified sequence groups where localizations are within `match_tolerance` of each other.

    Example:
        ::

            annotator = DelocalizationAnnotator()
            df = annotator.annotate(df)
    """

    def __init__(
        self,
        match_tolerance: int = 2,
    ) -> None:
        """
        Initialize the options for DelocalizationAnnotator.

        Args:
            match_tolerance: group all modifiable positions within n positions of modified sites.

        """
        self.match_tolerance = match_tolerance

    def load_annotations(self) -> None:
        pass

    @check_columns(["Modified sequence"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        r"""Group delocalized phospho-forms.

        This function identifies peptide sequences that differ only by the position
        of their phosphorylation (`(ph)`) group and collapses them into
        "delocalized" groups. Each group contains all modified sequence variants
        that represent the same underlying peptide backbone.

        The following columns are added to the dataframe\:

        - 'Delocalized sequence' = Canonical unmodified backbone with an index
        suffix to distinguish the number of modifications.
        - 'Modified sequence group' = All peptide variants belonging to the same
        delocalized group, concatenated with semicolons.

        Args:
            df: Input dataframe with:
                - `"Modified sequence"` column containing peptide strings with `(ph)` annotations
            inplace: add the new column to df in place

        Returns:
            pd.DataFrame: Dataframe with Modified sequence group column
        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        # Add delocalized sequences.
        annotated_df["Delocalized sequence"] = delocalize_phospho_sequence(
            annotated_df["Modified sequence"]
        )

        # Add modified sequence group clusters.
        annotated_df["Modified sequence group"] = aggregate_phospho_groups(
            annotated_df, self.match_tolerance
        )

        if not inplace:
            return annotated_df


def extract_phos_positions(mod_seq: str, pattern: re.Pattern) -> np.array:
    """
    Parses a modified sequence and reports all positions of the pattern in the aa sequence as numpy array.
    """
    return np.array(
        tuple(
            match.start() - (4 * i) - 1
            for i, match in enumerate(pattern.finditer(mod_seq))
        )
    )


def positional_distance(a: int, b: int) -> int:
    """
    Calculates the positional distance of two position ptm arrays a and b.
    If multiple positions exist, its the maximal distance that defines the distance.
    """
    return max(abs(a - b))


def find_clusters(seqs_pos: list[int], max_distance: int) -> np.array:
    """
    Clusters a group of position ptm arrays if they are closer than max_distance.

    Parameters
    ----------
    seqs_pos : array-like [positions of sequence A, positions of sequence B, ...]
        A list of positional sequences
    max_distance : int >= 0
        The maximal distance two ptm positions can be apart to be considered similar.

    Returns
    -------
    cluster_ids : array-like
        a list of cluster integer ids in the same order as the seqs_pos input list.
    """
    if len(seqs_pos) > 1:
        distance_matrix = [
            positional_distance(a, b) for a, b in itertools.combinations(seqs_pos, 2)
        ]
        linkage_matrix = hierarchy.linkage(
            distance_matrix, method="single", metric=None
        )
        cluster_ids = hierarchy.fcluster(
            linkage_matrix, t=max_distance, criterion="distance"
        )
        return cluster_ids
    return np.array([0])


def aggregate_phospho_groups(df: pd.DataFrame, match_tolerance: int) -> pd.Series:
    """
    This function delocalizes ptm-positions in modified sequences by match_tolerance and combines them if they are present in the data.

    Parameters
    ----------
    df : pd.DataFrame
        a DataFrame with columns <'Modified sequence', 'Delocalized sequence'>
    match_tolerance : int >= 0
        the matching tolerance that specifies how close two ptm sites can be to be considered the same.

    Returns
    -------
    mod_seqs_clusters : pd.Series(<seqs>)
    """
    assert ("Modified sequence" in df) and ("Delocalized sequence" in df)
    assert match_tolerance > 0

    # dont work on input
    df = df.copy()

    # Make a positional array using tqdm progress apply else standard pandas apply
    if hasattr(df, "progress_transform"):
        df["Positional array"] = df["Modified sequence"].progress_apply(
            extract_phos_positions, pattern=PHOSPHORYLATION_PATTERN
        )
    else:
        df["Positional array"] = df["Modified sequence"].apply(
            extract_phos_positions, pattern=PHOSPHORYLATION_PATTERN
        )

    # cluster sequences groups using tqdm progress transform else standard pandas transform
    if hasattr(df, "progress_transform"):
        clusters = df.groupby("Delocalized sequence")[
            "Positional array"
        ].progress_transform(find_clusters, max_distance=int(match_tolerance))
        df["Modified sequence group"] = (
            df["Delocalized sequence"] + "_" + clusters.astype(str)
        )
        mod_seqs_clusters = df.groupby("Modified sequence group")[
            "Modified sequence"
        ].progress_transform(lambda seqs: ";".join(sorted(set(seqs))))
    else:
        clusters = df.groupby("Delocalized sequence")["Positional array"].transform(
            find_clusters, max_distance=int(match_tolerance)
        )
        df["Modified sequence group"] = (
            df["Delocalized sequence"] + "_" + clusters.astype(str)
        )
        mod_seqs_clusters = df.groupby("Modified sequence group")[
            "Modified sequence"
        ].transform(lambda seqs: ";".join(sorted(set(seqs))))

    return mod_seqs_clusters


def delocalize_phospho_sequence(mod_seqs: pd.Series) -> pd.Series:
    """
    Removes the phospho position and adds a _N at the end of the sequence to indicate the number of phosphorylations.
    All other modifications remain untouched.
    Columns-wise operation is 10x faster than apply.

    Parameters
    ----------
    mod_seqs : pd.Series(<sequences>)

    Returns
    -------
    mod_seqs : pd.Series(<seqs>)
    """
    # Count
    ph_count = mod_seqs.str.count("(ph)").replace(np.nan, 0)
    # De-localize
    mod_seqs = (
        mod_seqs.str.replace(r"\(ph\)", "", regex=True)
        + "_"
        + ph_count.astype(int).astype(str)
    )
    return mod_seqs


def make_monophos_versions(
    mod_seq: str, pattern: re.Pattern = PHOSPHORYLATION_PATTERN
) -> list[str]:
    """
    This function returns a list of all mono-phosphorylated peptide versions given the available positions in the input sequence.
    The order of the output is sorted by the modification position.

    Parameters
    ----------
    mod_seqs : str

    Returns
    -------
    out : list(<seqs>, <seqs>, ...)
    """
    base_seq = pattern.sub("", mod_seq)
    out = []
    for pos in extract_phos_positions(mod_seq, pattern):
        out.append(base_seq[: (pos + 1)] + "(ph)" + base_seq[(pos + 1) :])
    return out
