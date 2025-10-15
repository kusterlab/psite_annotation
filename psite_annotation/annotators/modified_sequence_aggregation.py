# adapted from phospho_delocalization.py (Florian P. Bayer - 2025)
from typing import Any

import pandas as pd
import numpy as np

from .annotator_base import check_columns


class ModifiedSequenceAggregatorAnnotator:
    """Annotate and aggregate pandas dataframe with representative modified sequence from a modified sequence group.

    Example:
        ::

            annotator = ReprentativeModifiedSequenceAnnotator()
            df = annotator.annotate(df)
    """

    def __init__(
        self,
        experiment_cols: list[str],
        agg_func: str = "mean",
        agg_cols: dict[str, Any] = None,
    ) -> None:
        """
        Initialize the options for ReprentativeModifiedSequenceAnnotator.

        Args:
            experiment_cols: list of column names with quantitative values.
            agg_func: function to aggregate quantitative values within each group, e.g. 'mean', 'sum', etc.

        """
        self.experiment_cols = experiment_cols
        self.agg_func = agg_func
        self.agg_cols = {}
        if agg_cols:
            self.agg_cols = agg_cols

    def load_annotations(self) -> None:
        pass

    @check_columns(["Modified sequence", "Delocalized sequence", "Modified sequence group"])
    def annotate(self, df: pd.DataFrame, ) -> pd.DataFrame:
        """
        Group delocalized phospho-forms and aggregate their quantitative values.

        This function identifies peptide sequences that differ only by the position
        of their phosphorylation (`(ph)`) group and collapses them into
        "delocalized" groups. Each group contains all modified sequence variants
        that represent the same underlying peptide backbone.

        The following columns are added to the dataframe:

        - 'Modified sequence representative' = A single representative sequence
        selected from the group, i.e. the most frequently measured across experiments.
        - 'Modified sequence representative degree' = Fraction of summed observation
        frequency contributed by the representative peptide.

        All experiment columns (e.g. `"Experiment 1"`, `"Experiment 2"`, …) are aggregated
        per group by summing the intensities of member sequences.

        Args:
            df: Input dataframe with:
                - `"Modified sequence"` column containing peptide strings with `(ph)` annotations
                - 'Delocalized sequence' = Canonical unmodified backbone with an index
                suffix to distinguish the number of modifications.
                - 'Modified sequence group' = All peptide variants belonging to the same
                delocalized group, concatenated with semicolons.

        Returns:
            pd.DataFrame: Dataframe with grouped phospho-forms and aggregated intensities.
        """
        # TODO: implement inplace option. does not work currently because groupby().agg() cannot be done inplace
        annotated_df = df

        # Determine representative sequence for each cluster based on the observations in the experiments.
        # The Modified sequence representative degree gives the proportion of the representative relative to all observations.
        df_representative = find_representative_modified_sequence(
            annotated_df, self.experiment_cols
        )

        # Aggregate experiments per group.
        df_agg = (
            annotated_df.groupby(["Delocalized sequence", "Modified sequence group"])
            .agg({x: self.agg_func for x in self.experiment_cols} | self.agg_cols)
            .reset_index()
        )

        df_agg = df_representative.merge(df_agg, on="Modified sequence group")

        return df_agg


def find_representative_modified_sequence(
    df: pd.DataFrame, observation_cols: list[str]
) -> pd.DataFrame:
    """
    This function counts the number of observations of a specific modified sequence and defines the most representative sequence as
    the one with most observations. Missing is indicated as NaN. Any other value is considered an observation.

    Parameters
    ----------
    df : pd.DataFrame
        a DataFrame with columns <'Modified sequence', 'Delocalized sequence'>, observation_cols
    observation_cols : list of cols
        the names of the columns that are used for counting if a peptide was observed.

    Returns
    -------
    df : pd.DataFrame
        A new DataFrame with cols ['Modified sequence representative', 'Modified sequence representative degree']
    """
    col = "Modified sequence"
    group_col = "Modified sequence group"
    count_col = "Modified sequence count"
    assert (
        (col in df)
        and (group_col in df)
        and all(c in df.columns for c in observation_cols)
    )

    # Copy so the original df is not modified
    df = df[[col, group_col] + list(observation_cols)].copy()

    # Do the grouping and counting
    df[count_col] = df[observation_cols].isna().apply(np.logical_not).sum(axis=1)
    reprentitive_idx = (
        df.groupby(group_col)[count_col].transform(lambda x: x.idxmax()).values
    )
    df[f"{col} representative"] = df.loc[reprentitive_idx, col].values
    df[f"{col} representative degree"] = df.groupby(group_col)[count_col].transform(
        lambda x: max(x) / sum(x)
    )
    out = df.groupby(group_col)[
        [f"{col} representative", f"{col} representative degree"]
    ].first()
    return out
