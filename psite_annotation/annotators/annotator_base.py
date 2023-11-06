from typing import Protocol

import pandas as pd


class Annotator(Protocol):
    """Protocol for annotating a pandas dataframe."""

    def load_annotations(self) -> None:
        """Load annotations, e.g. from a csv file."""

    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add annotation columns to dataframe.

        Args:
            df: pandas dataframe to be annotated

        """
