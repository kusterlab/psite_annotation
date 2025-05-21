import pandas as pd

from .annotator_base import check_columns
from .separated_strings import merge_on_separated_string


class PSPStudiesAnnotator:
    """Annotate pandas dataframe with number of high and low-throughput studies according to PhosphositePlus.

    Example:
        ::

            annotator = PSPStudiesAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: str, organism: str = "human"):
        """
        Initialize the input files and options for PSPStudiesAnnotator.

        Args:
            annotation_file: tab separated file with PhosphositePlus annotations

        """
        self.annotation_file = annotation_file
        self.psp_df = None
        self.organism = organism

    def load_annotations(self) -> None:
        """Reads in tab separated file with PhosphositePlus annotations and stores it as a dictionary."""
        self.psp_df = pd.read_csv(
            self.annotation_file, sep="\t", skiprows=3, encoding="utf-8"
        )

        self.psp_df = self.psp_df[self.psp_df["ORGANISM"] == self.organism]

        self.psp_df[["LT_LIT", "MS_LIT", "MS_CST"]] = (
            self.psp_df[["LT_LIT", "MS_LIT", "MS_CST"]].fillna(0).astype(int)
        )

        # format: O75822_S11
        self.psp_df["Site positions"] = (
            self.psp_df["ACC_ID"]
            + "_"
            + self.psp_df["MOD_RSD"].apply(lambda x: x.split("-")[0])
        )

        self.psp_df = self.psp_df[["Site positions", "LT_LIT", "MS_LIT", "MS_CST"]]

        # in case there are multiple entries with the same "Site positions" identifier, take the maximum of studies
        self.psp_df = self.psp_df.groupby("Site positions", sort=False)[
            ["LT_LIT", "MS_LIT", "MS_CST"]
        ].agg({"LT_LIT": "max", "MS_LIT": "max", "MS_CST": "max"})

        self.psp_df = self.psp_df.rename(columns=lambda x: f"PSP_{x}")

        # make the "Site positions" an ordinary column again
        self.psp_df = self.psp_df.reset_index()

    @check_columns(["Site positions"])
    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds columns with number of studies.

        Adds the following annotation columns to dataframe\:

        - LT_LIT = number of low-throughput studies
        - MS_LIT = number of high-throughput Mass Spec studies
        - MS_CST = number of high-throughput Mass Spec studies by CellSignalingTechnologies

        Args:
            df: pandas dataframe with "Site positions" column

        Returns:
            pd.DataFrame: annotated dataframe

        """
        agg_func = {k: ";".join for k in self.psp_df.columns}
        annotated_df = merge_on_separated_string(
            df, self.psp_df, "Site positions", agg_func=agg_func
        )

        return annotated_df
