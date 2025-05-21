import pandas as pd

from .annotator_base import check_columns
from .separated_strings import combine_separated_strings, merge_on_separated_string


class PSPRegulatoryAnnotator:
    """Annotate pandas dataframe with regulatory functions according to PhosphositePlus.

    Example:
        ::

            annotator = PSPRegulatoryAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: str, organism: str = "human"):
        """
        Initialize the input files and options for PSPRegulatoryAnnotator.

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

        # format: O75822_S11
        self.psp_df["Site positions"] = (
            self.psp_df["ACC_ID"]
            + "_"
            + self.psp_df["MOD_RSD"].apply(lambda x: x.split("-")[0])
        )

        regulatory_columns = [
            "ON_FUNCTION",
            "ON_PROCESS",
            "ON_PROT_INTERACT",
            "ON_OTHER_INTERACT",
            "NOTES",
        ]
        self.psp_df = (
            self.psp_df[["Site positions"] + regulatory_columns].fillna("").astype(str)
        )
        self.psp_df = self.psp_df.groupby("Site positions", sort=False)[
            regulatory_columns
        ].agg({k: "; ".join for k in regulatory_columns})

        self.psp_df = self.psp_df.rename(columns=lambda x: f"PSP_{x}")

        # make the "Site positions" an ordinary column again
        self.psp_df = self.psp_df.reset_index()

    @check_columns(["Site positions"])
    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds columns with number of studies.

        Adds the following annotation columns to dataframe\:

        - PSP_ON_FUNCTION = functional annotations for downstream regulation
        - PSP_ON_PROCESS = process annotations for downstream regulation
        - PSP_ON_PROT_INTERACT = protein interactions
        - PSP_ON_OTHER_INTERACT = other interactions
        - PSP_NOTES = regulatory site notes

        Args:
            df: pandas dataframe with "Site positions" column

        Returns:
            pd.DataFrame: annotated dataframe

        """
        agg_func = {
            k: lambda x: combine_separated_strings(x, sep="; ")
            for k in self.psp_df.columns
        }
        agg_func["Site positions"] = ";".join

        annotated_df = merge_on_separated_string(
            df, self.psp_df, "Site positions", agg_func=agg_func
        )
        return annotated_df
