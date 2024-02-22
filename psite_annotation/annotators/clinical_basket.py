import pandas as pd

from .annotator_base import check_columns
from .separated_strings import (
    combine_separated_strings,
    explode_on_separated_string,
    keep_unique_and_join,
    merge_on_separated_string,
)


class ClinicalBasketAnnotator:
    """Annotate pandas dataframe with clinical baskets from Annika.

    Requires `Gene names` column in the dataframe to be annotated.

    Example:
        ::

            annotator = ClinicalBasketAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: str):
        """
        Initialize the input files and options for ClinicalBasketAnnotator.

        Args:
            annotation_file: excel file with basket-gene annotations

        """
        self.annotation_file = annotation_file
        self.basket_df = None

    def load_annotations(self) -> None:
        """Reads in excel file with basket-gene annotations.

        Creates a dataframe `basket_df` with two columns\:

        - `Gene names` contains a single gene name
        - `Clinical baskets` contains a semicolon-separated list of basket identifiers the gene is featured in

        Returns:
            None

        """
        self.basket_df = pd.read_excel(self.annotation_file)
        self.basket_df = self.basket_df.rename(
            columns={"GENE NAME": "Gene names", "BASKET": "Clinical baskets"}
        )

        # for phosphosites, the `Gene names` column can contain semicolon-separated strings
        self.basket_df = explode_on_separated_string(self.basket_df, "Gene names")

        self.basket_df = self.basket_df.groupby("Gene names", sort=False)[
            "Clinical baskets"
        ].agg(keep_unique_and_join)

        # make the "Gene names" column an ordinary column again
        self.basket_df = self.basket_df.reset_index()

    @check_columns(["Gene names"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds column with baskets the gene names correspond to.

        Adds the following annotation columns to dataframe\:
        
        - `Clinical baskets` = semicolon separated list of clinical baskets the gene name corresponds to

        Args:
            df: pandas dataframe with 'Gene names' column
            inplace: Whether to modify the DataFrame rather than creating a new one.

        Returns:
            pd.DataFrame: annotated dataframe

        """
        agg_func = {k: combine_separated_strings for k in self.basket_df.columns}
        agg_func["Gene names"] = ";".join

        annotated_df = merge_on_separated_string(
            df, self.basket_df, "Gene names", agg_func=agg_func
        )

        if inplace:
            df["Clinical baskets"] = annotated_df["Clinical baskets"]
        else:
            return annotated_df
