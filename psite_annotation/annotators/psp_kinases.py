import pandas as pd

from .annotator_base import check_columns
from .separated_strings import combine_separated_strings, merge_on_separated_string


class PSPKinasesAnnotator:
    """Annotate pandas dataframe with upstream kinases according to PhosphositePlus.

    Example:
        ::

            annotator = PSPKinasesAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(
        self,
        annotation_file: str,
        output_gene_names: bool = False,
        organism: str = "human",
    ):
        """
        Initialize the input files and options for PSPKinasesAnnotator.

        Args:
            annotation_file: tab separated file with PhosphositePlus kinase annotations
            output_gene_names: set to True to output the gene names instead of the kinase names

        """
        self.annotation_file = annotation_file
        self.output_gene_names = output_gene_names
        self.psp_df = None
        self.organism = organism

    def load_annotations(self) -> None:
        """Reads in tab separated file with PhosphositePlus annotations."""
        self.psp_df = pd.read_csv(
            self.annotation_file, sep="\t", skiprows=3, encoding="utf-8"
        )
        self.psp_df = self.psp_df[
            (self.psp_df["KIN_ORGANISM"] == self.organism)
            & (self.psp_df["SUB_ORGANISM"] == self.organism)
        ]

        # format: O75822_S11
        # in contrast to the other two PSP annotation files, the MOD_RSD does not contain the PTM type, because all
        # sites in this file are phosphorylations
        self.psp_df["Site positions"] = (
            self.psp_df["SUB_ACC_ID"] + "_" + self.psp_df["SUB_MOD_RSD"]
        )

        col = "GENE" if self.output_gene_names else "KINASE"
        self.psp_df = self.psp_df.groupby(["Site positions"], sort=False)[[col]].agg(
            ";".join
        )

        self.psp_df = self.psp_df.rename(columns={col: "PSP Kinases"})

        # make the "Site positions" an ordinary column again
        self.psp_df = self.psp_df.reset_index()

    @check_columns(["Site positions"])
    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds column with phosphorylating kinases.

        Adds the following annotation columns to dataframe\:

        - PSP Kinases = all phosphorylating kinases according to PhosphoSitePlus, no distinction is made between
          in vivo and in vitro evidence (this can be added in the future, if necessary)

        Args:
            df: pandas dataframe with "Site positions" column

        Returns:
            pd.DataFrame: annotated dataframe

        """
        agg_func = {k: combine_separated_strings for k in self.psp_df.columns}
        agg_func["Site positions"] = ";".join

        annotated_df = merge_on_separated_string(
            df, self.psp_df, "Site positions", agg_func=agg_func
        )
        return annotated_df
