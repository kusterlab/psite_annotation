import pandas as pd

from .annotator_base import check_columns
from .separated_strings import combine_separated_strings, merge_on_separated_string


class InVitroKinasesAnnotator:
    """Annotate pandas dataframe with upstream in vitro kinases according to Sugiyama et al (2019).

    https://www.nature.com/articles/s41598-019-46385-4

    Example:
        ::
            
            annotator = InVitroKinasesAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: str):
        """
        Initialize the input files and options for InVitroKinasesAnnotator.

        Args:
            annotation_file: tab separated file with in vitro kinase annotations

        """
        self.annotation_file = annotation_file
        self.kinase_df = None

    def load_annotations(self) -> None:
        """Reads in tab separated file with in vitro annotations.

        Columns: Type	Kinase	Uniprot ID	Protein description	Position	SIDIC	PTMscore UniprotId UniprotId_Site
        """
        self.kinase_df = pd.read_csv(self.annotation_file, sep="\t")

        # format: O75822_S11
        self.kinase_df = self.kinase_df[["Kinase", "UniprotId_Site"]]

        self.kinase_df = self.kinase_df.rename(
            columns={"Kinase": "In Vitro Kinases", "UniprotId_Site": "Site positions"}
        )

    @check_columns(["Site positions"])
    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds column with phosphorylating kinases.

        Adds the following annotation columns to dataframe\:
        
        - In Vitro Kinases = all phosphorylating kinases according to the Sugiyama in vitro kinase-substrate study

        Args:
            df: pandas dataframe with "Site positions" column

        Returns:
            pd.DataFrame: annotated dataframe

        """
        agg_func = {k: combine_separated_strings for k in self.kinase_df.columns}
        agg_func["Site positions"] = ";".join

        annotated_df = merge_on_separated_string(
            df, self.kinase_df, "Site positions", agg_func=agg_func
        )
        return annotated_df
