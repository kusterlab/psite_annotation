import pandas as pd

from .annotator_base import check_columns


class PTMTurnoverAnnotator:
    """Annotate pandas dataframe with Jana's PTM Turnover data.

    https://www.nature.com/articles/s41467-021-27639-0

    Example:
        ::

            annotator = PTMTurnoverAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: str):
        """
        Initialize the input files and options for PTMTurnoverAnnotator.

        Args:
            annotation_file: comma separated file with PTM turnover annotations

        """
        self.annotation_file = annotation_file
        self.turnover_df = None

    def load_annotations(self) -> None:
        """Reads in comma separated file with PTM turnover annotations.

        The "PTM_Turnover" column can contain one of the following values:
        - slower
        - faster
        - FALSE
        - ""
        FALSE indicates there as "no difference" in PTM turnover
        A missing value indicates there were "not enough replicates"

        """
        self.turnover_df = pd.read_csv(self.annotation_file, sep=",")

        self.turnover_df = self.turnover_df[["Significant global", "Peptidoforms"]]
        self.turnover_df.columns = ["PTM_Turnover", "Peptidoforms"]

        self.turnover_df["PTM_Turnover"] = (
            self.turnover_df["PTM_Turnover"]
            .replace("FALSE", "no difference")
            .fillna("not enough replicates")
        )

        self.turnover_df["Peptidoforms"] = self.turnover_df["Peptidoforms"].apply(
            lambda x: str(x).split(";")
        )
        self.turnover_df = self.turnover_df.explode("Peptidoforms")
        self.turnover_df = self.turnover_df.set_index("Peptidoforms")

    @check_columns(["Modified sequence"])
    def annotate(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds column regarding the PTM turnover behavior.

        Adds the following annotation columns to dataframe\:
        
        - 'PTM Turnover' = rate of turnover for the modification sites according to Jana's PTM Turnover data

        Args:
            df: pandas dataframe with 'Modified sequence' column

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df.set_index("Modified sequence")
        annotated_df = annotated_df.join(self.turnover_df)
        annotated_df = annotated_df.reset_index().rename(
            columns={"index": "Modified sequence"}
        )
        return annotated_df
