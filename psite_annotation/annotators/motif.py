import collections
import csv
import re
from typing import IO, Dict, List, Union

import pandas as pd

from .annotator_base import check_columns
from .separated_strings import combine_separated_strings


class MotifAnnotator:
    """Annotate pandas dataframe with motifs from uniprot.

    Requires "Site sequence context" column in the dataframe to be present.
    The "Site sequence context" column can be generated with PeptidePositionAnnotator().

    Example:
        ::

            annotator = MotifAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: Union[str, IO]):
        """
        Initialize the input files and options for MotifAnnotator.

        Args:
            annotation_file: tab separated file with motifs and their identifiers

        """
        self.annotation_file = annotation_file
        self.motif_dict = None

    def load_annotations(self) -> None:
        """Reads in tab separated file with motif annotations compiled from multiple databases."""
        self.motif_dict = collections.defaultdict(list)
        # Uniprot_ACC, Motif, Start_pos, End_pos

        annotation_file = self.annotation_file
        if isinstance(annotation_file, str):
            annotation_file = open(self.annotation_file)

        self.motif_dict = dict()

        # Database, Accession, Name, Description, Regex
        reader = csv.reader(annotation_file, delimiter="\t")
        header = next(reader)
        identifier_col = header.index("Identifier")
        regex_col = header.index("Regex")
        for row in reader:
            self.motif_dict[row[identifier_col]] = re.compile(
                ".*" + row[regex_col] + ".*"
            )

    @check_columns(["Site sequence context"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds column with motifs the site sequence context matches with.

        Adds the following annotation columns to dataframe\:
        
        - Motifs = semicolon separated list of motifs that match with the site sequence contexts

        Args:
            df: pandas dataframe with "Site sequence context" column
            inplace: Whether to modify the DataFrame rather than creating a new one.

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        annotated_df["Motifs"] = annotated_df["Site sequence context"].apply(
            lambda x: _get_motifs(self.motif_dict, x)
        )
        return annotated_df


def _get_motifs(
    motifDict: Dict[str, List[str]],
    siteContexts: str,
) -> str:
    """Get semicolon separated list of motifs that mathc the site sequence contexts.

    Args:
        motifDict: dictionary mapping UniProt identifiers to a list of motifs with start and end positions
        siteContexts: semicolon separated list of site sequence contexts (+/- 15 amino acids) around the site

    Returns:
        str: semicolon separated list of motifs that match the site sequence contexts

    """
    motifNames = list()
    for sc in siteContexts.split(";"):
        for motifName, motifRegex in motifDict.items():
            if re.match(motifRegex, sc):
                motifNames.append(motifName)

    return combine_separated_strings(motifNames)
