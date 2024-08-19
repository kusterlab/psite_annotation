import collections
import csv
from typing import IO, Dict, List, Union

import pandas as pd

from .annotator_base import check_columns


class DomainAnnotator:
    """Annotate pandas dataframe with domains from uniprot.

    Requires 'Matched proteins', 'Start positions', 'End positions' columns in the dataframe to be annotated.
    The 'Matched proteins', 'Start positions', 'End positions' columns can be generated with PeptidePositionAnnotator().

    Example:
        ::

            annotator = DomainAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(self, annotation_file: Union[str, IO]):
        """
        Initialize the input files and options for DomainAnnotator.

        Args:
            annotation_file: comma separated file with domains and their positions within the protein

        """
        self.annotation_file = annotation_file
        self.domain_dict = None

    def load_annotations(self) -> None:
        """Reads in comma separated file with domain annotations extracted from ProteomicsDB."""
        self.domain_dict = collections.defaultdict(list)
        # Uniprot_ACC, Domain, Start_pos, End_pos

        annotation_file = self.annotation_file
        if isinstance(annotation_file, str):
            annotation_file = open(self.annotation_file)

        reader = csv.reader(annotation_file)
        _ = next(reader)
        for row in reader:
            self.domain_dict[row[0]].append((int(row[2]), int(row[3]), row[1]))

    @check_columns(["Matched proteins", "Start positions", "End positions"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds column with domains the peptide overlaps with.

        Adds the following annotation columns to dataframe\:

        - Domains = semicolon separated list of domains that overlap with the peptide

        Args:
            df: pandas dataframe with 'Proteins', 'Start positions' and 'End positions' columns
            inplace: Whether to modify the DataFrame rather than creating a new one.

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        annotated_df["Domains"] = df[
            ["Matched proteins", "Start positions", "End positions"]
        ].apply(
            lambda x: _get_domains(
                x["Matched proteins"],
                self.domain_dict,
                x["Start positions"],
                x["End positions"],
            ),
            axis=1,
        )

        return annotated_df


def _get_domains(
    proteinIds: str,
    domainDict: Dict[str, List[str]],
    startPositions: str,
    endPositions: str,
) -> str:
    """Get semicolon separated list of domains that overlap the stretch between startPositions and endPositions.

    Args:
        proteinIds: semicolon separated string with UniProt identifiers
        domainDict: dictionary mapping UniProt identifiers to a list of domains with start and end positions
        startPositions: semicolon separated string with starting positions corresponding to the proteinId in the same
            position in the semicolon separated string
        endPositions: semicolon separated string with ending positions corresponding to the proteinId in the same
            position in the semicolon separated string

    Returns:
        str: semicolon separated list of domains that overlap the stretch between startPositions and endPositions.

    """
    if str(proteinIds) == "nan":
        return ""

    domains = list()
    for proteinId, startPos, endPos in zip(  # noqa: B905
        proteinIds.split(";"),
        startPositions.split(";"),
        endPositions.split(";"),
    ):
        startPos, endPos = int(startPos), int(endPos)
        for domainStartPos, domainEndPos, domainName in domainDict.get(proteinId, []):
            if startPos < domainEndPos and endPos > domainStartPos:
                domains.append(domainName)

    return ";".join(sorted(set(domains)))
