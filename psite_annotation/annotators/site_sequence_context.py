import collections
from typing import Dict, Tuple

import pandas as pd

from .annotator_base import check_columns
from .peptide_position import _read_fasta_maxquant, _read_fasta_phosphositeplus


class SiteSequenceContextAnnotator:
    """Annotate pandas dataframe with +/- 15 amino acids around each of the modified sites, separated by semicolons.

    Typical usage example:
      annotator = SiteSequenceContextAnnotator(<path_to_annotation_file>)
      annotator.load_annotations()
      df = annotator.annotate(df)
    """

    def __init__(
        self,
        annotation_file: str,
        pspInput: bool = False,
        context_left: int = 15,
        context_right: int = 15,
    ):
        """
        Initialize the input files and options for PeptidePositionAnnotator.

        Args:
            annotation_file: fasta file containing protein sequences
            pspInput: set to True if fasta file was obtained from PhosphositePlus
            context_left: number of amino acids to the left of the modification to include
            context_right: number of amino acids to the right of the modification to include

        """
        self.annotation_file = annotation_file
        self.pspInput = pspInput
        self.protein_sequences = None
        self.context_left = context_left
        self.context_right = context_right

    def load_annotations(self) -> None:
        """Reads in protein sequences from fasta file."""
        readFasta = _read_fasta_maxquant
        if self.pspInput:
            readFasta = _read_fasta_phosphositeplus

        self.protein_sequences = collections.defaultdict(str)
        for proteinId, seq in readFasta(self.annotation_file):
            self.protein_sequences[proteinId] = seq

    @check_columns(["Site positions"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds columns regarding the peptide position within the protein to a pandas dataframe.

        Adds the following annotation columns to dataframe:
        - 'Site sequence context' = +/- 15 amino acids around each of the modified sites, separated by semicolons

        Args:
            df: pandas dataframe to be annotated which contains a column "Site positions"
            inplace: add the new column to df in place

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        annotated_df["Site sequence context"] = annotated_df["Site positions"].apply(
            lambda x: _get_site_sequence_contexts(
                x,
                self.protein_sequences,
                context_left=self.context_left,
                context_right=self.context_right,
            )
        )
        return annotated_df


def _get_site_sequence_contexts(
    site_position_string: str, protein_sequences: Dict[str, str], **kwargs
) -> str:
    if len(site_position_string) == 0:
        return ""

    site_position_strings = site_position_string.split(";")
    contexts = map(
        lambda x: _get_site_sequence_context(x, protein_sequences, **kwargs),
        site_position_strings,
    )
    return ";".join(sorted(set(contexts)))


def _get_site_sequence_context(
    site_position_string: str,
    protein_sequences: Dict[str, str],
    context_left: int = 15,
    context_right: int = 15,
) -> str:
    """Get sequence context with +/-15 amino acids around the modification site.

    Args:
        site_position_string: UniProt protein identifier with its modified amino acid and position, e.g. Q86U42_S19
        protein_sequences: dictionary of UniProt protein identifiers to protein sequence
        context_left: number of amino acids to the left of the modification to include
        context_right: number of amino acids to the right of the modification to include

    Returns:
        str: sequence context with +/-15 amino acids around the modification site
    """
    proteinId, sitePos, mod = _unpack_site_position_string(site_position_string)
    prefix, suffix = "", ""

    proteinLength = len(protein_sequences[proteinId])
    if (
        proteinLength == 0
        or sitePos > proteinLength
        or protein_sequences[proteinId][sitePos].lower() != mod
    ):
        return ""

    contextStart = sitePos - context_left
    contextEnd = sitePos + context_right + 1
    if sitePos - context_left < 0:
        contextStart = 0
        prefix = "_" * (context_left - sitePos)

    if sitePos + context_right + 1 > proteinLength:
        contextEnd = proteinLength
        suffix = "_" * (context_right + 1 + sitePos - proteinLength)
    return (
        prefix
        + protein_sequences[proteinId][contextStart:sitePos]
        + mod
        + protein_sequences[proteinId][sitePos + 1 : contextEnd]
        + suffix
    )


def _unpack_site_position_string(site_position_string: str) -> Tuple[str, str, str]:
    proteinId, modString = site_position_string.split("_")
    sitePos = int(modString[1:]) - 1
    mod = modString[0].lower()
    return proteinId, sitePos, mod
