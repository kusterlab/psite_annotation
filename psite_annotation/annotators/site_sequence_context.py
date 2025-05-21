from typing import Dict, List, Optional, Tuple
import collections
import re

import pandas as pd

from .annotator_base import check_columns
from .peptide_position import _read_fasta_maxquant, _read_fasta_phosphositeplus


SITE_POSITION_PATTERN = re.compile(r"[a-zA-Z\-0-9]*_[A-Z][0-9]*")


class SiteSequenceContextAnnotator:
    """Annotate pandas dataframe with +/- 15 amino acids around each of the modified sites, separated by semicolons.

    Example:
        ::
        
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
        retain_other_mods: bool = False,
        return_unique: bool=False,
        return_sorted: bool=False,
        organism: str = "human",
    ):
        """
        Initialize the input files and options for PeptidePositionAnnotator.

        Args:
            annotation_file: fasta file containing protein sequences
            pspInput: set to True if fasta file was obtained from PhosphositePlus
            context_left: number of amino acids to the left of the modification to include
            context_right: number of amino acids to the right of the modification to include
            retain_other_mods: retain other modifications from the modified peptide in the sequence context in lower case

        """
        self.annotation_file = annotation_file
        self.pspInput = pspInput
        self.protein_sequences = None
        self.context_left = context_left
        self.context_right = context_right
        self.retain_other_mods = retain_other_mods
        self.return_unique = return_unique
        self.return_sorted = return_sorted
        self.organism = organism

    def load_annotations(self) -> None:
        """Reads in protein sequences from fasta file."""
        readFasta = _read_fasta_maxquant
        if self.pspInput:
            readFasta = _read_fasta_phosphositeplus

        self.protein_sequences = collections.defaultdict(str)
        for proteinId, seq in readFasta(self.annotation_file, organism=self.organism):
            self.protein_sequences[proteinId] = seq

    @check_columns(["Site positions"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds columns regarding the peptide position within the protein to a pandas dataframe.

        Adds the following annotation columns to dataframe\:
        
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
                retain_other_mods=self.retain_other_mods,
                return_unique=self.return_unique,
                return_sorted=self.return_sorted,
            )
        )
        return annotated_df


def _get_site_sequence_contexts(
        site_position_string: str,
        protein_sequences: Dict[str, str],
        return_unique: bool = False,
        return_sorted: bool = False,
        **kwargs
) -> str:
    if len(site_position_string) == 0:
        return ""

    site_position_strings = site_position_string.split(";")
    contexts = map(
        lambda x: _get_site_sequence_context(
            x, protein_sequences, site_position_strings, **kwargs
        ),
        site_position_strings,
    )

    if return_unique:
        contexts = set(contexts)

    if return_sorted:
        contexts = sorted(contexts)

    return ";".join(contexts)


def _get_site_sequence_context(
    site_position_string: str,
    protein_sequences: Dict[str, str],
    all_site_position_strings: Optional[List[str]] = None,
    context_left: int = 15,
    context_right: int = 15,
    retain_other_mods: bool = False,
) -> str:
    """Get sequence context with +/-15 amino acids around the modification site.

    Args:
        site_position_string: UniProt protein identifier with its modified amino acid and position, e.g. Q86U42_S19
        protein_sequences: dictionary of UniProt protein identifiers to protein sequence
        all_site_position_strings: list of all site position strings for modifications, necessary for retain_other_mods
        context_left: number of amino acids to the left of the modification to include
        context_right: number of amino acids to the right of the modification to include
        retain_other_mods: retain other modifications from the modified peptide in the sequence context in lower case

    Returns:
        str: sequence context with +/-15 amino acids around the modification site
    """
    proteinId, sitePos, mod = _unpack_site_position_string(site_position_string)
    prefix, suffix = "", ""

    proteinLength = len(protein_sequences[proteinId])
    if (
        proteinLength == 0
        or sitePos >= proteinLength
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

    sequence_context_string = (
        prefix
        + protein_sequences[proteinId][contextStart:sitePos]
        + mod
        + protein_sequences[proteinId][sitePos + 1 : contextEnd]
        + suffix
    )

    if retain_other_mods:
        for other_site_position_string in all_site_position_strings:
            sequence_context_string = _add_modification_to_sequence_context(
                sequence_context_string,
                site_position_string,
                other_site_position_string,
                context_left,
            )

    return sequence_context_string


def _add_modification_to_sequence_context(
    sequence_context_string,
    site_position_string,
    other_site_position_string,
    context_left,
):
    if other_site_position_string == site_position_string:
        return sequence_context_string

    proteinId, sitePos, _ = _unpack_site_position_string(site_position_string)
    other_proteinId, other_sitePos, other_mod = _unpack_site_position_string(
        other_site_position_string
    )
    if other_proteinId != proteinId:
        return sequence_context_string

    relative_pos = other_sitePos - sitePos + context_left
    if relative_pos < 0 or relative_pos >= len(sequence_context_string):
        return sequence_context_string

    if other_mod != sequence_context_string[relative_pos].lower():
        raise ValueError(
            f"Incorrect modified amino acid at position {relative_pos} in {sequence_context_string}. "
            f"Expected {other_mod.upper()}, encountered {sequence_context_string[relative_pos]}"
        )

    return (
        sequence_context_string[:relative_pos]
        + other_mod
        + sequence_context_string[relative_pos + 1 :]
    )


def _unpack_site_position_string(site_position_string: str) -> Tuple[str, int, str]:
    if not re.match(SITE_POSITION_PATTERN, site_position_string):
        raise ValueError(
            f"Invalid format for site_position_string: {site_position_string}"
        )
    proteinId, modString = site_position_string.split("_")
    sitePos = int(modString[1:]) - 1
    mod = modString[0].lower()
    return proteinId, sitePos, mod
