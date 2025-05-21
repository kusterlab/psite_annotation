import collections
import itertools
import re
from typing import Dict, Pattern, Tuple

import pandas as pd

from .annotator_base import check_columns

MOD_DICT = {
    "S(ph)": "s",
    "T(ph)": "t",
    "Y(ph)": "y",
    "S(Phospho (STY))": "s",
    "T(Phospho (STY))": "t",
    "Y(Phospho (STY))": "y",
    "pS": "s",
    "pT": "t",
    "pY": "y",
}


class PeptidePositionAnnotator:
    """Annotate pandas dataframe with positions of the peptide within the protein sequence based on a fasta file.

    Example:
        ::

            annotator = PeptidePositionAnnotator(<path_to_annotation_file>)
            annotator.load_annotations()
            df = annotator.annotate(df)
    """

    def __init__(
        self,
        annotation_file: str,
        pspInput: bool = False,
        returnAllPotentialSites: bool = False,
        localization_uncertainty: int = 0,
        mod_dict: Dict[str, str] = MOD_DICT,
        return_unique: bool = False,
        return_sorted: bool = False,
        organism: str = "human",
    ) -> None:
        """
        Initialize the input files and options for PeptidePositionAnnotator.

        Args:
            annotation_file: fasta file containing protein sequences
            pspInput: set to True if fasta file was obtained from PhosphositePlus
            returnAllPotentialSites: return all modifiable positions within the peptide as potential p-sites.
            localization_uncertainty: return all modifiable positions within n positions of modified sites as potential p-sites.
            mod_regex: regex to capture all modification strings

        """
        self.annotation_file = annotation_file
        self.pspInput = pspInput
        self.returnAllPotentialSites = returnAllPotentialSites
        self.localization_uncertainty = localization_uncertainty
        self.protein_sequences = None
        self.mod_dict = mod_dict
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

    @check_columns(["Proteins", "Modified sequence"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds columns regarding the peptide position within the protein to a pandas dataframe.

        Adds the following annotation columns to dataframe\:

        - 'Matched proteins' = subset of 'Proteins' in the input column in which the protein could indeed be found. If
          the same peptide is found multiple times, the protein identifier will be repeated.
        - 'Start positions' = starting positions of the modified peptide in the protein sequence (1-based, methionine is
          counted). If multiple isoforms/proteins contain the sequence, the starting positions are separated by
          semicolons in the same order as they are listed in the 'Matched proteins' column
        - 'End positions' = end positions of the modified peptide in the protein sequence (see above for details)
        - 'Site positions' = position of the modification (see 'Start positions' above for details on how the position
          is counted)

        Args:
            df: pandas dataframe to be annotated with "Proteins" and "Modified sequence" columns
            inplace: add the new column to df in place

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        mod_regex = _get_mod_regex(self.mod_dict)
        mod_pattern = _get_mod_pattern(self.mod_dict)
        potential_mods = _get_single_letter_mods(self.mod_dict)
        annotated_df[
            [
                "Matched proteins",
                "Start positions",
                "End positions",
                "Site positions",
            ]
        ] = annotated_df[["Proteins", "Modified sequence"]].apply(
            lambda x: _get_peptide_positions(
                proteinIds=x["Proteins"],
                protein_sequences=self.protein_sequences,
                mod_peptide_sequence=x["Modified sequence"],
                return_unique=self.return_unique,
                return_sorted=self.return_sorted,
                returnAllPotentialSites=self.returnAllPotentialSites,
                localization_uncertainty=self.localization_uncertainty,
                mod_dict=self.mod_dict,
                mod_regex=mod_regex,
                mod_pattern=mod_pattern,
                potential_mods=potential_mods,
            ),
            axis=1,
            result_type="expand",
        )
        return annotated_df


def _get_single_letter_mods(mod_dict: Dict[str, str]) -> str:
    """Get a string of all single letter modifications

    Args:
        mod_regex (Dict[str, str]): _description_
    """
    return "".join([x for x in mod_dict.values() if len(x) == 1])


def _get_mod_pattern(mod_dict: Dict[str, str]) -> Pattern:
    """Create regex to capture all single letter modifications

    Args:
        mod_regex (Dict[str, str]): _description_
    """
    single_letter_mods = _get_single_letter_mods(mod_dict)
    return re.compile(r"([" + single_letter_mods + "])")


def _get_mod_regex(mod_dict: Dict[str, str]) -> Pattern:
    return re.compile("|".join(map(re.escape, mod_dict.keys())))


def _remove_modifications(mod_peptide_sequence: str) -> str:
    # Define the regex pattern to match any content between the outermost parentheses
    pattern = re.compile(r"\(([^()]|\([^()]*\))*\)")

    while re.search(pattern, mod_peptide_sequence):
        # Use the sub function to replace the matched pattern with an empty string iteratively
        mod_peptide_sequence = re.sub(pattern, "", mod_peptide_sequence)

    return mod_peptide_sequence.replace("_", "")


def _make_maxquant_mods_consistent(
    text: str, mod_dict: Dict[str, str], mod_regex: Pattern
) -> str:
    return mod_regex.sub(lambda match: mod_dict[match.group(0)], text)


def _get_peptide_positions(
    proteinIds: str,
    protein_sequences: Dict[str, str],
    mod_peptide_sequence: str,
    returnAllPotentialSites: bool = False,
    return_unique: bool = False,
    return_sorted: bool = False,
    localization_uncertainty: int = 0,
    mod_dict: Dict[str, str] = MOD_DICT,
    mod_regex: Pattern = _get_mod_regex(MOD_DICT),
    mod_pattern: Pattern = _get_mod_pattern(MOD_DICT),
    potential_mods: str = _get_single_letter_mods(MOD_DICT),
) -> Tuple[str, str, str, str]:
    if str(proteinIds) == "nan" or len(mod_peptide_sequence) == 0:
        return ("", "", "", "")

    mod_peptide_sequence = _make_maxquant_mods_consistent(
        mod_peptide_sequence, mod_dict, mod_regex
    )
    mod_peptide_sequence = _remove_modifications(mod_peptide_sequence)

    if returnAllPotentialSites:
        localization_uncertainty = len(mod_peptide_sequence)

    if localization_uncertainty > 0:
        mod_peptide_sequence = _apply_localization_uncertainty(
            mod_peptide_sequence, localization_uncertainty, potential_mods
        )

    matchedProteins, startPositions, endPositions, proteinPositions = (
        list(),
        list(),
        list(),
        list(),
    )
    for proteinId in proteinIds.split(";"):
        if len(proteinId) == 0:
            continue

        if "|" in proteinId:
            proteinId = proteinId.split("|")[1]

        start_position = -1
        while True:
            start_position = protein_sequences[proteinId].find(
                mod_peptide_sequence.upper(), start_position + 1
            )
            if start_position == -1:
                break

            end_position = start_position + len(mod_peptide_sequence)

            matchedProteins.append(proteinId)
            startPositions.append(start_position)
            endPositions.append(end_position)

            for mod in mod_pattern.finditer(mod_peptide_sequence):
                sitePos = start_position + mod.start()
                site_position_string = (
                    proteinId + "_" + mod.group(0).upper() + str(sitePos + 1)
                )

                proteinPositions.append(site_position_string)

    if return_unique:
        proteinPositions = set(proteinPositions)

    if return_sorted:
        proteinPositions = sorted(proteinPositions)

    return (
        ";".join(map(str, matchedProteins)),
        ";".join(map(str, startPositions)),
        ";".join(map(str, endPositions)),
        ";".join(map(str, proteinPositions)),
    )


def parse_until_first_space(fasta_id: str) -> str:
    """Parses fasta identifier by taking everything before the first space.

    Args:
        fasta_id (str): protein identifier in fasta file with ">" stripped off

    Returns:
        str: protein identifier up to first space
    """
    return fasta_id.split(" ")[0]


def parse_uniprot_id(fasta_id: str) -> str:
    """Parses uniprot identifier from fasta file identifier.

    Args:
        fasta_id (str): protein identifier in fasta file with ">" stripped off.

    Returns:
        str: UniProt protein identifier.
    """
    protein_id = parse_until_first_space(fasta_id)
    if "|" in protein_id:
        return protein_id.split("|")[1]
    else:
        return protein_id


def _read_fasta_maxquant(file_path, parse_id=parse_uniprot_id, organism: str = "human"):
    name, seq = None, []
    with open(file_path) as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield (name, "".join(seq))

                if len(line) > 1:
                    name, seq = parse_id(line[1:]), []
            else:
                seq.append(line)


def _read_fasta_phosphositeplus(
    file_path, parse_id=lambda x: x.split("|")[3], organism: str = "human"
):
    name, seq = None, []
    with open(file_path, encoding="latin-1") as fp:
        next(fp)
        next(fp)
        next(fp)
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield (name, "".join(seq))

                if len(line) > 1:
                    name, seq = parse_id(line[1:]), []
                    if line[1:].split("|")[2] != organism:
                        name = False
            else:
                seq.append(line)


def _apply_localization_uncertainty(
    mod_peptide_sequence: str, localization_uncertainty: int, potential_mods: str
) -> str:
    mod_positions = [idx for idx, aa in enumerate(mod_peptide_sequence) if aa.islower()]
    potential_mod_positions = [
        idx + x
        for idx in mod_positions
        for x in range(localization_uncertainty * -1, localization_uncertainty + 1)
    ]
    return "".join(
        [
            (
                aa.lower()
                if aa in potential_mods.upper() and x in potential_mod_positions
                else aa
            )
            for x, aa in enumerate(mod_peptide_sequence)
        ]
    )
