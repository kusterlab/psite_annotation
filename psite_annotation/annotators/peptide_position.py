import collections
import itertools
import re

import pandas as pd

from .annotator_base import check_columns

MOD_DICT = {
    "K(ac)": "k",
    "K( GlyGly (K) without TMT)": "k",
    "K( GlyGly )": "k",
    "K(GlyGly)": "k",
    "C(dbia)": "c",
    "K(Acetyl (K))": "k",
    "S(ph)": "s",
    "T(ph)": "t",
    "Y(ph)": "y",
    "S(Phospho (STY))": "s",
    "T(Phospho (STY))": "t",
    "Y(Phospho (STY))": "y",
    "pS": "s",
    "pT": "t",
    "pY": "y",
    "(ac)": "",
    "(Acetyl (Protein N-term))": "",
    "(ox)": "",
    "(Oxidation (M))": "",
    "_": "",
}
MOD_REGEX = re.compile("|".join(map(re.escape, MOD_DICT.keys())))
MOD_PATTERN = re.compile(r"([cksty])")


class PeptidePositionAnnotator:
    """Annotate pandas dataframe with positions of the peptide within the protein sequence based on a fasta file.

    Typical usage example:
      annotator = PeptidePositionAnnotator(<path_to_annotation_file>)
      annotator.load_annotations()
      df = annotator.annotate(df)
    """

    def __init__(
        self,
        annotation_file: str,
        pspInput: bool = False,
        returnAllPotentialSites: bool = False,
        mod_regex=MOD_REGEX,
        mod_pattern=MOD_PATTERN,
    ):
        """
        Initialize the input files and options for PeptidePositionAnnotator.

        Args:
            annotation_file: fasta file containing protein sequences
            pspInput: set to True if fasta file was obtained from PhosphositePlus
            returnAllPotentialSites: set to True if all S, T and Y positions should be returned as potention p-sites.

        """
        self.annotation_file = annotation_file
        self.pspInput = pspInput
        self.returnAllPotentialSites = returnAllPotentialSites
        self.protein_sequences = None
        self.MOD_REGEX = mod_regex
        self.MOD_PATTERN = mod_pattern

    def load_annotations(self) -> None:
        """Reads in protein sequences from fasta file."""
        readFasta = _read_fasta_maxquant
        if self.pspInput:
            readFasta = _read_fasta_phosphositeplus

        self.protein_sequences = collections.defaultdict(str)
        for proteinId, seq in readFasta(self.annotation_file):
            self.protein_sequences[proteinId] = seq

    @check_columns(["Proteins", "Modified sequence"])
    def annotate(self, df: pd.DataFrame, inplace: bool = False) -> pd.DataFrame:
        """Adds columns regarding the peptide position within the protein to a pandas dataframe.

        Adds the following annotation columns to dataframe:
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

        Returns:
            pd.DataFrame: annotated dataframe

        """
        annotated_df = df
        if not inplace:
            annotated_df = df.copy()

        annotated_df[
            [
                "Matched proteins",
                "Start positions",
                "End positions",
                "Site positions",
            ]
        ] = annotated_df[["Proteins", "Modified sequence"]].apply(
            lambda x: _get_peptide_positions(
                x["Proteins"],
                self.protein_sequences,
                x["Modified sequence"],
                self.returnAllPotentialSites,
                self.MOD_REGEX,
                self.MOD_PATTERN,
            ),
            axis=1,
            result_type="expand",
        )
        return annotated_df


def _make_maxquant_mods_consistent(text, mod_regex=MOD_REGEX):
    return mod_regex.sub(lambda match: MOD_DICT[match.group(0)], text)


def _get_peptide_positions(
    proteinIds,
    protein_sequences,
    mod_peptide_sequence,
    returnAllPotentialSites=False,
    mod_regex=MOD_REGEX,
    mod_pattern=MOD_PATTERN,
):
    if str(proteinIds) == "nan" or len(mod_peptide_sequence) == 0:
        return ("", "", "", "")

    mod_peptide_sequence = _make_maxquant_mods_consistent(
        mod_peptide_sequence, mod_regex
    )

    if returnAllPotentialSites:
        mod_peptide_sequence = (
            mod_peptide_sequence.replace("S", "s")
            .replace("T", "t")
            .replace("Y", "y")
            .replace("K", "k")
            .replace("c", "C")
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

    return (
        ";".join(map(str, matchedProteins)),
        ";".join(map(str, startPositions)),
        ";".join(map(str, endPositions)),
        ";".join(sorted(set(proteinPositions))),
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


def _read_fasta_maxquant(file_path, parse_id=parse_uniprot_id):
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


def _read_fasta_phosphositeplus(file_path, parse_id=lambda x: x.split("|")[3]):
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
                    if line[1:].split("|")[2] != "human":
                        name = False
            else:
                seq.append(line)
