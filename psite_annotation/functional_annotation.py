import logging
import sys
from typing import Dict

import pandas as pd

from . import annotators
from .config import _getConfigDicts, _getConfigSetting

logger = logging.getLogger(__name__)

__all__ = [
    "domainMappingFile",
    "inVitroKinaseSubstrateMappingFile",
    "motifsFile",
    "turnoverFile",
    "pspAnnotationFile",
    "pspRegulatoryFile",
    "pspKinaseSubstrateFile",
    "pspFastaFile",
    "kinaseLibraryMotifsFile",
    "kinaseLibraryQuantilesFile",
    "addPeptideAndPsitePositions",
    "addSiteSequenceContext",
    "addPSPAnnotations",
    "addPSPRegulatoryAnnotations",
    "addPSPKinaseSubstrateAnnotations",
    "addDomains",
    "addMotifs",
    "addInVitroKinases",
    "addTurnoverRates",
    "addKinaseLibraryAnnotations",
]

defaults, user = _getConfigDicts()

domainMappingFile = _getConfigSetting("domainMappingFile", user, defaults)
inVitroKinaseSubstrateMappingFile = _getConfigSetting(
    "inVitroKinaseSubstrateMappingFile", user, defaults
)
motifsFile = _getConfigSetting("motifsFile", user, defaults)
turnoverFile = _getConfigSetting("turnoverFile", user, defaults)
pspAnnotationFile = _getConfigSetting("pspAnnotationFile", user, defaults)
pspRegulatoryFile = _getConfigSetting("pspRegulatoryFile", user, defaults)
pspKinaseSubstrateFile = _getConfigSetting("pspKinaseSubstrateFile", user, defaults)
pspFastaFile = _getConfigSetting("pspFastaFile", user, defaults)
kinaseLibraryMotifsFile = _getConfigSetting("kinaseLibraryMotifsFile", user, defaults)
kinaseLibraryQuantilesFile = _getConfigSetting(
    "kinaseLibraryQuantilesFile", user, defaults
)


def addPeptideAndPsitePositions(
    df: pd.DataFrame,
    fastaFile: str,
    pspInput: bool = False,
    returnAllPotentialSites: bool = False,
    localization_uncertainty: int = 0,
    context_left: int = 15,
    context_right: int = 15,
    retain_other_mods: bool = False,
    mod_dict: Dict[str, str] = None,
    return_unique: bool = False,
    return_sorted: bool = False,
    organism: str = "human",
) -> pd.DataFrame:
    """Annotate pandas dataframe with positions of the peptide within the protein sequence based on a fasta file.

    Adds the following annotation columns to dataframe\:

    - 'Matched proteins' = subset of 'Proteins' in the input column in which the protein could indeed be found. If the same peptide is found multiple times in the same protein sequence, the protein identifier will be repeated.
    - 'Start positions' = starting positions of the modified peptide in the protein sequence (1-based, methionine is counted). If multiple isoforms/proteins contain the sequence, the starting positions are separated by semicolons in the same order as they are listed in the 'Matched proteins' column
    - 'End positions' = end positions of the modified peptide in the protein sequence (see above for details)
    - 'Site positions' = position of the modification (see 'Start positions' above for details on how the position is counted)
    - 'Site sequence context' = +/- 15 amino acids around each of the modified sites, separated by semicolons

    Example:
        Annotate with psite positions as given by PhosphoSitePlus::

            import psite_annotation as pa
            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)

        Annotate a custom modification::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True, mod_dict={'R[0.9840]': 'r'})


    Required columns:
        :code:`Proteins`, :code:`Modified sequence`

    Args:
        df: pandas dataframe with "Proteins" and "Modified sequence" columns
        fastaFile: fasta file containing protein sequences
        pspInput: set to True if fasta file was obtained from PhosphositePlus
        returnAllPotentialSites: return all modifiable positions within the peptide as potential p-sites.
        localization_uncertainty: return all modifiable positions within n positions of modified sites as potential p-sites.
        context_left: number of amino acids to the left of the modification to include
        context_right: number of amino acids to the right of the modification to include
        retain_other_mods: retain other modifications from the modified peptide in the sequence context in lower case
        mod_dict: dictionary of modifications to single amino acid replacements, e.g. :code:`{"S(ph)": "s", "T(ph)": "t", "Y(ph)": "y"}`. If set to :code:`None`, uses the default annotations for S, T and Y phosphorylation.
        return_unique: eliminate duplicates from the 'Site sequence context' and Site positions' columns, not preserving the order between the them and the rest of the data frame
        return_sorted: sort the 'Site sequence context' and Site positions' columns alphabetically, not preserving the order between the them and the rest of the data frame

    Returns:
        pd.DataFrame: annotated dataframe

    """
    if mod_dict is None:
        mod_dict = annotators.peptide_position.MOD_DICT

    peptide_position_annotator = annotators.PeptidePositionAnnotator(
        fastaFile,
        pspInput=pspInput,
        returnAllPotentialSites=returnAllPotentialSites,
        localization_uncertainty=localization_uncertainty,
        mod_dict=mod_dict,
        return_unique=return_unique,
        return_sorted=return_sorted,
        organism=organism,
    )
    peptide_position_annotator.load_annotations()
    df = peptide_position_annotator.annotate(df)

    site_seq_context_annotator = annotators.SiteSequenceContextAnnotator(
        fastaFile,
        pspInput=pspInput,
        context_left=context_left,
        context_right=context_right,
        retain_other_mods=retain_other_mods,
        return_unique=return_unique,
        return_sorted=return_sorted,
        organism=organism,
    )
    site_seq_context_annotator.load_annotations()
    df = site_seq_context_annotator.annotate(df)

    return df


def addSiteSequenceContext(
    df: pd.DataFrame,
    fastaFile: str,
    pspInput: bool = False,
    context_left: int = 15,
    context_right: int = 15,
    retain_other_mods: bool = False,
    return_unique: bool = False,
    return_sorted: bool = False,
    organism: str = "human",
) -> pd.DataFrame:
    """Annotate pandas dataframe with sequence context of a p-site.

    Adds the following annotation columns to dataframe\:

    - 'Site sequence context' = +/- 15 amino acids around each of the modified sites, separated by semicolons

    Required columns:
        :code:`Site positions`

    Args:
        df: pandas dataframe with 'Site positions' column
        fastaFile: fasta file containing protein sequences
        pspInput: set to True if fasta file was obtained from PhosphositePlus
        context_left: number of amino acids to the left of the modification to include
        context_right: number of amino acids to the right of the modification to include
        retain_other_mods: retain other modifications from the modified peptide in the sequence context in lower case
        return_unique: eliminate duplicated sequences from the 'Site sequence context' column, not preserving the order between the this column and the rest of the data frame
        return_sorted: sort the sequences from the 'Site sequence context' column alphabetically, not preserving the order between the this column and the rest of the data frame

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.SiteSequenceContextAnnotator(
        fastaFile,
        pspInput=pspInput,
        context_left=context_left,
        context_right=context_right,
        retain_other_mods=retain_other_mods,
        return_unique=return_unique,
        return_sorted=return_sorted,
        organism=organism,
    )
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addTurnoverRates(df: pd.DataFrame, turnoverFile: str) -> pd.DataFrame:
    """Annotate pandas dataframe with PTM turnover behavior.

    Adds column regarding the PTM turnover behavior.

    Adds the following annotation columns to dataframe\:

    - 'PTM Turnover' = rate of turnover for the modification sites according to Jana's PTM Turnover data

    Example:
        ::

            df = pa.addTurnoverRates(df, pa.turnoverFile)

    Required columns:
        :code:`Modified sequence`

    Args:
        df: pandas dataframe with 'Modified sequence' column
        turnoverFile: comma separated file with mapping from phosphosites to turnover information

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PTMTurnoverAnnotator(turnoverFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addPSPAnnotations(
    df: pd.DataFrame, phosphoSitePlusFile: str, organism: str = "human"
) -> pd.DataFrame:
    """Annotate pandas dataframe with number of high and low-throughput studies according to PhosphositePlus.

    Adds the following annotation columns to dataframe\:

    - LT_LIT = number of low-throughput studies
    - MS_LIT = number of high-throughput Mass Spec studies
    - MS_CST = number of high-throughput Mass Spec studies by CellSignalingTechnologies

    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addPSPAnnotations(df, pa.pspAnnotationFile)

    Required columns:
        :code:`Site positions`

    Args:
        df: pandas dataframe with 'Site positions' column
        phosphoSitePlusFile: tab separated file with PhosphositePlus annotations

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PSPStudiesAnnotator(phosphoSitePlusFile, organism=organism)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addPSPRegulatoryAnnotations(
    df: pd.DataFrame, phosphoSitePlusRegulatoryFile: str, organism: str = "human"
) -> pd.DataFrame:
    """Annotate pandas dataframe with regulatory functions according to PhosphositePlus.

    Adds the following annotation columns to dataframe\:

    - PSP_ON_FUNCTION = functional annotations for downstream regulation
    - PSP_ON_PROCESS = process annotations for downstream regulation
    - PSP_ON_PROT_INTERACT = protein interactions
    - PSP_ON_OTHER_INTERACT = other interactions
    - PSP_NOTES = regulatory site notes

    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addPSPRegulatoryAnnotations(df, pa.pspRegulatoryFile)

    Required columns:
        :code:`Site positions`

    Args:
        df: pandas dataframe with 'Site positions' column
        phosphoSitePlusRegulatoryFile: tab separated file with PhosphositePlus regulatory annotations

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PSPRegulatoryAnnotator(
        phosphoSitePlusRegulatoryFile, organism=organism
    )
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addPSPKinaseSubstrateAnnotations(
    df: pd.DataFrame,
    phosphoSitePlusKinaseSubstrateFile: str,
    gene_name: bool = False,
    organism: str = "human",
) -> pd.DataFrame:
    """Annotate pandas dataframe with upstream kinases according to PhosphositePlus.

    Adds the following annotation columns to dataframe\:

    - PSP Kinases = all phosphorylating kinases according to PhosphoSitePlus, no distinction is made between in vivo and in vitro evidence (this can be added in the future, if necessary)

    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addPSPKinaseSubstrateAnnotations(df, pa.pspKinaseSubstrateFile)

    Required columns:
        :code:`Site positions`

    Args:
        df: pandas dataframe with 'Site positions' column
        phosphoSitePlusKinaseSubstrateFile: tab separated file with PhosphositePlus kinase substrate relations
        gene_name: set to True to output the gene names instead of the kinase names

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PSPKinasesAnnotator(
        phosphoSitePlusKinaseSubstrateFile,
        output_gene_names=gene_name,
        organism=organism,
    )
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addDomains(df: pd.DataFrame, domainMappingFile: str) -> pd.DataFrame:
    """Adds column with domains the peptide overlaps with.

    Adds the following annotation columns to dataframe\:

    - Domains = semicolon separated list of domains that overlap with the peptide

    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addDomains(df, pa.domainMappingFile)

    Required columns:
        :code:`Matched proteins`, :code:`Start positions`, :code:`End positions`

    Args:
        df: pandas dataframe with 'Matched proteins', 'Start positions' and 'End positions' columns
        domainMappingFile: comma separated file with domains and their positions within the protein

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.DomainAnnotator(domainMappingFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addMotifs(df: pd.DataFrame, motifsFile: str) -> pd.DataFrame:
    """Adds column with motifs the site sequence context matches with.

    Adds the following annotation columns to dataframe\:

    - Motifs = semicolon separated list of motifs that match with the site sequence contexts

    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addMotifs(df, pa.motifsFile)

    Required columns:
        :code:`Site sequence context`

    Args:
        df: pandas dataframe with 'Site sequence context' column
        motifsFile: tab separated file with motifs and their identifiers

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.MotifAnnotator(motifsFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addInVitroKinases(
    df: pd.DataFrame, inVitroKinaseSubstrateMappingFile: str
) -> pd.DataFrame:
    """Annotate pandas dataframe with upstream in vitro kinases according to Sugiyama et al (2019).

    https://www.nature.com/articles/s41598-019-46385-4

    Adds the following annotation columns to dataframe\:

    - In Vitro Kinases = all phosphorylating kinases according to the Sugiyama in vitro kinase-substrate study


    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addInVitroKinases(df, pa.inVitroKinaseSubstrateMappingFile)

    Required columns:
        :code:`Site positions`

    Args:
        df: pandas dataframe with 'Site positions' column
        inVitroKinaseSubstrateMappingFile: tab separated file with in vitro kinase annotations

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.InVitroKinasesAnnotator(inVitroKinaseSubstrateMappingFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addKinaseLibraryAnnotations(
    df: pd.DataFrame,
    motifs_file: str,
    quantiles_file: str,
    top_n: int = 5,
    sort_type="total",
    threshold_type="total",
    score_cutoff: float = 3,
    split_sequences: bool = False,
) -> pd.DataFrame:
    """Annotate pandas dataframe with highest scoring kinases from the kinase library.

    Johnson et al. 2023, https://doi.org/10.1038/s41586-022-05575-3

    Requires "Site sequence context" column in the dataframe to be present.
    The "Site sequence context" column can be generated with PeptidePositionAnnotator()
    followed by SiteSequenceContextAnnotator().

    Adds the following annotation columns to dataframe\:

    - Motif Kinases = semicolon separated list of kinases that match with the site sequence contexts
    - Motif Scores = semicolon separated list of scores corresponding to Motif Kinases
    - Motif Percentiles = semicolon separated list of percentiles corresponding to Motif Kinases
    - Motif Totals = semicolon separated list of score*percentile corresponding to Motif Kinases

    Example:
        ::

            df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
            df = pa.addKinaseLibraryAnnotations(df, pa.kinaseLibraryMotifsFile, pa.kinaseLibraryQuantilesFile)

    Required columns:
        :code:`Site sequence context`

    Args:
        df: pandas dataframe with 'Site sequence context' column
        motifs_file: tab separated file with in odds ratios for each kinase, AA and position
        quantiles_file: tab separated file with quantile score for each kinase
        top_n: maximum number of returned kinases (default: 5)
        sort_type: score by which to sort the kinases, one of "percentile", "score" or "total" (default: "total")
        threshold_type: score to which to apply the cutoff, one of "percentile", "score" or "total" (total=score*percentile) (default: "total")
        score_cutoff: do not report kinases with a score below this cutoff (default: 3.0)
        split_sequences: if set to True, the 'Site sequence context' column is split by ';' and exploded before annotating

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.KinaseLibraryAnnotator(
        motifs_file,
        quantiles_file,
        top_n=top_n,
        sort_type=sort_type,
        threshold_type=threshold_type,
        score_cutoff=score_cutoff,
        split_sequences=split_sequences,
    )
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def main(argv):
    df = pd.read_csv(argv[0], sep="\t")
    addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True)

    from io import StringIO

    output = StringIO()
    df.to_csv(output)

    print(output.getvalue())


if __name__ == "__main__":
    main(sys.argv[1:])
