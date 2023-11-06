import logging
import sys

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
    "addPeptideAndPsitePositions",
    "addSiteSequenceContext",
    "addPSPAnnotations",
    "addPSPRegulatoryAnnotations",
    "addPSPKinaseSubstrateAnnotations",
    "addDomains",
    "addMotifs",
    "addInVitroKinases",
    "addTurnoverRates",
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


def addPeptideAndPsitePositions(
    df: pd.DataFrame,
    fastaFile: str,
    pspInput: bool = False,
    returnAllPotentialSites: bool = False,
) -> pd.DataFrame:
    """Annotate pandas dataframe with positions of the peptide within the protein sequence based on a fasta file.

    Adds the following annotation columns to dataframe:
    - 'Start positions' = starting positions of the modified peptide in the protein sequence (1-based, methionine is
      counted). If multiple isoforms/proteins contain the sequence, the starting positions are separated by
      semicolons in the same order as they are listed in the 'Proteins' input column
    - 'End positions' = end positions of the modified peptide in the protein sequence (see above for details)
    - 'Site positions' = position of the modification (see 'Start positions' above for details on how the position
      is counted)

    Args:
        df: pandas dataframe with "Proteins" and "Modified sequence" columns
        fastaFile: fasta file containing protein sequences
        pspInput: set to True if fasta file was obtained from PhosphositePlus
        returnAllPotentialSites: set to True if all S, T and Y positions should be returned as potention p-sites.

    Returns:
        pd.DataFrame: annotated dataframe

    """
    peptide_position_annotator = annotators.PeptidePositionAnnotator(
        fastaFile, pspInput, returnAllPotentialSites
    )
    peptide_position_annotator.load_annotations()
    df = peptide_position_annotator.annotate(df)

    annotator = annotators.SiteSequenceContextAnnotator(fastaFile, pspInput)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addSiteSequenceContext(
    df: pd.DataFrame, fastaFile: str, pspInput: bool = False
) -> pd.DataFrame:
    """Annotate pandas dataframe with sequence context of a p-site.

    Adds the following annotation columns to dataframe:
    - 'Site sequence context' = +/- 15 amino acids around each of the modified sites, separated by semicolons

    Args:
        df: pandas dataframe with 'Site positions' column
        fastaFile: fasta file containing protein sequences
        pspInput: set to True if fasta file was obtained from PhosphositePlus

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.SiteSequenceContextAnnotator(fastaFile, pspInput)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addTurnoverRates(df: pd.DataFrame, turnoverFile: str) -> pd.DataFrame:
    """Annotate pandas dataframe with PTM turnover behavior.

    Adds column regarding the PTM turnover behavior.

    Adds the following annotation columns to dataframe:
    - 'PTM Turnover' = rate of turnover for the modification sites according to Jana's PTM Turnover data

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


def addPSPAnnotations(df: pd.DataFrame, phosphoSitePlusFile: str) -> pd.DataFrame:
    """Annotate pandas dataframe with number of high and low-throughput studies according to PhosphositePlus.

    Adds the following annotation columns to dataframe:
    - LT_LIT = number of low-throughput studies
    - MS_LIT = number of high-throughput Mass Spec studies
    - MS_CST = number of high-throughput Mass Spec studies by CellSignalingTechnologies

    Args:
        df: pandas dataframe with 'Site positions' column
        phosphoSitePlusFile: tab separated file with PhosphositePlus annotations

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PSPStudiesAnnotator(phosphoSitePlusFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addPSPRegulatoryAnnotations(
    df: pd.DataFrame, phosphoSitePlusRegulatoryFile: str
) -> pd.DataFrame:
    """Annotate pandas dataframe with regulatory functions according to PhosphositePlus.

    Adds the following annotation columns to dataframe:
    - PSP_ON_FUNCTION = functional annotations for downstream regulation
    - PSP_ON_PROCESS = process annotations for downstream regulation
    - PSP_ON_PROT_INTERACT = protein interactions
    - PSP_ON_OTHER_INTERACT = other interactions
    - PSP_NOTES = regulatory site notes

    Args:
        df: pandas dataframe with 'Site positions' column
        phosphoSitePlusRegulatoryFile: tab separated file with PhosphositePlus regulatory annotations

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PSPRegulatoryAnnotator(phosphoSitePlusRegulatoryFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addPSPKinaseSubstrateAnnotations(
    df: pd.DataFrame, phosphoSitePlusKinaseSubstrateFile: str, gene_name: bool = False
):
    """Annotate pandas dataframe with upstream kinases according to PhosphositePlus.

    Adds the following annotation columns to dataframe:
    - PSP Kinases = all phosphorylating kinases according to PhosphoSitePlus, no distinction is made between in vivo
      and in vitro evidence (this can be added in the future, if necessary)

    Args:
        df: pandas dataframe with 'Site positions' column
        phosphoSitePlusKinaseSubstrateFile: tab separated file with PhosphositePlus kinase substrate relations
        gene_name: set to True to output the gene names instead of the kinase names

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.PSPKinasesAnnotator(
        phosphoSitePlusKinaseSubstrateFile, output_gene_names=gene_name
    )
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addDomains(df: pd.DataFrame, domainMappingFile: str):
    """Adds column with domains the peptide overlaps with.

    Adds the following annotation columns to dataframe:
    - Domains = semicolon separated list of domains that overlap with the peptide

    Args:
        df: pandas dataframe with 'Proteins', 'Start positions' and 'End positions' columns
        domainMappingFile: comma separated file with domains and their positions within the protein

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.DomainAnnotator(domainMappingFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addMotifs(df: pd.DataFrame, motifsFile: str):
    """Adds column with motifs the site sequence context matches with.

    Adds the following annotation columns to dataframe:
    - Motifs = semicolon separated list of motifs that match with the site sequence contexts

    Args:
        df: pandas dataframe with 'Proteins', 'Start positions' and 'End positions' columns
        motifsFile: tab separated file with motifs and their identifiers

    Returns:
        pd.DataFrame: annotated dataframe

    """
    annotator = annotators.MotifAnnotator(motifsFile)
    annotator.load_annotations()
    df = annotator.annotate(df)

    return df


def addInVitroKinases(df: pd.DataFrame, inVitroKinaseSubstrateMappingFile: str):
    """Annotate pandas dataframe with upstream in vitro kinases according to Sugiyama et al (2019).

    https://www.nature.com/articles/s41598-019-46385-4

    Adds the following annotation columns to dataframe:
    - In Vitro Kinases = all phosphorylating kinases according to the Sugiyama in vitro kinase-substrate study

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


def main(argv):
    df = pd.read_csv(argv[0], sep="\t")
    addPeptideAndPsitePositions(df, pspFastaFile, pspInput=True)

    from io import StringIO

    output = StringIO()
    df.to_csv(output)

    print(output.getvalue())


if __name__ == "__main__":
    main(sys.argv[1:])
