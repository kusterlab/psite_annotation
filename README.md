# Psite annotation

Python module for annotating a pandas dataframe with phosphosites, e.g. PhosphoSitePlus annotations, kinase-substrate relations, domain information, etc.

## Installation

If you have setup SSH keys in gitlab, you can easily install this package with:

```
pip install psite-annotation
```

Otherwise, you can clone this repository and install it with pip manually:

```
git clone https://www.github.com/kusterlab/psite_annotation.git
cd psite_annotation
pip install .
```

### Installing annotation files

The easiest way to use the package is to supply the annotation files each time you call the respective functions.
Alternatively, it is also possible to install a configuration file to automatically point to the annotation files.

To set the paths to the annotation files, create a `config.json` file with the following content:

```
{
    "domainMappingFile": "/path/to/uniprot_to_domain.csv",
    "inVitroKinaseSubstrateMappingFile": "/path/to/yasushi_supp_table2_kinase_substrate_relations_mapped_ids.tsv",
    "motifsFile": "/path/to/motifs_all.tsv",
    "turnoverFile": "/path/to/TurnoverSites.csv",
    "kinaseLibraryMotifsFile": "/path/to/Motif_Odds_Ratios.txt",
    "kinaseLibraryQuantilesFile": "/path/to/Kinase_Score_Quantile_Matrix.txt",
    "pspFastaFile": "/path/to/Phosphosite_seq.fasta",
    "pspKinaseSubstrateFile": "/path/to/Kinase_Substrate_Dataset",
    "pspAnnotationFile": "/path/to/Phosphorylation_site_dataset",
    "pspRegulatoryFile": "/path/to/Regulatory_sites"
}
```

Where `/path/to` should be replaced by the absolute path to the annotation files.

For Kusterlab internal users:

- ask Matthew for the config file.

For external users:

- PhosphoSitePlus annotation files can be downloaded from https://www.phosphosite.org/staticDownloads.action (account needed).
- The other annotation files are available from the `annotations.zip` file in this repository.

You can then install this config file with:

```
python -c "import psite_annotation.config as c; c.setUserConfig('./config.json')"
```

## Usage

To add upstream kinases to a pandas dataframe `df` with columns `Proteins` (UniProt identifiers separated by semicolons, e.g. `Q86U42-2;Q86U42`) and `Modified sequence` (standard MaxQuant notation, e.g. `(ac)AAAAAAAAAAGAAGGRGS(ph)GPGR`):

```
import psite_annotation as pa

df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
df = pa.addPSPKinaseSubstrateAnnotations(df, pa.pspKinaseSubstrateFile)
```

We first need to annotate the peptide and modification positions within the protein using `addPeptideAndPsitePositions()`. This adds a column with an identifier for the phosphosite, which can then be mapped to the phosphorylating kinase using `addPSPKinaseSubstrateAnnotations()`.

## Functions

### pa.addPeptideAndPsitePositions()

```
input: pandas dataframe with 'Proteins' (Usually the Uniprot ID) and 'Modified sequence' columns
output: pandas dataframe with the following added columns:
- 'Start positions' = starting positions of the modified peptide in the protein sequence (1-based, methionine is counted). If multiple isoforms/proteins contain the sequence, the starting positions are separated by semicolons in the same order as they are listed in the 'Proteins' input column
- 'End positions' = end positions of the modified peptide in the protein sequence (see above for details)
- 'Site sequence context' = +/- 15 amino acids around each of the modified sites, separated by semicolons
- 'Site positions' = position of the modification (see 'Start positions' above for details on how the position is counted)
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, fastaFile)
```

### pa.addPSPAnnotations()

```
input: pandas dataframe with 'Site positions' column (this column can be obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- PSP_LT_LIT = number of low-throughput studies
- PSP_MS_LIT = number of high-throughput Mass Spec studies
- PSP_MS_CST = number of high-throughput Mass Spec studies by CellSignalingTechnologies
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
df = pa.addPSPAnnotations(df, pa.pspAnnotationFile)
```

### pa.addPSPRegulatoryAnnotations()

```
input: pandas dataframe with 'Site positions' column (this column can be obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- PSP_ON_FUNCTION = functional annotations for downstream regulation
- PSP_ON_PROCESS = process annotations for downstream regulation
- PSP_ON_PROT_INTERACT = protein interactions
- PSP_ON_OTHER_INTERACT = other interactions
- PSP_NOTES = regulatory site notes
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
df = pa.addPSPRegulatoryAnnotations(df, pa.pspRegulatoryFile)
```

### pa.addPSPKinaseSubstrateAnnotations()

```
input: pandas dataframe with 'Site positions' column (this column can be obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- 'PSP Kinases' = all phosphorylating kinases according to PhosphoSitePlus, no distinction is made between in vivo and in vitro evidence (this can be added in the future, if necessary)
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
df = pa.addPSPKinaseSubstrateAnnotations(df, pa.pspKinaseSubstrateFile)
```

### pa.addDomains()

```
input: pandas dataframe with 'Proteins', 'Start positions' and 'End positions' columns (the latter two are obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- Domains = domains that overlap with the modified peptide sequence
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, fastaFile)
df = pa.addDomains(df, pa.domainMappingFile)
```

### pa.addMotifs()

```
input: pandas dataframe with 'Site sequence context' columns (can be obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- Motifs = matching motif identifiers for all of the modified sites
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, fastaFile)
df = pa.addMotifs(df, pa.motifsFile)
```

### pa.addInVitroKinases()

```
input: pandas dataframe with 'Site positions' column (this column can be obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- 'In Vitro Kinases' = all phosphorylating kinases according to the Yasushi in vitro kinase-substrate study
```

Usage:

```
df = pa.addPeptideAndPsitePositions(df, fastaFile)
df = pa.addInVitroKinases(df, pa.inVitroKinaseSubstrateMappingFile)
```

### pa.addTurnoverRates()

```
input: pandas dataframe with 'Modified sequence' column
output: pandas dataframe with the following added columns:
- 'PTM Turnover' = rate of turnover for the modification sites according to Jana's PTM Turnover data, e.g. slower, faster
```

Usage:

```
df = pa.addTurnoverRates(df, pa.turnoverFile)
```

### pa.addKinaseLibraryAnnotations()

```
input: pandas dataframe with 'Site sequence context' columns (can be obtained from the addPeptideAndPsitePositions() function)
output: pandas dataframe with the following added columns:
- Motif Kinases = semicolon separated list of kinases that match with the site sequence contexts
- Motif Scores = semicolon separated list of scores corresponding to Motif Kinases
- Motif Percentiles = semicolon separated list of percentiles corresponding to Motif Kinases
- Motif Totals = semicolon separated list of score*percentile corresponding to Motif Kinases
```

Usage:

```
df = pa.addTurnoverRates(df, pa.turnoverFile)
```
