import json
import logging

from config_path import ConfigPath

logger = logging.getLogger(__name__)

default_config_json = """{
    "domainMappingFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/uniprot_to_domain.csv",
    "inVitroKinaseSubstrateMappingFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/yasushi_supp_table2_kinase_substrate_relations_mapped_ids.tsv",
    "motifsFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/motifs_all.tsv",
    "turnoverFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/TurnoverSites.csv",
    "pspFastaFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/PhosphoSitePlus/Phosphosite_seq.fasta",
    "pspKinaseSubstrateFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/PhosphoSitePlus/Kinase_Substrate_Dataset",
    "pspAnnotationFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/PhosphoSitePlus/Phosphorylation_site_dataset",
    "pspRegulatoryFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/PhosphoSitePlus/Regulatory_sites",
    "kinaseLibraryMotifsFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/Motif_Odds_Ratios.txt",
    "kinaseLibraryQuantilesFile": "Z:/line_functions/bioinformatics/Databases/psite_annotation_mapping/Kinase_Score_Quantile_Matrix.txt"
}"""


def _getUserConfigPath():
    return ConfigPath("psite_annotation", "mls.ls.tum.de", ".json")


def _getConfigDicts():
    config = _getUserConfigPath()
    config_file = config.readFilePath()
    if config_file is None:
        json_config = default_config_json
        logger.info("Using default builtin config")
    else:
        logger.info(f"Using config from {config_file}")
        with open(config_file) as f:
            json_config = f.read()

    defaults = json.loads(default_config_json)
    user = json.loads(json_config)
    
    return defaults, user


def setUserConfig(config_path_user: str):
    with open(config_path_user) as f:
        config_user = json.load(f)

    config = _getUserConfigPath()
    config_path = config.saveFilePath(mkdir=True)

    logger.info(f"Overwriting user's config.json located at {config_path}")

    with open(config_path, "w") as f:
        json.dump(config_user, f, indent=4)


def _getConfigSetting(name, user, defaults):
    if name in user:
        return user[name]
    else:
        return defaults[name]
