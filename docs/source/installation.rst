.. highlight:: shell

Installation
============

Prerequisites
~~~~~~~~~~~~~

Picked Group FDR requires python >=3.9,<=3.11.

Using pip (recommended)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install psite-annotation

Alternatively, you can clone this repository and install it with pip manually:

.. code-block:: bash

    git clone https://github.com/kusterlab/psite_annotation.git
    cd psite_annotation
    pip install .


Installing annotation files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to use the package is to supply the annotation files each time you call the respective functions.
Alternatively, it is also possible to install a configuration file to automatically point to the annotation files.

To set the paths to the annotation files, create a `config.json` file with the following content:

.. code-block:: json

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

Where `/path/to` should be replaced by the absolute path to the annotation files.

For Kusterlab internal users:

- ask Matthew for the config file.

For external users:

- PhosphoSitePlus annotation files can be downloaded from https://www.phosphosite.org/staticDownloads.action (account needed).
- The other annotation files are available from the `annotations.zip` file in this repository.

You can then install this config file with:

.. code-block:: bash

    python -c "import psite_annotation.config as c; c.setUserConfig('./config.json')"
