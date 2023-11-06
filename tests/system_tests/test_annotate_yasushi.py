import pandas as pd

import psite_annotation as pa


def test_annotateYasushi():
    """Test if the site position strings in the Yasushi data table can be mapped to site sequence contexts."""
    df = pd.read_csv(
        "./tests/data/yasushi_subset.tsv",
        sep="\t",
    )
    df = df.rename(columns={"UniprotId_Site": "Site positions"})

    pa.addSiteSequenceContext(df, pa.pspFastaFile, pspInput=True)

    assert df.loc[0, "Site sequence context"] == "AFDEAIAELDTLNEEsYKDSTLIMQLLRDNL"


if __name__ == "__main__":
    test_annotateYasushi()
