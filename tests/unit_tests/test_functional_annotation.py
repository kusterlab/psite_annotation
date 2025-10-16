import unittest
import unittest.mock

import pytest
import pandas as pd
import numpy as np

import psite_annotation.functional_annotation as pa


mock_input_file_psp = """;110220
;Data extracted from PhosphoSitePlus(R), created by Cell Signaling Technology Inc. PhosphoSitePlus is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. Attribution must be given in written, oral and digital presentations to PhosphoSitePlus, www.phosphosite.org. Written documents should additionally cite Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, Sullivan M (2012) PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse. Nucleic Acids Res. 40, D261<D0>70.; www.phosphosite.org.

>GN:Cbln1|CBLN1|mouse|Q9R171
MLGVVELLLLGTAWLAGPARGQNETEPIVLEGKCLVVCDSNPTSDPTGTALGISVRSGSA
KVAFSAIRSTNHEPSEMSNRTMIIYFDQVLVNIGNNFDSERSTFIAPRKGIYSFNFHVVK
VYNRQTIQVSLMLNGWPVISAFAGDQDVTREAASNGVLIQMEKGDRAYLKLERGNLMGGW
KYSTFSGFLVFPL
>GN:Cox7a2|COX7A2|mouse|P48771
MLRNLLALRQIAQRTISTTSRRHFENKVPEKQKLFQEDNGMPVHLKGGASDALLYRATMA
LTLGGTAYAIYLLAMAAFPKKQN
>GN:RPS20|RPS20|human|P60866
MAFKDTGKTPVEPEVAIHRIRITLTSRNVKSLEKVCADLIRGAKEKNLKVKGPVRMPTKT
LRITTRKTPCGEGSKTWDRFQMRIHKRLIDLHSPSEIVKQITSISIEPGVEVEVTIADA
>GN:UFC1|UFC1|human|Q9Y3C8
MADEATRRVVSEIPVLKTNAGPRDRELWVQRLKEEYQSLIRYVENNKNADNDWFRLESNK
EGTRWFGKCWYIHDLLKYEFDIEFDIPITYPTTAPEIAVPELDGKTAKMYRGGKICLTDH
FKPLWARNVPKFGLAHLMALGLGPWLAVEIPDLIQKGVIQHKEKCNQ
"""  # noqa: E501, B950


class TestAddSiteSequenceContext:
    """Test the addSiteSequenceContext function."""

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file_psp),
        create=True,
    )
    def test_add_site_sequence_context(self):
        df = pd.DataFrame(
            {"Site positions": ["Q9Y3C8_T6", "Q9Y3C8_Q167", "Q9Y3C8_Q168"]}
        )
        result_df = pa.addSiteSequenceContext(df, "mock_input_file_psp", pspInput=True)

        assert (
            result_df["Site sequence context"].iloc[0]
            == "__________MADEAtRRVVSEIPVLKTNAG"
        )
        assert (
            result_df["Site sequence context"].iloc[1]
            == "DLIQKGVIQHKEKCNq_______________"
        )
        assert result_df["Site sequence context"].iloc[2] == ""


@pytest.fixture
def modified_sequence_df() -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "Modified sequence": [
                "(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK",
                "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
                "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK",
            ],
            "Experiment 1": [10.0, np.nan, np.nan, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 2": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 3": [np.nan, 10.0, np.nan, 10.0, np.nan, np.nan, np.nan, 5.0],
            "Experiment 4": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 5": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 6": [np.nan, 10.0, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 7": [np.nan, np.nan, 10.0, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 8": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 9": [10.0, np.nan, np.nan, np.nan, np.nan, 10.0, 5.0, np.nan],
            "Experiment 10": [np.nan, 10.0, np.nan, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 11": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 12": [10.0, np.nan, np.nan, np.nan, np.nan, 10.0, 5.0, np.nan],
            "Experiment 13": [np.nan, 10.0, np.nan, 10.0, 10.0, np.nan, 5.0, np.nan],
            "Experiment 14": [np.nan, np.nan, 10.0, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 15": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 16": [10.0, 10.0, np.nan, 10.0, np.nan, np.nan, np.nan, 5.0],
            "Experiment 17": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 18": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 19": [np.nan, np.nan, 10.0, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 20": [10.0, 10.0, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
        }
    )
    return df


@pytest.fixture
def annotated_expected_df() -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "Modified sequence": [
                "(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK",
                "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
                "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK",
            ],
            "Experiment 1": [10.0, np.nan, np.nan, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 2": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 3": [np.nan, 10.0, np.nan, 10.0, np.nan, np.nan, np.nan, 5.0],
            "Experiment 4": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 5": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 6": [np.nan, 10.0, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 7": [np.nan, np.nan, 10.0, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 8": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 9": [10.0, np.nan, np.nan, np.nan, np.nan, 10.0, 5.0, np.nan],
            "Experiment 10": [np.nan, 10.0, np.nan, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 11": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 12": [10.0, np.nan, np.nan, np.nan, np.nan, 10.0, 5.0, np.nan],
            "Experiment 13": [np.nan, 10.0, np.nan, 10.0, 10.0, np.nan, 5.0, np.nan],
            "Experiment 14": [np.nan, np.nan, 10.0, np.nan, 10.0, np.nan, 5.0, np.nan],
            "Experiment 15": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 16": [10.0, 10.0, np.nan, 10.0, np.nan, np.nan, np.nan, 5.0],
            "Experiment 17": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 18": [10.0, np.nan, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 19": [np.nan, np.nan, 10.0, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Experiment 20": [10.0, 10.0, np.nan, 10.0, np.nan, np.nan, 5.0, np.nan],
            "Delocalized sequence": [
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
            ],
            "Modified sequence group": [
                "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
                "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
            ],
        }
    )
    return df


@pytest.fixture
def aggregated_expected_df():
    return pd.DataFrame(
        {
            "Modified sequence group": [
                "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
            ],
            "Modified sequence representative": [
                "(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
                "(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
                "(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
            ],
            "Modified sequence representative degree": [
                0.9,
                0.6666666666666666,
                0.5909090909090909,
            ],
            "Delocalized sequence": [
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
                "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
            ],
            "Experiment 1": [5.0, 10.0, 10.0],
            "Experiment 2": [5.0, 10.0, 10.0],
            "Experiment 3": [5.0, 10.0, 10.0],
            "Experiment 4": [5.0, 10.0, 10.0],
            "Experiment 5": [5.0, 10.0, 10.0],
            "Experiment 6": [5.0, 10.0, 10.0],
            "Experiment 7": [5.0, 10.0, 10.0],
            "Experiment 8": [5.0, 10.0, 10.0],
            "Experiment 9": [5.0, 10.0, 10.0],
            "Experiment 10": [5.0, 10.0, 10.0],
            "Experiment 11": [5.0, 10.0, 10.0],
            "Experiment 12": [5.0, 10.0, 10.0],
            "Experiment 13": [5.0, 10.0, 10.0],
            "Experiment 14": [5.0, 10.0, 10.0],
            "Experiment 15": [5.0, 10.0, 10.0],
            "Experiment 16": [5.0, 10.0, 10.0],
            "Experiment 17": [5.0, 10.0, 10.0],
            "Experiment 18": [5.0, 10.0, 10.0],
            "Experiment 19": [5.0, 10.0, 10.0],
            "Experiment 20": [5.0, 10.0, 10.0],
        }
    )


class TestAddModifiedSequenceGroups:
    def test_aggregate_modified_sequence_groups(
        self, modified_sequence_df: pd.DataFrame, annotated_expected_df: pd.DataFrame
    ):
        annotated_df = pa.addModifiedSequenceGroups(
            modified_sequence_df,
        )
        pd.testing.assert_frame_equal(annotated_df, annotated_expected_df)


class TestAggregateModifiedSequenceGroups:
    def test_aggregate_modified_sequence_groups(
        self, modified_sequence_df: pd.DataFrame, aggregated_expected_df: pd.DataFrame
    ):
        annotated_df = pa.aggregateModifiedSequenceGroups(
            modified_sequence_df,
            experiment_cols=modified_sequence_df.columns[
                modified_sequence_df.columns.str.startswith("Experiment")
            ],
        )
        pd.testing.assert_frame_equal(annotated_df, aggregated_expected_df)

    def test_aggregate_modified_sequence_groups_extra_columns(
        self, modified_sequence_df: pd.DataFrame, aggregated_expected_df: pd.DataFrame
    ):
        modified_sequence_df["Gene Names"] = [
            "GeneA",
            "GeneA;GeneB",
            "GeneB",
            "GeneB",
            "GeneA;GeneB",
            "GeneB",
            "GeneC",
            "GeneA;GeneB",
        ]
        annotated_df = pa.aggregateModifiedSequenceGroups(
            modified_sequence_df,
            experiment_cols=modified_sequence_df.columns[
                modified_sequence_df.columns.str.startswith("Experiment")
            ],
            agg_cols={"Gene Names": "first"},
        )
        aggregated_expected_df["Gene Names"] = ["GeneC", "GeneB", "GeneA"]
        pd.testing.assert_frame_equal(annotated_df, aggregated_expected_df)
