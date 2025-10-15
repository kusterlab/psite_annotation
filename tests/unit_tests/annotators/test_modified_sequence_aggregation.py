import pytest
import pandas as pd
import numpy as np

from psite_annotation.annotators import modified_sequence_aggregation


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
def delocalized_sequence_series():
    deloc_seqs = pd.Series(
        [
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_2",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
            "(ac)ASNSWNASSSPGEAREDGPEGLDK_1",
        ],
        name="Modified sequence",
    )
    return deloc_seqs


@pytest.fixture
def delocalized_sequence_df(
    modified_sequence_df: pd.DataFrame, delocalized_sequence_series: pd.Series
):
    modified_sequence_df.insert(
        loc=1, column="Delocalized sequence", value=delocalized_sequence_series
    )
    return modified_sequence_df


@pytest.fixture
def sequence_groups_series():
    delocalized_sequence_groups = pd.Series(
        [
            "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
            "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
            "(ac)ASNSWNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNSWNASS(ph)SPGEAREDGPEGLDK;(ac)ASNSWNASSS(ph)PGEAREDGPEGLDK",
            "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
            "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
            "(ac)ASNS(ph)WNAS(ph)SSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASS(ph)SPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSS(ph)PGEAREDGPEGLDK",
            "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
            "(ac)AS(ph)NSWNASSSPGEAREDGPEGLDK;(ac)ASNS(ph)WNASSSPGEAREDGPEGLDK",
        ],
        name="Modified sequence",
    )
    return delocalized_sequence_groups


@pytest.fixture
def sequence_groups_df(
    delocalized_sequence_df: pd.DataFrame, sequence_groups_series: pd.Series
):
    delocalized_sequence_df.insert(
        loc=2, column="Modified sequence group", value=sequence_groups_series
    )
    return delocalized_sequence_df


@pytest.fixture
def representative_sequence_df():
    representative_df = pd.DataFrame(
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
        }
    )
    return representative_df.set_index("Modified sequence group")


@pytest.fixture
def annotated_expected_df():
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


class TestFindRepresentativeModifiedSequence:
    def test_find_representative_modified_sequence(
        self, sequence_groups_df: pd.DataFrame, representative_sequence_df: pd.DataFrame
    ):
        experiment_cols = sequence_groups_df.columns[
            sequence_groups_df.columns.str.contains("Experiment")
        ]
        df_representative = (
            modified_sequence_aggregation.find_representative_modified_sequence(
                sequence_groups_df, experiment_cols
            )
        )
        pd.testing.assert_frame_equal(df_representative, representative_sequence_df)


class TestDelocalizationAnnotator:
    def test_delocalization_annotator(
        self, sequence_groups_df: pd.DataFrame, annotated_expected_df: pd.DataFrame
    ):
        annotator = modified_sequence_aggregation.ModifiedSequenceAggregatorAnnotator(
            experiment_cols=sequence_groups_df.columns[
                sequence_groups_df.columns.str.startswith("Experiment")
            ]
        )
        annotated_df = annotator.annotate(
            sequence_groups_df,
        )
        pd.testing.assert_frame_equal(
            annotated_df, annotated_expected_df, check_like=True
        )
