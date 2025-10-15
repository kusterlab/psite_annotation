import pytest
import pandas as pd
import numpy as np

from psite_annotation.annotators import modified_sequence_group


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


class TestDelocalizeSequence:
    def test_delocalize_sequence(
        self, modified_sequence_df: pd.DataFrame, delocalized_sequence_series: pd.Series
    ):
        deloc_seqs = modified_sequence_group.delocalize_phospho_sequence(
            modified_sequence_df["Modified sequence"]
        )
        pd.testing.assert_series_equal(deloc_seqs, delocalized_sequence_series)


class TestAggregateGroups:
    def test_aggregate_ptm_groups(
        self, delocalized_sequence_df: pd.DataFrame, sequence_groups_series: pd.Series
    ):
        sequence_groups = modified_sequence_group.aggregate_phospho_groups(
            delocalized_sequence_df, match_tolerance=2
        )
        pd.testing.assert_series_equal(sequence_groups, sequence_groups_series)


class TestDelocalizationAnnotator:
    def test_delocalization_annotator(
        self, modified_sequence_df: pd.DataFrame, sequence_groups_df: pd.DataFrame
    ):
        annotator = modified_sequence_group.ModifiedSequenceGroupAnnotator(
            match_tolerance=2
        )
        annotated_df = annotator.annotate(modified_sequence_df)
        pd.testing.assert_frame_equal(annotated_df, sequence_groups_df)


class TestMakeMonophosVersions:
    """Unit tests for make_monophos_versions function."""

    def test_multiple_sites(self):
        mod_seq = "(ph)ACDE(ph)FGHIK"
        expected = [
            "(ph)ACDEFGHIK",  # phosphorylation at position 0
            "ACDE(ph)FGHIK",  # phosphorylation at position 4
        ]
        assert modified_sequence_group.make_monophos_versions(mod_seq) == expected

    def test_no_sites(self):
        mod_seq = "ACDEFGHIK"
        assert modified_sequence_group.make_monophos_versions(mod_seq) == []

    def test_single_site(self):
        mod_seq = "ACDE(ph)FGHIK"
        expected = ["ACDE(ph)FGHIK"]
        assert modified_sequence_group.make_monophos_versions(mod_seq) == expected

    def test_ordering(self):
        mod_seq = "A(ph)C(ph)D(ph)E"
        expected = [
            "A(ph)CDE",  # position 1
            "AC(ph)DE",  # position 2
            "ACD(ph)E",  # position 3
        ]
        assert modified_sequence_group.make_monophos_versions(mod_seq) == expected
