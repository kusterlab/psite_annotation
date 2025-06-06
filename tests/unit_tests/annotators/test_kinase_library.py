import io

import pandas as pd
import pytest
from psite_annotation.annotators.annotator_base import MissingColumnsError

from psite_annotation.annotators.kinase_library import KinaseLibraryAnnotator, _score


@pytest.fixture
def annotator() -> KinaseLibraryAnnotator:
    """Fixture to create a MotifAnnotator object using a mock input file.

    Returns:
        MotifAnnotator: Annotator object with mock file loaded
    """
    motifs_file = io.StringIO(mock_motifs_input_file)
    quantiles_file = io.StringIO(mock_quantiles_input_file)
    return KinaseLibraryAnnotator(motifs_file, quantiles_file, score_cutoff=-1)


@pytest.fixture
def input_df() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE",
                "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG",
                "AAAAAAAGAAGGRGStYGPGRRRHLVPGAGG",
                "AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE",
                "AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE",
            ],
        }
    )


@pytest.fixture
def expected_output_df() -> pd.DataFrame:
    """Fixture to create an example ouput dataframe to the annotate() method given the input_df fixture.

    Returns:
        pd.DataFrame: Output dataframe corresponding to input_df fixture

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE",
                "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG",
                "AAAAAAAGAAGGRGStYGPGRRRHLVPGAGG",
                "AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE",
                "AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE",
            ],
            "Motif Kinases": [
                "AAK1;ACVR2A",
                "ACVR2A",
                "ACVR2A",
                "",
                "",
            ],
            "Motif Scores": [
                "0.515;-2.848",
                "-1.365",
                "-1.174",
                "",
                "",
            ],
            "Motif Percentiles": [
                "0.946;0.261",
                "0.638",
                "0.675",
                "",
                "",
            ],
            "Motif Totals": [
                "0.487;-0.745",
                "-0.871",
                "-0.792",
                "",
                "",
            ],
        }
    )

@pytest.fixture
def input_df_exact_context() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "AGGRGsGPGRR",
                "AGGRGsTYGPG",
                "GGRGStYGPGR",
                "GRGSTyGPGRR",
                "GRGSTyGPGRR",
            ],
        }
    )


@pytest.fixture
def expected_output_df_exact_context() -> pd.DataFrame:
    """Fixture to create an example ouput dataframe to the annotate() method given the input_df fixture.

    Returns:
        pd.DataFrame: Output dataframe corresponding to input_df fixture

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "AGGRGsGPGRR",
                "AGGRGsTYGPG",
                "GGRGStYGPGR",
                "GRGSTyGPGRR",
                "GRGSTyGPGRR",
            ],
            "Motif Kinases": [
                "AAK1;ACVR2A",
                "ACVR2A",
                "ACVR2A",
                "",
                "",
            ],
            "Motif Scores": [
                "0.515;-2.848",
                "-1.365",
                "-1.174",
                "",
                "",
            ],
            "Motif Percentiles": [
                "0.946;0.261",
                "0.638",
                "0.675",
                "",
                "",
            ],
            "Motif Totals": [
                "0.487;-0.745",
                "-0.871",
                "-0.792",
                "",
                "",
            ],
        }
    )

@pytest.fixture
def input_df_too_short_context() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "GRGsGPG",
                "GRGsTYG",
                "RGStYGP",
                "GSTyGPG",
                "GSTyGPG",
            ],
        }
    )


@pytest.fixture
def expected_output_df_too_short_context() -> pd.DataFrame:
    """Fixture to create an example ouput dataframe to the annotate() method given the input_df fixture.

    Returns:
        pd.DataFrame: Output dataframe corresponding to input_df fixture

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "GRGsGPG",
                "GRGsTYG",
                "RGStYGP",
                "GSTyGPG",
                "GSTyGPG",
            ],
            "Motif Kinases": ["AAK1;ACVR2A", "ACVR2A", "ACVR2A", "", ""],
            "Motif Scores": ["0.782;-2.476", "-1.619", "-0.67", "", ""],
            "Motif Percentiles": ["0.952;0.354", "0.581", "0.752", "", ""],
            "Motif Totals": ["0.745;-0.876", "-0.94", "-0.504", "", ""],
        }
    )

@pytest.fixture
def input_df_multiple_seq() -> pd.DataFrame:
    """Fixture to create an example input dataframe with multiple ';'-separated site sequence contexts

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE;AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG;AAAAAAAGAAGGRGStYGPGRRRHLVPGAGG",
                "AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE;AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE",
            ],
        }
    )


class TestKinaseLibraryAnnotator:
    """Test the KinaseLibraryAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file into `odds_df`.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        assert len(annotator.odds_dict.keys()) == 420
        assert ("ACVR2A", 4, "y") in annotator.odds_dict.keys()

        assert annotator.odds_dict[("ACVR2A", 4, "y")] == 1.0015

        assert len(annotator.quantiles.keys()) == 2

    def test_raise_error_on_missing_site_sequence_context_column(self, annotator):
        """Test that a value error is raised if the "Site sequence context" column is missing.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        input_df = pd.DataFrame(
            {
                "Sequence": [
                    "AAAAAAAAG",
                    "AAAAAAAAG",
                ],
            }
        )
        with pytest.raises(MissingColumnsError):
            annotator.annotate(input_df)

    def test_annotate(self, annotator, input_df, expected_output_df):
        """Test that the annotate method correctly adds the motif annotations to the input dataframe as new column.

        Args:
            annotator: Annotator object with mock file loaded
            input_df: Example input dataframe
            expected_output_df: Expected output dataframe

        """
        annotator.load_annotations()

        # Annotate the input dataframe
        output_df = annotator.annotate(input_df)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Site sequence context",
            "Motif Kinases",
            "Motif Scores",
            "Motif Percentiles",
            "Motif Totals",
        }

        # Assert that the output dataframe has the expected values

        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)

    def test_annotate_exact_context(
        self,
        annotator,
        input_df_exact_context,
        expected_output_df_exact_context,
    ):
        """Test that the annotate method correctly adds the motif annotations to the input dataframe as new column.

        Args:
            annotator: Annotator object with mock file loaded
            input_df: Example input dataframe
            expected_output_df: Expected output dataframe

        """
        annotator.load_annotations()

        # Annotate the input dataframe
        output_df = annotator.annotate(input_df_exact_context)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Site sequence context",
            "Motif Kinases",
            "Motif Scores",
            "Motif Percentiles",
            "Motif Totals",
        }

        # Assert that the output dataframe has the expected values
        pd.testing.assert_frame_equal(
            output_df, expected_output_df_exact_context, check_like=True
        )

    def test_annotate_too_short_context(
        self,
        annotator,
        input_df_too_short_context,
        expected_output_df_too_short_context,
    ):
        """Test that the annotate method correctly adds the motif annotations to the input dataframe as new column.

        Args:
            annotator: Annotator object with mock file loaded
            input_df: Example input dataframe
            expected_output_df: Expected output dataframe

        """
        annotator.load_annotations()

        # Annotate the input dataframe
        output_df = annotator.annotate(input_df_too_short_context)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Site sequence context",
            "Motif Kinases",
            "Motif Scores",
            "Motif Percentiles",
            "Motif Totals",
        }

        # Assert that the output dataframe has the expected values
        pd.testing.assert_frame_equal(
            output_df, expected_output_df_too_short_context, check_like=True
        )

    def test_annotate_input_df_unchanged(self, annotator, input_df):
        """Test that the annotate method does not alter the input dataframe.

        Args:
            annotator: Annotator object with mock file loaded
            input_df: Example input dataframe

        """
        input_df_copy = input_df.copy()

        annotator.load_annotations()

        # Annotate the input dataframe
        _ = annotator.annotate(input_df)

        pd.testing.assert_frame_equal(input_df, input_df_copy, check_like=True)

    def test_annotate_inplace(self, annotator, input_df, expected_output_df):
        """Test that the annotate method can annotate the input dataframe inplace.

        Args:
            annotator: Annotator object with mock file loaded
            input_df: Example input dataframe
            expected_output_df: Expected output dataframe

        """
        annotator.load_annotations()

        # Annotate the input dataframe
        annotator.annotate(input_df, inplace=True)

        pd.testing.assert_frame_equal(input_df, expected_output_df, check_like=True)

    def test_raise_error_on_illegal_char(self, annotator, input_df_multiple_seq):
        """Test that an error is raised if the "Site sequence context" column contains non-AA characters.
        If 'split_sequences' is not set to True, the ';' is considered an illegal character.
        Args:
            annotator: Annotator object with mock file loaded
            input_df_multiple_seq: Example input dataframe with ';' in the "Site sequence context" column
        """
        annotator.load_annotations()

        with pytest.raises(ValueError):
            annotator.annotate(input_df_multiple_seq)

    def test_annotate_split_multiple_seq(self, annotator, input_df_multiple_seq, expected_output_df):
        """Test that the annotate method correctly splits up ';'-separated sequences before annotation, if requested.

        Args:
            annotator: Annotator object with mock file loaded
            input_df_multiple_seq: Example input dataframe with ';'-split "Site sequence context"
            expected_output_df: Expected output dataframe

        """
        annotator.split_sequences = True

        annotator.load_annotations()

        output_df = annotator.annotate(input_df_multiple_seq)

        #In order for the output to be equal to the expected output, we need to drop the index
        output_df.reset_index(inplace=True, drop=True)

        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)

    def test_raise_error_on_too_short_sequence(self):
        """Test that the _score function throws an error if the submitted sequence does not have the length:
         (2 * motif_size + 1)"""
        with pytest.raises(AssertionError):
            _score('SHORT', None, None, 5)

    def test_annotate_empty_seq(self, annotator):
        """Test that the annotate method correctly annotated empty input sequences with no scores

        Args:
            annotator: Annotator object with mock file loaded

        """

        df = pd.DataFrame({'Site sequence context': ['', '_']})
        annotator.load_annotations()

        output_df = annotator.annotate(df)

        expected_empty_output_df = pd.DataFrame({
            "Site sequence context": ['', '_'],
            "Motif Kinases": ['', ''],
            "Motif Scores": ['', ''],
            "Motif Percentiles": ['', ''],
            "Motif Totals": ['', ''],
        })

        pd.testing.assert_frame_equal(output_df, expected_empty_output_df, check_like=True)


# Define the mock input file as a string
mock_motifs_input_file = """Kinase	Position	AA	Odds Ratio
AAK1	-1	A	0.6544
AAK1	-1	C	0.5328
AAK1	-1	D	0.2744
AAK1	-1	E	0.2484
AAK1	-1	F	0.9257
AAK1	-1	G	1.2002
AAK1	-1	H	1.17
AAK1	-1	I	0.7387
AAK1	-1	K	2.0893
AAK1	-1	L	1.1116
AAK1	-1	M	0.9412
AAK1	-1	N	0.7371
AAK1	-1	P	1.9305
AAK1	-1	Q	0.6751
AAK1	-1	R	1.5694
AAK1	-1	S	0.9257
AAK1	-1	T	0.9257
AAK1	-1	V	0.6659
AAK1	-1	W	0.7476
AAK1	-1	Y	1.3205
AAK1	-1	s	0.2057
AAK1	-1	t	0.2057
AAK1	-1	y	0.3033
AAK1	-2	A	0.9151
AAK1	-2	C	1.0
AAK1	-2	D	1.3679
AAK1	-2	E	1.0816
AAK1	-2	F	0.5713
AAK1	-2	G	0.4803
AAK1	-2	H	0.9634
AAK1	-2	I	0.6424
AAK1	-2	K	0.9984
AAK1	-2	L	0.9678
AAK1	-2	M	1.9828
AAK1	-2	N	1.3738
AAK1	-2	P	0.3848
AAK1	-2	Q	2.5792
AAK1	-2	R	0.9432
AAK1	-2	S	0.9432
AAK1	-2	T	0.9432
AAK1	-2	V	0.5784
AAK1	-2	W	0.4917
AAK1	-2	Y	0.678
AAK1	-2	s	0.356
AAK1	-2	t	0.356
AAK1	-2	y	0.8261
AAK1	-3	A	1.9015
AAK1	-3	C	1.1138
AAK1	-3	D	0.4616
AAK1	-3	E	0.6293
AAK1	-3	F	0.8278
AAK1	-3	G	0.87
AAK1	-3	H	0.6756
AAK1	-3	I	1.0183
AAK1	-3	K	0.8788
AAK1	-3	L	0.8963
AAK1	-3	M	1.3291
AAK1	-3	N	0.9892
AAK1	-3	P	1.8434
AAK1	-3	Q	1.1817
AAK1	-3	R	0.9871
AAK1	-3	S	0.8963
AAK1	-3	T	0.8963
AAK1	-3	V	1.0893
AAK1	-3	W	0.7219
AAK1	-3	Y	0.699
AAK1	-3	s	0.5145
AAK1	-3	t	0.5145
AAK1	-3	y	0.4657
AAK1	-4	A	1.2007
AAK1	-4	C	0.9518
AAK1	-4	D	0.5645
AAK1	-4	E	0.9512
AAK1	-4	F	0.8834
AAK1	-4	G	1.0916
AAK1	-4	H	0.8741
AAK1	-4	I	1.0562
AAK1	-4	K	1.3751
AAK1	-4	L	1.2609
AAK1	-4	M	1.1786
AAK1	-4	N	0.7297
AAK1	-4	P	0.9078
AAK1	-4	Q	1.0655
AAK1	-4	R	1.2151
AAK1	-4	S	1.0516
AAK1	-4	T	1.0516
AAK1	-4	V	1.0516
AAK1	-4	W	0.6859
AAK1	-4	Y	0.9083
AAK1	-4	s	0.5639
AAK1	-4	t	0.5639
AAK1	-4	y	0.5761
AAK1	-5	A	0.4825
AAK1	-5	C	0.7757
AAK1	-5	D	0.2725
AAK1	-5	E	0.2601
AAK1	-5	F	0.7224
AAK1	-5	G	0.4165
AAK1	-5	H	0.5635
AAK1	-5	I	2.6421
AAK1	-5	K	0.4449
AAK1	-5	L	1.6881
AAK1	-5	M	1.4685
AAK1	-5	N	0.4677
AAK1	-5	P	1.2242
AAK1	-5	Q	0.9524
AAK1	-5	R	1.6248
AAK1	-5	S	0.7224
AAK1	-5	T	0.7224
AAK1	-5	V	1.6169
AAK1	-5	W	0.535
AAK1	-5	Y	1.6179
AAK1	-5	s	0.3422
AAK1	-5	t	0.3422
AAK1	-5	y	1.0388
AAK1	0	s	0.1013
AAK1	0	t	1.0
AAK1	0	y	0.0
AAK1	1	A	0.5297
AAK1	1	C	0.379
AAK1	1	D	0.1745
AAK1	1	E	0.1668
AAK1	1	F	0.2188
AAK1	1	G	12.267
AAK1	1	H	0.2721
AAK1	1	I	0.1733
AAK1	1	K	0.3299
AAK1	1	L	0.2223
AAK1	1	M	0.21
AAK1	1	N	0.2894
AAK1	1	P	0.7872
AAK1	1	Q	0.2346
AAK1	1	R	0.4235
AAK1	1	S	0.2429
AAK1	1	T	0.2429
AAK1	1	V	0.1821
AAK1	1	W	0.2757
AAK1	1	Y	0.2429
AAK1	1	s	0.2094
AAK1	1	t	0.2094
AAK1	1	y	0.1705
AAK1	2	A	1.2744
AAK1	2	C	1.1636
AAK1	2	D	0.9386
AAK1	2	E	0.89
AAK1	2	F	0.8257
AAK1	2	G	1.2734
AAK1	2	H	1.0532
AAK1	2	I	0.6291
AAK1	2	K	0.9804
AAK1	2	L	0.7032
AAK1	2	M	0.8174
AAK1	2	N	1.4172
AAK1	2	P	0.8956
AAK1	2	Q	1.3109
AAK1	2	R	1.3158
AAK1	2	S	0.9386
AAK1	2	T	0.9386
AAK1	2	V	0.9222
AAK1	2	W	0.9728
AAK1	2	Y	0.7801
AAK1	2	s	0.6957
AAK1	2	t	0.6957
AAK1	2	y	0.6972
AAK1	3	A	0.9894
AAK1	3	C	1.2613
AAK1	3	D	0.6882
AAK1	3	E	0.6315
AAK1	3	F	0.7306
AAK1	3	G	1.5691
AAK1	3	H	1.1451
AAK1	3	I	0.6591
AAK1	3	K	1.2557
AAK1	3	L	0.8306
AAK1	3	M	0.7427
AAK1	3	N	1.2498
AAK1	3	P	1.1582
AAK1	3	Q	1.059
AAK1	3	R	1.5315
AAK1	3	S	0.9894
AAK1	3	T	0.9894
AAK1	3	V	1.037
AAK1	3	W	0.8171
AAK1	3	Y	0.9053
AAK1	3	s	0.5695
AAK1	3	t	0.5695
AAK1	3	y	0.6096
AAK1	4	A	1.0974
AAK1	4	C	1.0252
AAK1	4	D	0.6619
AAK1	4	E	0.7777
AAK1	4	F	0.8896
AAK1	4	G	1.1928
AAK1	4	H	0.9513
AAK1	4	I	0.7061
AAK1	4	K	1.4123
AAK1	4	L	0.7835
AAK1	4	M	0.7889
AAK1	4	N	1.0056
AAK1	4	P	1.0673
AAK1	4	Q	1.0798
AAK1	4	R	1.5775
AAK1	4	S	0.9513
AAK1	4	T	0.9513
AAK1	4	V	0.7181
AAK1	4	W	1.4046
AAK1	4	Y	0.8857
AAK1	4	s	0.4274
AAK1	4	t	0.4274
AAK1	4	y	0.4597
ACVR2A	-1	A	0.8222
ACVR2A	-1	C	1.7711
ACVR2A	-1	D	1.2305
ACVR2A	-1	E	1.2939
ACVR2A	-1	F	1.1695
ACVR2A	-1	G	0.5701
ACVR2A	-1	H	1.0702
ACVR2A	-1	I	0.8043
ACVR2A	-1	K	0.5321
ACVR2A	-1	L	1.2744
ACVR2A	-1	M	1.44
ACVR2A	-1	N	0.8163
ACVR2A	-1	P	0.6759
ACVR2A	-1	Q	0.8997
ACVR2A	-1	R	0.5773
ACVR2A	-1	S	1.0702
ACVR2A	-1	T	1.0702
ACVR2A	-1	V	1.1047
ACVR2A	-1	W	1.3875
ACVR2A	-1	Y	1.3314
ACVR2A	-1	s	1.3725
ACVR2A	-1	t	1.3725
ACVR2A	-1	y	1.8113
ACVR2A	-2	A	0.5303
ACVR2A	-2	C	1.234
ACVR2A	-2	D	4.3007
ACVR2A	-2	E	5.2302
ACVR2A	-2	F	0.5044
ACVR2A	-2	G	0.5
ACVR2A	-2	H	0.4874
ACVR2A	-2	I	0.4537
ACVR2A	-2	K	0.3544
ACVR2A	-2	L	0.4308
ACVR2A	-2	M	0.3986
ACVR2A	-2	N	0.6759
ACVR2A	-2	P	0.321
ACVR2A	-2	Q	0.7687
ACVR2A	-2	R	0.5045
ACVR2A	-2	S	0.5022
ACVR2A	-2	T	0.5022
ACVR2A	-2	V	0.5022
ACVR2A	-2	W	0.5397
ACVR2A	-2	Y	0.4976
ACVR2A	-2	s	1.0109
ACVR2A	-2	t	1.0109
ACVR2A	-2	y	0.7207
ACVR2A	-3	A	0.9123
ACVR2A	-3	C	1.0
ACVR2A	-3	D	1.7852
ACVR2A	-3	E	2.0478
ACVR2A	-3	F	0.8767
ACVR2A	-3	G	1.11
ACVR2A	-3	H	0.9
ACVR2A	-3	I	0.8548
ACVR2A	-3	K	0.6027
ACVR2A	-3	L	0.8141
ACVR2A	-3	M	0.8754
ACVR2A	-3	N	0.8538
ACVR2A	-3	P	0.8538
ACVR2A	-3	Q	0.8634
ACVR2A	-3	R	0.783
ACVR2A	-3	S	0.8767
ACVR2A	-3	T	0.8767
ACVR2A	-3	V	0.9278
ACVR2A	-3	W	1.0041
ACVR2A	-3	Y	0.9349
ACVR2A	-3	s	2.0483
ACVR2A	-3	t	2.0483
ACVR2A	-3	y	1.6573
ACVR2A	-4	A	0.943
ACVR2A	-4	C	0.9373
ACVR2A	-4	D	1.768
ACVR2A	-4	E	1.4772
ACVR2A	-4	F	1.0519
ACVR2A	-4	G	0.9328
ACVR2A	-4	H	0.9371
ACVR2A	-4	I	0.8946
ACVR2A	-4	K	0.69
ACVR2A	-4	L	0.8839
ACVR2A	-4	M	0.9799
ACVR2A	-4	N	0.8519
ACVR2A	-4	P	0.7915
ACVR2A	-4	Q	0.8197
ACVR2A	-4	R	0.7771
ACVR2A	-4	S	0.9328
ACVR2A	-4	T	0.9328
ACVR2A	-4	V	0.9239
ACVR2A	-4	W	1.2818
ACVR2A	-4	Y	0.9957
ACVR2A	-4	s	1.3745
ACVR2A	-4	t	1.3745
ACVR2A	-4	y	1.1583
ACVR2A	-5	A	0.992
ACVR2A	-5	C	0.8313
ACVR2A	-5	D	1.5097
ACVR2A	-5	E	1.3419
ACVR2A	-5	F	1.0206
ACVR2A	-5	G	0.8178
ACVR2A	-5	H	0.9693
ACVR2A	-5	I	1.0624
ACVR2A	-5	K	0.8675
ACVR2A	-5	L	1.0124
ACVR2A	-5	M	0.886
ACVR2A	-5	N	0.9115
ACVR2A	-5	P	0.7057
ACVR2A	-5	Q	0.7316
ACVR2A	-5	R	0.8081
ACVR2A	-5	S	0.982
ACVR2A	-5	T	0.982
ACVR2A	-5	V	1.0164
ACVR2A	-5	W	1.3652
ACVR2A	-5	Y	0.982
ACVR2A	-5	s	1.3316
ACVR2A	-5	t	1.3316
ACVR2A	-5	y	1.2928
ACVR2A	0	s	0.9833
ACVR2A	0	t	1.0
ACVR2A	0	y	0.0
ACVR2A	1	A	0.6471
ACVR2A	1	C	1.2723
ACVR2A	1	D	1.6836
ACVR2A	1	E	3.2095
ACVR2A	1	F	0.9316
ACVR2A	1	G	0.5753
ACVR2A	1	H	0.7001
ACVR2A	1	I	1.1755
ACVR2A	1	K	0.327
ACVR2A	1	L	0.8732
ACVR2A	1	M	0.9958
ACVR2A	1	N	0.6383
ACVR2A	1	P	0.3404
ACVR2A	1	Q	1.2292
ACVR2A	1	R	0.4162
ACVR2A	1	S	0.9316
ACVR2A	1	T	0.9316
ACVR2A	1	V	1.325
ACVR2A	1	W	0.9614
ACVR2A	1	Y	0.9708
ACVR2A	1	s	4.6492
ACVR2A	1	t	4.6492
ACVR2A	1	y	2.8075
ACVR2A	2	A	1.1214
ACVR2A	2	C	1.0194
ACVR2A	2	D	1.2375
ACVR2A	2	E	1.0262
ACVR2A	2	F	0.9832
ACVR2A	2	G	1.3578
ACVR2A	2	H	1.1221
ACVR2A	2	I	0.9147
ACVR2A	2	K	0.8046
ACVR2A	2	L	0.6742
ACVR2A	2	M	0.8788
ACVR2A	2	N	0.8627
ACVR2A	2	P	0.9887
ACVR2A	2	Q	1.0049
ACVR2A	2	R	0.9352
ACVR2A	2	S	0.9887
ACVR2A	2	T	0.9887
ACVR2A	2	V	0.9992
ACVR2A	2	W	0.9827
ACVR2A	2	Y	1.1059
ACVR2A	2	s	1.1044
ACVR2A	2	t	1.1044
ACVR2A	2	y	0.8772
ACVR2A	3	A	0.8995
ACVR2A	3	C	0.8179
ACVR2A	3	D	1.0568
ACVR2A	3	E	1.4797
ACVR2A	3	F	1.0949
ACVR2A	3	G	1.0064
ACVR2A	3	H	1.0024
ACVR2A	3	I	0.9422
ACVR2A	3	K	0.7851
ACVR2A	3	L	1.0408
ACVR2A	3	M	1.0014
ACVR2A	3	N	0.8913
ACVR2A	3	P	1.1382
ACVR2A	3	Q	0.8481
ACVR2A	3	R	0.6576
ACVR2A	3	S	1.0024
ACVR2A	3	T	1.0024
ACVR2A	3	V	0.9631
ACVR2A	3	W	1.08
ACVR2A	3	Y	1.1126
ACVR2A	3	s	0.8818
ACVR2A	3	t	0.8818
ACVR2A	3	y	1.3859
ACVR2A	4	A	0.8471
ACVR2A	4	C	0.8793
ACVR2A	4	D	1.0875
ACVR2A	4	E	1.0888
ACVR2A	4	F	0.8899
ACVR2A	4	G	0.9241
ACVR2A	4	H	0.9741
ACVR2A	4	I	0.9572
ACVR2A	4	K	0.8961
ACVR2A	4	L	0.8704
ACVR2A	4	M	1.1252
ACVR2A	4	N	0.9455
ACVR2A	4	P	1.2892
ACVR2A	4	Q	1.0467
ACVR2A	4	R	0.8352
ACVR2A	4	S	0.9572
ACVR2A	4	T	0.9572
ACVR2A	4	V	0.8773
ACVR2A	4	W	1.3606
ACVR2A	4	Y	0.985
ACVR2A	4	s	1.1958
ACVR2A	4	t	1.1958
ACVR2A	4	y	1.0015
"""

mock_quantiles_input_file = """Score	AAK1	ACVR2A
-inf	0.0	0.0
-12.0	0.0007	0.0
-11.99	0.0007	0.0
-11.98	0.0007	0.0
-11.97	0.0007	0.0
-11.96	0.0007	0.0
-11.95	0.0007	0.0
-11.94	0.0007	0.0
-11.93	0.0008	0.0
-11.92	0.0008	0.0
-11.91	0.0008	0.0
-11.9	0.0008	0.0
-11.89	0.0009	0.0
-11.88	0.0009	0.0
-11.87	0.0009	0.0
-11.86	0.0009	0.0
-11.85	0.0009	0.0
-11.84	0.001	0.0
-11.83	0.001	0.0
-11.82	0.001	0.0
-11.81	0.0011	0.0
-11.8	0.0011	0.0
-11.79	0.0011	0.0
-11.78	0.0011	0.0
-11.77	0.0012	0.0
-11.76	0.0012	0.0
-11.75	0.0012	0.0
-11.74	0.0012	0.0
-11.73	0.0012	0.0
-11.72	0.0013	0.0
-11.71	0.0013	0.0
-11.7	0.0013	0.0
-11.69	0.0013	0.0
-11.68	0.0013	0.0
-11.67	0.0014	0.0
-11.66	0.0014	0.0
-11.65	0.0014	0.0
-11.64	0.0015	0.0
-11.63	0.0015	0.0
-11.62	0.0015	0.0
-11.61	0.0015	0.0
-11.6	0.0015	0.0
-11.59	0.0015	0.0
-11.58	0.0016	0.0
-11.57	0.0016	0.0
-11.56	0.0016	0.0
-11.55	0.0017	0.0
-11.54	0.0017	0.0
-11.53	0.0017	0.0
-11.52	0.0017	0.0
-11.51	0.0018	0.0
-11.5	0.0018	0.0
-11.49	0.0018	0.0
-11.48	0.0019	0.0
-11.47	0.0019	0.0
-11.46	0.0019	0.0
-11.45	0.002	0.0
-11.44	0.002	0.0
-11.43	0.002	0.0
-11.42	0.002	0.0
-11.41	0.002	0.0
-11.4	0.0021	0.0
-11.39	0.0021	0.0
-11.38	0.0021	0.0
-11.37	0.0022	0.0
-11.36	0.0022	0.0
-11.35	0.0023	0.0
-11.34	0.0023	0.0
-11.33	0.0023	0.0
-11.32	0.0024	0.0
-11.31	0.0024	0.0
-11.3	0.0024	0.0
-11.29	0.0025	0.0
-11.28	0.0026	0.0
-11.27	0.0026	0.0
-11.26	0.0027	0.0
-11.25	0.0027	0.0
-11.24	0.0027	0.0
-11.23	0.0027	0.0
-11.22	0.0028	0.0
-11.21	0.0029	0.0
-11.2	0.0029	0.0
-11.19	0.003	0.0
-11.18	0.003	0.0
-11.17	0.0031	0.0
-11.16	0.0031	0.0
-11.15	0.0032	0.0
-11.14	0.0032	0.0
-11.13	0.0032	0.0
-11.12	0.0033	0.0
-11.11	0.0033	0.0
-11.1	0.0034	0.0
-11.09	0.0034	0.0
-11.08	0.0035	0.0
-11.07	0.0035	0.0
-11.06	0.0036	0.0
-11.05	0.0036	0.0
-11.04	0.0037	0.0
-11.03	0.0037	0.0
-11.02	0.0038	0.0
-11.01	0.0038	0.0
-11.0	0.0039	0.0
-10.99	0.0039	0.0
-10.98	0.004	0.0
-10.97	0.0041	0.0
-10.96	0.0041	0.0
-10.95	0.0042	0.0
-10.94	0.0042	0.0
-10.93	0.0043	0.0
-10.92	0.0043	0.0
-10.91	0.0044	0.0
-10.9	0.0044	0.0
-10.89	0.0045	0.0
-10.879999999999999	0.0046	0.0
-10.87	0.0047	0.0
-10.86	0.0047	0.0
-10.85	0.0048	0.0
-10.84	0.0049	0.0
-10.83	0.0049	0.0
-10.82	0.005	0.0
-10.81	0.005	0.0
-10.8	0.0051	0.0
-10.79	0.0052	0.0
-10.78	0.0053	0.0
-10.77	0.0054	0.0
-10.76	0.0055	0.0
-10.75	0.0056	0.0
-10.74	0.0057	0.0
-10.73	0.0058	0.0
-10.72	0.0059	0.0
-10.71	0.006	0.0
-10.7	0.006	0.0
-10.69	0.0062	0.0
-10.68	0.0063	0.0
-10.67	0.0063	0.0
-10.66	0.0064	0.0
-10.65	0.0065	0.0
-10.64	0.0066	0.0
-10.629999999999999	0.0067	0.0
-10.62	0.0068	0.0
-10.61	0.0069	0.0
-10.6	0.007	0.0
-10.59	0.0071	0.0
-10.58	0.0072	0.0
-10.57	0.0073	0.0
-10.56	0.0074	0.0
-10.55	0.0075	0.0
-10.54	0.0076	0.0
-10.53	0.0077	0.0
-10.52	0.0078	0.0
-10.51	0.0079	0.0
-10.5	0.008	0.0
-10.49	0.0081	0.0
-10.48	0.0082	0.0
-10.47	0.0084	0.0
-10.46	0.0085	0.0
-10.45	0.0086	0.0
-10.44	0.0087	0.0
-10.43	0.0088	0.0
-10.42	0.0089	0.0
-10.41	0.009	0.0
-10.4	0.0092	0.0
-10.39	0.0092	0.0
-10.379999999999999	0.0093	0.0
-10.37	0.0094	0.0
-10.36	0.0095	0.0
-10.35	0.0096	0.0
-10.34	0.0098	0.0
-10.33	0.0098	0.0
-10.32	0.0099	0.0
-10.31	0.01	0.0
-10.3	0.0101	0.0
-10.29	0.0102	0.0
-10.28	0.0103	0.0
-10.27	0.0104	0.0
-10.26	0.0105	0.0
-10.25	0.0106	0.0
-10.24	0.0107	0.0
-10.23	0.0108	0.0
-10.22	0.0108	0.0
-10.21	0.011	0.0
-10.2	0.0111	0.0
-10.19	0.0113	0.0
-10.18	0.0114	0.0
-10.17	0.0116	0.0
-10.16	0.0117	0.0
-10.15	0.0118	0.0
-10.14	0.012	0.0
-10.129999999999999	0.0122	0.0
-10.12	0.0123	0.0
-10.11	0.0125	0.0
-10.1	0.0127	0.0
-10.09	0.0128	0.0
-10.08	0.0129	0.0
-10.07	0.0131	0.0
-10.06	0.0133	0.0
-10.05	0.0134	0.0
-10.04	0.0136	0.0
-10.03	0.0139	0.0
-10.02	0.014	0.0
-10.01	0.0141	0.0
-10.0	0.0143	0.0
-9.99	0.0144	0.0
-9.98	0.0145	0.0
-9.969999999999999	0.0147	0.0
-9.96	0.0149	0.0
-9.95	0.0151	0.0
-9.94	0.0152	0.0
-9.93	0.0154	0.0
-9.92	0.0155	0.0
-9.91	0.0156	0.0
-9.9	0.0158	0.0
-9.89	0.0161	0.0
-9.879999999999999	0.0163	0.0
-9.870000000000001	0.0164	0.0
-9.86	0.0166	0.0
-9.85	0.0167	0.0
-9.84	0.0168	0.0
-9.83	0.017	0.0
-9.82	0.0172	0.0
-9.81	0.0174	0.0
-9.8	0.0175	0.0
-9.79	0.0176	0.0
-9.78	0.0178	0.0
-9.77	0.0179	0.0
-9.76	0.0181	0.0
-9.75	0.0183	0.0
-9.74	0.0185	0.0
-9.73	0.0187	0.0
-9.719999999999999	0.0189	0.0
-9.71	0.0191	0.0
-9.7	0.0193	0.0
-9.69	0.0195	0.0
-9.68	0.0198	0.0
-9.67	0.02	0.0
-9.66	0.0202	0.0
-9.65	0.0204	0.0
-9.64	0.0205	0.0
-9.629999999999999	0.0208	0.0
-9.620000000000001	0.0209	0.0
-9.61	0.0211	0.0
-9.6	0.0214	0.0
-9.59	0.0216	0.0
-9.58	0.0217	0.0
-9.57	0.022	0.0
-9.56	0.0222	0.0
-9.55	0.0225	0.0
-9.54	0.0227	0.0
-9.53	0.0229	0.0
-9.52	0.0232	0.0
-9.51	0.0233	0.0
-9.5	0.0235	0.0
-9.49	0.0237	0.0
-9.48	0.0239	0.0
-9.469999999999999	0.0241	0.0
-9.46	0.0243	0.0
-9.45	0.0246	0.0
-9.44	0.0248	0.0
-9.43	0.025	0.0
-9.42	0.0253	0.0
-9.41	0.0255	0.0
-9.4	0.0259	0.0
-9.39	0.0261	0.0
-9.379999999999999	0.0263	0.0
-9.370000000000001	0.0266	0.0
-9.36	0.0271	0.0
-9.35	0.0274	0.0
-9.34	0.0276	0.0
-9.33	0.0279	0.0
-9.32	0.0281	0.0
-9.31	0.0284	0.0
-9.3	0.0287	0.0
-9.29	0.029	0.0
-9.28	0.0294	0.0
-9.27	0.0297	0.0
-9.26	0.03	0.0
-9.25	0.0304	0.0
-9.24	0.0308	0.0
-9.23	0.0311	0.0
-9.219999999999999	0.0313	0.0
-9.21	0.0317	0.0
-9.2	0.032	0.0
-9.19	0.0322	0.0
-9.18	0.0324	0.0
-9.17	0.0327	0.0
-9.16	0.0329	0.0
-9.15	0.0332	0.0
-9.14	0.0335	0.0
-9.129999999999999	0.0339	0.0
-9.120000000000001	0.0343	0.0
-9.11	0.0346	0.0
-9.1	0.0349	0.0
-9.09	0.0352	0.0
-9.08	0.0356	0.0
-9.07	0.0359	0.0
-9.06	0.0363	0.0
-9.05	0.0367	0.0
-9.04	0.0371	0.0
-9.03	0.0374	0.0
-9.02	0.0378	0.0
-9.01	0.0382	0.0
-9.0	0.0385	0.0
-8.99	0.0388	0.0
-8.98	0.0392	0.0
-8.969999999999999	0.0395	0.0
-8.96	0.0398	0.0
-8.95	0.0402	0.0
-8.94	0.0404	0.0
-8.93	0.0408	0.0
-8.92	0.0411	0.0
-8.91	0.0415	0.0
-8.9	0.042	0.0
-8.89	0.0423	0.0
-8.879999999999999	0.0427	0.0
-8.870000000000001	0.0431	0.0
-8.86	0.0435	0.0
-8.85	0.0439	0.0
-8.84	0.0443	0.0
-8.83	0.0447	0.0
-8.82	0.0451	0.0
-8.81	0.0456	0.0
-8.8	0.0459	0.0
-8.79	0.0464	0.0
-8.78	0.0469	0.0
-8.77	0.0473	0.0
-8.76	0.0477	0.0
-8.75	0.0481	0.0
-8.74	0.0484	0.0
-8.73	0.0488	0.0
-8.719999999999999	0.0491	0.0
-8.71	0.0496	0.0
-8.7	0.0501	0.0
-8.69	0.0504	0.0
-8.68	0.051	0.0
-8.67	0.0516	0.0
-8.66	0.052	0.0
-8.65	0.0524	0.0
-8.64	0.0528	0.0
-8.629999999999999	0.0533	0.0
-8.620000000000001	0.0537	0.0
-8.61	0.0542	0.0
-8.6	0.0548	0.0
-8.59	0.0551	0.0
-8.58	0.0556	0.0
-8.57	0.0561	0.0
-8.56	0.0566	0.0
-8.55	0.057	0.0
-8.54	0.0575	0.0
-8.53	0.058	0.0
-8.52	0.0585	0.0
-8.51	0.0589	0.0
-8.5	0.0594	0.0
-8.49	0.06	0.0
-8.48	0.0605	0.0
-8.469999999999999	0.0609	0.0
-8.46	0.0614	0.0
-8.45	0.0619	0.0
-8.44	0.0624	0.0
-8.43	0.063	0.0
-8.42	0.0636	0.0
-8.41	0.0642	0.0
-8.4	0.0647	0.0
-8.39	0.0652	0.0
-8.379999999999999	0.0656	0.0
-8.370000000000001	0.0661	0.0
-8.36	0.0667	0.0
-8.35	0.0672	0.0
-8.34	0.0677	0.0
-8.33	0.0682	0.0
-8.32	0.0687	0.0
-8.31	0.0692	0.0
-8.3	0.0696	0.0
-8.29	0.0702	0.0
-8.28	0.0707	0.0
-8.27	0.0713	0.0
-8.26	0.0719	0.0
-8.25	0.0724	0.0
-8.24	0.0731	0.0
-8.23	0.0738	0.0
-8.219999999999999	0.0744	0.0
-8.21	0.0751	0.0
-8.2	0.0756	0.0
-8.19	0.0762	0.0
-8.18	0.0769	0.0
-8.17	0.0776	0.0
-8.16	0.0782	0.0
-8.15	0.0788	0.0
-8.14	0.0795	0.0
-8.129999999999999	0.0801	0.0
-8.120000000000001	0.0807	0.0
-8.11	0.0812	0.0
-8.1	0.082	0.0
-8.09	0.0825	0.0
-8.08	0.0831	0.0
-8.07	0.0836	0.0
-8.06	0.0845	0.0
-8.05	0.0851	0.0
-8.04	0.0858	0.0
-8.03	0.0866	0.0
-8.02	0.0872	0.0
-8.01	0.0878	0.0
-8.0	0.0883	0.0
-7.99	0.0891	0.0
-7.9799999999999995	0.0899	0.0
-7.97	0.0908	0.0
-7.96	0.0915	0.0
-7.95	0.0921	0.0
-7.9399999999999995	0.0928	0.0
-7.93	0.0934	0.0
-7.92	0.0941	0.0
-7.91	0.0949	0.0
-7.9	0.0956	0.0
-7.89	0.0965	0.0
-7.88	0.0974	0.0
-7.87	0.0979	0.0
-7.86	0.0985	0.0
-7.85	0.0994	0.0
-7.84	0.1001	0.0
-7.83	0.1008	0.0
-7.82	0.1016	0.0
-7.81	0.1023	0.0
-7.8	0.1031	0.0
-7.79	0.1038	0.0
-7.78	0.1047	0.0
-7.77	0.1056	0.0
-7.76	0.1064	0.0
-7.75	0.1072	0.0
-7.74	0.108	0.0
-7.7299999999999995	0.1087	0.0
-7.72	0.1096	0.0
-7.71	0.1104	0.0
-7.7	0.1112	0.0
-7.6899999999999995	0.1119	0.0
-7.68	0.1126	0.0
-7.67	0.1135	0.0
-7.66	0.1144	0.0
-7.6499999999999995	0.1152	0.0
-7.64	0.1162	0.0
-7.63	0.117	0.0
-7.62	0.1179	0.0
-7.61	0.1187	0.0
-7.6	0.1195	0.0
-7.59	0.1206	0.0
-7.58	0.1215	0.0
-7.57	0.1224	0.0
-7.56	0.1232	0.0
-7.55	0.1241	0.0
-7.54	0.125	0.0
-7.53	0.126	0.0
-7.52	0.1269	0.0
-7.51	0.1278	0.0
-7.5	0.1287	0.0
-7.49	0.1295	0.0
-7.4799999999999995	0.1307	0.0
-7.47	0.1317	0.0
-7.46	0.1325	0.0
-7.45	0.1335	0.0
-7.4399999999999995	0.1345	0.0
-7.43	0.1353	0.0
-7.42	0.1361	0.0
-7.41	0.1368	0.0
-7.3999999999999995	0.1378	0.0
-7.39	0.1388	0.0
-7.38	0.1398	0.0
-7.37	0.1406	0.0
-7.36	0.1415	0.0
-7.35	0.1425	0.0
-7.34	0.1437	0.0
-7.33	0.1447	0.0
-7.32	0.1457	0.0
-7.31	0.1469	0.0
-7.3	0.1477	0.0
-7.29	0.1488	0.0
-7.28	0.15	0.0
-7.27	0.1511	0.0
-7.26	0.152	0.0
-7.25	0.1529	0.0
-7.24	0.154	0.0
-7.2299999999999995	0.155	0.0
-7.22	0.156	0.0
-7.21	0.157	0.0
-7.2	0.1581	0.0
-7.1899999999999995	0.1592	0.0
-7.18	0.1605	0.0
-7.17	0.1614	0.0
-7.16	0.1623	0.0
-7.1499999999999995	0.1635	0.0
-7.14	0.1645	0.0
-7.13	0.1657	0.0
-7.12	0.1666	0.0
-7.11	0.1675	0.0
-7.1	0.1687	0.0
-7.09	0.17	0.0
-7.08	0.171	0.0
-7.07	0.1721	0.0
-7.06	0.1734	0.0
-7.05	0.1743	0.0
-7.04	0.1755	0.0
-7.03	0.1768	0.0
-7.02	0.1779	0.0
-7.01	0.1791	0.0
-7.0	0.1801	0.0
-6.99	0.1812	0.0
-6.9799999999999995	0.1823	0.0
-6.97	0.1832	0.0
-6.96	0.1843	0.0
-6.95	0.1853	0.0
-6.9399999999999995	0.1865	0.0
-6.93	0.1877	0.0
-6.92	0.189	0.0
-6.91	0.1901	0.0
-6.8999999999999995	0.1913	0.0
-6.89	0.1925	0.0
-6.88	0.1938	0.0
-6.87	0.195	0.0
-6.86	0.1962	0.0
-6.85	0.1973	0.0
-6.84	0.1987	0.0
-6.83	0.1999	0.0
-6.82	0.2008	0.0
-6.81	0.2022	0.0
-6.8	0.2032	0.0
-6.79	0.2046	0.0
-6.78	0.2058	0.0
-6.77	0.207	0.0
-6.76	0.2082	0.0
-6.75	0.2095	0.0
-6.74	0.2108	0.0
-6.7299999999999995	0.2121	0.0
-6.72	0.2133	0.0
-6.71	0.2144	0.0
-6.7	0.2156	0.0
-6.6899999999999995	0.2171	0.0
-6.68	0.2182	0.0
-6.67	0.2193	0.0
-6.66	0.2204	0.0
-6.6499999999999995	0.2215	0.0
-6.64	0.2231	0.0
-6.63	0.2244	0.0
-6.62	0.2256	0.0
-6.61	0.227	0.0
-6.6	0.2283	0.0
-6.59	0.2294	0.0
-6.58	0.2304	0.0
-6.57	0.2318	0.0
-6.56	0.2334	0.0
-6.55	0.2346	0.0
-6.54	0.2359	0.0
-6.53	0.2371	0.0
-6.52	0.2386	0.0
-6.51	0.2398	0.0
-6.5	0.241	0.0
-6.49	0.2422	0.0
-6.4799999999999995	0.2435	0.0
-6.47	0.2449	0.0
-6.46	0.2461	0.0
-6.45	0.2474	0.0
-6.4399999999999995	0.2488	0.0
-6.43	0.2502	0.0
-6.42	0.2515	0.0
-6.41	0.2527	0.0
-6.3999999999999995	0.2539	0.0
-6.39	0.2552	0.0
-6.38	0.2565	0.0
-6.37	0.2579	0.0
-6.36	0.2595	0.0
-6.35	0.2606	0.0
-6.34	0.2618	0.0
-6.33	0.2631	0.0
-6.32	0.2646	0.0
-6.31	0.2659	0.0
-6.3	0.2674	0.0
-6.29	0.2689	0.0
-6.28	0.2704	0.0
-6.27	0.2725	0.0
-6.26	0.2739	0.0
-6.25	0.2752	0.0
-6.24	0.2765	0.0
-6.2299999999999995	0.2779	0.0
-6.22	0.2793	0.0
-6.21	0.2805	0.0
-6.2	0.2819	0.0
-6.1899999999999995	0.2832	0.0
-6.18	0.2846	0.0
-6.17	0.2862	0.0
-6.16	0.2879	0.0
-6.1499999999999995	0.2895	0.0
-6.14	0.291	0.0
-6.13	0.2925	0.0
-6.12	0.2939	0.0
-6.11	0.2954	0.0
-6.1	0.2968	0.0
-6.09	0.2983	0.0
-6.08	0.2996	0.0
-6.07	0.3008	0.0
-6.06	0.3021	0.0
-6.05	0.3038	0.0001
-6.04	0.3053	0.0001
-6.03	0.3068	0.0001
-6.02	0.3082	0.0001
-6.01	0.3098	0.0001
-6.0	0.3113	0.0001
-5.99	0.3128	0.0001
-5.9799999999999995	0.3142	0.0001
-5.97	0.3156	0.0001
-5.96	0.3172	0.0001
-5.95	0.3185	0.0001
-5.9399999999999995	0.3199	0.0001
-5.93	0.3216	0.0001
-5.92	0.3233	0.0001
-5.91	0.3249	0.0001
-5.8999999999999995	0.3265	0.0001
-5.89	0.3278	0.0001
-5.88	0.3295	0.0001
-5.87	0.3312	0.0001
-5.859999999999999	0.3327	0.0001
-5.85	0.3341	0.0001
-5.84	0.3357	0.0002
-5.83	0.3372	0.0002
-5.82	0.3387	0.0002
-5.81	0.34	0.0002
-5.8	0.3416	0.0002
-5.79	0.3432	0.0002
-5.78	0.3445	0.0002
-5.77	0.3459	0.0002
-5.76	0.3474	0.0002
-5.75	0.3491	0.0002
-5.74	0.3507	0.0003
-5.7299999999999995	0.3524	0.0003
-5.72	0.3542	0.0003
-5.71	0.3558	0.0003
-5.7	0.3572	0.0004
-5.6899999999999995	0.3585	0.0004
-5.68	0.3599	0.0004
-5.67	0.3615	0.0004
-5.66	0.363	0.0004
-5.6499999999999995	0.3646	0.0004
-5.64	0.3663	0.0004
-5.63	0.3679	0.0004
-5.62	0.3692	0.0004
-5.609999999999999	0.3707	0.0005
-5.6	0.3725	0.0005
-5.59	0.3739	0.0005
-5.58	0.3756	0.0005
-5.57	0.3772	0.0006
-5.56	0.3785	0.0006
-5.55	0.3803	0.0006
-5.54	0.3819	0.0007
-5.53	0.3834	0.0007
-5.52	0.3851	0.0007
-5.51	0.3868	0.0008
-5.5	0.3886	0.0008
-5.49	0.3901	0.0008
-5.4799999999999995	0.3914	0.0009
-5.47	0.393	0.0009
-5.46	0.3943	0.001
-5.45	0.3961	0.001
-5.4399999999999995	0.3976	0.0011
-5.43	0.3995	0.0011
-5.42	0.4013	0.0012
-5.41	0.4028	0.0012
-5.3999999999999995	0.4045	0.0012
-5.39	0.4061	0.0012
-5.38	0.4079	0.0013
-5.37	0.4094	0.0014
-5.359999999999999	0.411	0.0014
-5.35	0.4127	0.0016
-5.34	0.4146	0.0016
-5.33	0.4163	0.0017
-5.32	0.4178	0.0017
-5.31	0.4193	0.0018
-5.3	0.4208	0.0019
-5.29	0.4225	0.002
-5.28	0.4241	0.002
-5.27	0.4255	0.0021
-5.26	0.4272	0.0022
-5.25	0.4288	0.0023
-5.24	0.4303	0.0024
-5.2299999999999995	0.432	0.0025
-5.22	0.4335	0.0025
-5.21	0.4351	0.0026
-5.2	0.4369	0.0028
-5.1899999999999995	0.4385	0.0029
-5.18	0.4402	0.0029
-5.17	0.4416	0.003
-5.16	0.4431	0.0031
-5.1499999999999995	0.4447	0.0033
-5.14	0.4463	0.0034
-5.13	0.448	0.0036
-5.12	0.4497	0.0038
-5.109999999999999	0.4512	0.0039
-5.1	0.4528	0.004
-5.09	0.4545	0.0042
-5.08	0.456	0.0044
-5.07	0.4577	0.0045
-5.06	0.4594	0.0046
-5.05	0.4611	0.0047
-5.04	0.4627	0.0048
-5.03	0.4642	0.005
-5.02	0.4657	0.0053
-5.01	0.4674	0.0055
-5.0	0.469	0.0056
-4.99	0.4707	0.0058
-4.9799999999999995	0.4723	0.006
-4.97	0.4738	0.0062
-4.96	0.4757	0.0064
-4.95	0.4775	0.0066
-4.9399999999999995	0.4792	0.0068
-4.93	0.4806	0.007
-4.92	0.4823	0.0072
-4.91	0.4837	0.0074
-4.8999999999999995	0.4852	0.0075
-4.89	0.4869	0.0078
-4.88	0.4883	0.008
-4.87	0.49	0.0083
-4.859999999999999	0.4916	0.0084
-4.85	0.4933	0.0087
-4.84	0.495	0.009
-4.83	0.4963	0.0093
-4.82	0.4981	0.0096
-4.81	0.4995	0.0099
-4.8	0.5011	0.0102
-4.79	0.5024	0.0104
-4.78	0.5041	0.0108
-4.77	0.5059	0.011
-4.76	0.5073	0.0113
-4.75	0.5089	0.0117
-4.74	0.5105	0.012
-4.7299999999999995	0.5122	0.0123
-4.72	0.5133	0.0126
-4.71	0.5151	0.0129
-4.7	0.5165	0.0133
-4.6899999999999995	0.5181	0.0137
-4.68	0.5199	0.014
-4.67	0.5214	0.0144
-4.66	0.5229	0.0149
-4.6499999999999995	0.5244	0.0153
-4.64	0.5259	0.0158
-4.63	0.5274	0.0163
-4.62	0.529	0.0167
-4.609999999999999	0.5302	0.0172
-4.6	0.5317	0.0176
-4.59	0.5331	0.018
-4.58	0.5344	0.0183
-4.57	0.536	0.0187
-4.56	0.5377	0.0192
-4.55	0.5392	0.0196
-4.54	0.5406	0.0201
-4.53	0.542	0.0206
-4.52	0.5434	0.0211
-4.51	0.5451	0.0215
-4.5	0.5465	0.022
-4.49	0.5478	0.0226
-4.4799999999999995	0.5494	0.0232
-4.47	0.5509	0.0238
-4.46	0.5526	0.0243
-4.45	0.554	0.0248
-4.4399999999999995	0.5556	0.0255
-4.43	0.557	0.026
-4.42	0.5586	0.0267
-4.41	0.5602	0.0272
-4.3999999999999995	0.5615	0.0278
-4.39	0.563	0.0284
-4.38	0.5642	0.0291
-4.37	0.5658	0.0298
-4.359999999999999	0.5673	0.0304
-4.35	0.5687	0.0311
-4.34	0.5701	0.0318
-4.33	0.5715	0.0325
-4.32	0.5729	0.0331
-4.31	0.5742	0.0339
-4.3	0.5755	0.0347
-4.29	0.577	0.0356
-4.28	0.5783	0.0363
-4.27	0.5796	0.0371
-4.26	0.5812	0.0379
-4.25	0.5827	0.0386
-4.24	0.5842	0.0393
-4.2299999999999995	0.5856	0.04
-4.22	0.5869	0.0409
-4.21	0.588	0.0418
-4.2	0.5895	0.0426
-4.1899999999999995	0.5911	0.0435
-4.18	0.5924	0.0443
-4.17	0.5939	0.0451
-4.16	0.5952	0.0461
-4.1499999999999995	0.5966	0.0469
-4.14	0.598	0.0478
-4.13	0.5996	0.0488
-4.12	0.6009	0.0497
-4.109999999999999	0.6024	0.0507
-4.1	0.6038	0.0518
-4.09	0.6052	0.0525
-4.08	0.6066	0.0535
-4.069999999999999	0.608	0.0543
-4.06	0.609	0.0554
-4.05	0.6105	0.0562
-4.04	0.6118	0.0571
-4.03	0.6133	0.0582
-4.02	0.6146	0.0594
-4.01	0.6157	0.0605
-4.0	0.6172	0.0617
-3.99	0.6185	0.0628
-3.9800000000000004	0.6197	0.064
-3.9700000000000006	0.6214	0.0651
-3.959999999999999	0.6226	0.0662
-3.9499999999999993	0.6239	0.0672
-3.9399999999999995	0.6252	0.0681
-3.9299999999999997	0.6263	0.0691
-3.92	0.6277	0.0701
-3.91	0.6292	0.0713
-3.9000000000000004	0.6304	0.0726
-3.8900000000000006	0.632	0.0738
-3.879999999999999	0.6333	0.0751
-3.869999999999999	0.6348	0.0766
-3.8599999999999994	0.6361	0.0777
-3.8499999999999996	0.6373	0.0788
-3.84	0.6386	0.0801
-3.83	0.6401	0.0813
-3.8200000000000003	0.6415	0.0827
-3.8100000000000005	0.643	0.0838
-3.8000000000000007	0.6441	0.085
-3.789999999999999	0.6452	0.0862
-3.7799999999999994	0.6464	0.0874
-3.7699999999999996	0.6475	0.0888
-3.76	0.6488	0.0903
-3.75	0.6498	0.0917
-3.74	0.6511	0.0933
-3.7300000000000004	0.6525	0.0947
-3.7200000000000006	0.6538	0.0961
-3.709999999999999	0.6551	0.0977
-3.6999999999999993	0.6562	0.0993
-3.6899999999999995	0.6576	0.1009
-3.6799999999999997	0.6588	0.1023
-3.67	0.66	0.1039
-3.66	0.6615	0.1055
-3.6500000000000004	0.6628	0.1071
-3.6400000000000006	0.6641	0.1085
-3.629999999999999	0.6651	0.1098
-3.619999999999999	0.6664	0.1116
-3.6099999999999994	0.6675	0.1129
-3.5999999999999996	0.6688	0.1145
-3.59	0.6699	0.116
-3.58	0.6713	0.1179
-3.5700000000000003	0.6724	0.1195
-3.5600000000000005	0.6735	0.1214
-3.5500000000000007	0.6746	0.1231
-3.539999999999999	0.6756	0.1248
-3.5299999999999994	0.6772	0.1265
-3.5199999999999996	0.6784	0.1284
-3.51	0.6795	0.1301
-3.5	0.6806	0.1317
-3.49	0.6819	0.1333
-3.4800000000000004	0.6829	0.1349
-3.4700000000000006	0.684	0.1366
-3.459999999999999	0.6854	0.1382
-3.4499999999999993	0.6866	0.1401
-3.4399999999999995	0.6878	0.1418
-3.4299999999999997	0.6888	0.1436
-3.42	0.6898	0.1451
-3.41	0.691	0.147
-3.4000000000000004	0.6923	0.1487
-3.3900000000000006	0.6934	0.1506
-3.379999999999999	0.6945	0.1524
-3.369999999999999	0.6956	0.1543
-3.3599999999999994	0.6966	0.156
-3.3499999999999996	0.6977	0.1579
-3.34	0.6989	0.1597
-3.33	0.7001	0.1615
-3.3200000000000003	0.7011	0.1632
-3.3100000000000005	0.7023	0.1651
-3.299999999999999	0.7034	0.167
-3.289999999999999	0.7046	0.1687
-3.2799999999999994	0.7059	0.1704
-3.2699999999999996	0.7069	0.1724
-3.26	0.7082	0.1743
-3.25	0.7092	0.1761
-3.24	0.7101	0.1779
-3.2300000000000004	0.7112	0.1798
-3.2200000000000006	0.7121	0.1815
-3.209999999999999	0.7131	0.1833
-3.1999999999999993	0.7141	0.1852
-3.1899999999999995	0.7152	0.1871
-3.1799999999999997	0.7164	0.1893
-3.17	0.7175	0.1915
-3.16	0.7188	0.1933
-3.1500000000000004	0.7197	0.1955
-3.1400000000000006	0.7208	0.1976
-3.129999999999999	0.7219	0.1998
-3.119999999999999	0.723	0.2019
-3.1099999999999994	0.724	0.2039
-3.0999999999999996	0.725	0.206
-3.09	0.7262	0.2085
-3.08	0.7272	0.2107
-3.0700000000000003	0.7282	0.2125
-3.0600000000000005	0.7294	0.2147
-3.049999999999999	0.7306	0.2168
-3.039999999999999	0.7315	0.2192
-3.0299999999999994	0.7325	0.2214
-3.0199999999999996	0.7335	0.2234
-3.01	0.7346	0.2255
-3.0	0.7357	0.2276
-2.99	0.7366	0.2297
-2.9800000000000004	0.7377	0.2318
-2.9700000000000006	0.7387	0.2341
-2.959999999999999	0.7396	0.2361
-2.9499999999999993	0.7406	0.2384
-2.9399999999999995	0.7415	0.2403
-2.9299999999999997	0.7426	0.2423
-2.92	0.7435	0.2449
-2.91	0.7446	0.247
-2.9000000000000004	0.7457	0.2494
-2.8900000000000006	0.7466	0.2517
-2.879999999999999	0.7477	0.2539
-2.869999999999999	0.7486	0.2558
-2.8599999999999994	0.7497	0.2582
-2.8499999999999996	0.7507	0.2611
-2.84	0.7518	0.2636
-2.83	0.7528	0.2662
-2.8200000000000003	0.7538	0.2685
-2.8100000000000005	0.7547	0.2708
-2.799999999999999	0.7558	0.2732
-2.789999999999999	0.7569	0.2753
-2.7799999999999994	0.7577	0.2775
-2.7699999999999996	0.7587	0.28
-2.76	0.7596	0.2824
-2.75	0.7604	0.2848
-2.74	0.7614	0.2869
-2.7300000000000004	0.7623	0.2893
-2.7200000000000006	0.7631	0.2918
-2.709999999999999	0.7639	0.2956
-2.6999999999999993	0.7648	0.2978
-2.6899999999999995	0.7655	0.3002
-2.6799999999999997	0.7666	0.3026
-2.67	0.7675	0.305
-2.66	0.7685	0.3075
-2.6500000000000004	0.7694	0.31
-2.6400000000000006	0.7702	0.3122
-2.629999999999999	0.7709	0.3147
-2.619999999999999	0.7718	0.3173
-2.6099999999999994	0.7727	0.3199
-2.5999999999999996	0.7735	0.3225
-2.59	0.7745	0.3256
-2.58	0.7755	0.328
-2.5700000000000003	0.7764	0.3311
-2.5600000000000005	0.7774	0.3337
-2.549999999999999	0.7783	0.3365
-2.539999999999999	0.7789	0.3397
-2.5299999999999994	0.7799	0.342
-2.5199999999999996	0.7809	0.3446
-2.51	0.7817	0.3474
-2.5	0.7825	0.35
-2.49	0.7835	0.3526
-2.4800000000000004	0.7843	0.3553
-2.4700000000000006	0.7851	0.358
-2.459999999999999	0.7859	0.3645
-2.4499999999999993	0.7868	0.3668
-2.4399999999999995	0.7874	0.3695
-2.4299999999999997	0.7883	0.3718
-2.42	0.7891	0.3742
-2.41	0.7901	0.3768
-2.4000000000000004	0.7909	0.3796
-2.3900000000000006	0.7919	0.3823
-2.379999999999999	0.7929	0.3849
-2.369999999999999	0.7937	0.3875
-2.3599999999999994	0.7945	0.3903
-2.3499999999999996	0.7953	0.3932
-2.34	0.7961	0.3955
-2.33	0.7968	0.3982
-2.3200000000000003	0.7976	0.4022
-2.3100000000000005	0.7983	0.4046
-2.299999999999999	0.7991	0.4073
-2.289999999999999	0.7999	0.4105
-2.2799999999999994	0.8008	0.4131
-2.2699999999999996	0.8014	0.4161
-2.26	0.8024	0.4184
-2.25	0.8033	0.4211
-2.24	0.804	0.4239
-2.2300000000000004	0.8048	0.4267
-2.2200000000000006	0.8056	0.4293
-2.209999999999999	0.8062	0.4321
-2.1999999999999993	0.8068	0.4347
-2.1899999999999995	0.8078	0.4372
-2.1799999999999997	0.8086	0.4398
-2.17	0.8094	0.4424
-2.16	0.8104	0.4448
-2.1500000000000004	0.8112	0.4476
-2.1400000000000006	0.8119	0.4503
-2.129999999999999	0.8125	0.453
-2.119999999999999	0.8134	0.4555
-2.1099999999999994	0.8143	0.4579
-2.0999999999999996	0.815	0.4608
-2.09	0.8157	0.4637
-2.08	0.8164	0.4662
-2.0700000000000003	0.8171	0.4688
-2.0600000000000005	0.8179	0.4713
-2.049999999999999	0.8186	0.4737
-2.039999999999999	0.8192	0.4762
-2.0299999999999994	0.8202	0.4789
-2.0199999999999996	0.821	0.4813
-2.01	0.8217	0.4838
-2.0	0.8226	0.4866
-1.9900000000000002	0.8233	0.4894
-1.9800000000000004	0.8241	0.4919
-1.9700000000000006	0.8249	0.4944
-1.959999999999999	0.8256	0.4969
-1.9499999999999993	0.8263	0.4992
-1.9399999999999995	0.8271	0.5017
-1.9299999999999997	0.8279	0.5045
-1.92	0.8288	0.5072
-1.9100000000000001	0.8296	0.51
-1.9000000000000004	0.8302	0.5125
-1.8900000000000006	0.831	0.5151
-1.879999999999999	0.8317	0.5174
-1.8699999999999992	0.8326	0.5201
-1.8599999999999994	0.8335	0.5227
-1.8499999999999996	0.8341	0.5249
-1.8399999999999999	0.8348	0.527
-1.83	0.8356	0.5293
-1.8200000000000003	0.8364	0.5314
-1.8100000000000005	0.8372	0.5338
-1.799999999999999	0.8379	0.5362
-1.7899999999999991	0.8387	0.5387
-1.7799999999999994	0.8393	0.5407
-1.7699999999999996	0.84	0.5433
-1.7599999999999998	0.8407	0.5458
-1.75	0.8414	0.5483
-1.7400000000000002	0.8421	0.5507
-1.7300000000000004	0.8427	0.5536
-1.7200000000000006	0.8435	0.5561
-1.709999999999999	0.8441	0.5581
-1.6999999999999993	0.8448	0.5607
-1.6899999999999995	0.8455	0.5631
-1.6799999999999997	0.8462	0.5655
-1.67	0.8467	0.5682
-1.6600000000000001	0.8475	0.5706
-1.6500000000000004	0.8481	0.5731
-1.6400000000000006	0.8486	0.5757
-1.629999999999999	0.8493	0.5779
-1.6199999999999992	0.85	0.5804
-1.6099999999999994	0.8504	0.5828
-1.5999999999999996	0.8511	0.5851
-1.5899999999999999	0.8517	0.5875
-1.58	0.8524	0.5901
-1.5700000000000003	0.853	0.5922
-1.5600000000000005	0.8537	0.5945
-1.549999999999999	0.8542	0.5967
-1.5399999999999991	0.855	0.5992
-1.5299999999999994	0.8556	0.6014
-1.5199999999999996	0.8562	0.6038
-1.5099999999999998	0.8569	0.6059
-1.5	0.8576	0.6082
-1.4900000000000002	0.8582	0.6104
-1.4800000000000004	0.8588	0.6126
-1.4700000000000006	0.8595	0.6156
-1.459999999999999	0.86	0.6175
-1.4499999999999993	0.8606	0.6193
-1.4399999999999995	0.861	0.6213
-1.4299999999999997	0.8616	0.6235
-1.42	0.8621	0.6257
-1.4100000000000001	0.8626	0.6285
-1.4000000000000004	0.8632	0.6305
-1.3900000000000006	0.8638	0.6328
-1.379999999999999	0.8643	0.635
-1.3699999999999992	0.865	0.6371
-1.3599999999999994	0.8656	0.6391
-1.3499999999999996	0.8661	0.641
-1.3399999999999999	0.8668	0.6429
-1.33	0.8674	0.6448
-1.3200000000000003	0.868	0.6467
-1.3100000000000005	0.8686	0.6487
-1.299999999999999	0.8691	0.6506
-1.2899999999999991	0.8697	0.6529
-1.2799999999999994	0.8703	0.6548
-1.2699999999999996	0.8709	0.6567
-1.2599999999999998	0.8714	0.6584
-1.25	0.8719	0.6602
-1.2400000000000002	0.8723	0.662
-1.2300000000000004	0.8729	0.6639
-1.2200000000000006	0.8733	0.6658
-1.209999999999999	0.874	0.6676
-1.1999999999999993	0.8746	0.6694
-1.1899999999999995	0.875	0.6715
-1.1799999999999997	0.8755	0.6734
-1.17	0.8762	0.6754
-1.1600000000000001	0.8769	0.6773
-1.1500000000000004	0.8774	0.6792
-1.1400000000000006	0.8781	0.681
-1.129999999999999	0.8786	0.6826
-1.1199999999999992	0.8791	0.6841
-1.1099999999999994	0.8797	0.6857
-1.0999999999999996	0.8804	0.6873
-1.0899999999999999	0.8807	0.6893
-1.08	0.8812	0.6909
-1.0700000000000003	0.8817	0.6925
-1.0600000000000005	0.8824	0.6943
-1.049999999999999	0.8829	0.696
-1.0399999999999991	0.8836	0.6975
-1.0299999999999994	0.8843	0.6992
-1.0199999999999996	0.8848	0.7007
-1.0099999999999998	0.8854	0.7026
-1.0	0.8859	0.7041
-0.9900000000000002	0.8864	0.7055
-0.9800000000000004	0.887	0.7071
-0.9700000000000006	0.8875	0.7089
-0.9599999999999991	0.8883	0.7104
-0.9499999999999993	0.8888	0.7122
-0.9399999999999995	0.8896	0.714
-0.9299999999999997	0.8901	0.7155
-0.9199999999999999	0.8906	0.7169
-0.9100000000000001	0.8911	0.7183
-0.9000000000000004	0.8916	0.7199
-0.8900000000000006	0.8922	0.7213
-0.879999999999999	0.8927	0.7229
-0.8699999999999992	0.8933	0.7242
-0.8599999999999994	0.8938	0.7257
-0.8499999999999996	0.8943	0.7272
-0.8399999999999999	0.8948	0.7289
-0.8300000000000001	0.8953	0.7304
-0.8200000000000003	0.8958	0.7317
-0.8100000000000005	0.8962	0.7331
-0.7999999999999989	0.8967	0.7349
-0.7899999999999991	0.8972	0.7362
-0.7799999999999994	0.8976	0.7377
-0.7699999999999996	0.8982	0.7391
-0.7599999999999998	0.8987	0.7404
-0.75	0.8992	0.7417
-0.7400000000000002	0.8996	0.7431
-0.7300000000000004	0.9002	0.7445
-0.7200000000000006	0.9006	0.746
-0.7099999999999991	0.9012	0.7473
-0.6999999999999993	0.9017	0.7487
-0.6899999999999995	0.9022	0.7499
-0.6799999999999997	0.9027	0.7513
-0.6699999999999999	0.9031	0.7525
-0.6600000000000001	0.9035	0.7538
-0.6500000000000004	0.904	0.7554
-0.6400000000000006	0.9045	0.7568
-0.629999999999999	0.9051	0.7581
-0.6199999999999992	0.9055	0.7597
-0.6099999999999994	0.906	0.7612
-0.5999999999999996	0.9065	0.7627
-0.5899999999999999	0.9069	0.7642
-0.5800000000000001	0.9073	0.7655
-0.5700000000000003	0.9078	0.7669
-0.5600000000000005	0.9082	0.7681
-0.5499999999999989	0.9085	0.7692
-0.5399999999999991	0.9092	0.7703
-0.5299999999999994	0.9096	0.7718
-0.5199999999999996	0.9101	0.773
-0.5099999999999998	0.9105	0.7742
-0.5	0.911	0.7752
-0.4900000000000002	0.9113	0.7764
-0.4800000000000004	0.9117	0.7777
-0.47000000000000064	0.9122	0.7789
-0.4599999999999991	0.9126	0.7801
-0.4499999999999993	0.9134	0.7812
-0.4399999999999995	0.9139	0.7826
-0.4299999999999997	0.9143	0.7837
-0.41999999999999993	0.9148	0.7848
-0.41000000000000014	0.9152	0.7859
-0.40000000000000036	0.9157	0.7871
-0.39000000000000057	0.9161	0.7883
-0.379999999999999	0.9165	0.7893
-0.3699999999999992	0.9169	0.7904
-0.35999999999999943	0.9173	0.7916
-0.34999999999999964	0.9177	0.7926
-0.33999999999999986	0.918	0.7936
-0.33000000000000007	0.9183	0.7946
-0.3200000000000003	0.9187	0.7959
-0.3100000000000005	0.9191	0.7969
-0.29999999999999893	0.9195	0.7982
-0.28999999999999915	0.9199	0.7993
-0.27999999999999936	0.9203	0.8003
-0.2699999999999996	0.9206	0.8014
-0.2599999999999998	0.9211	0.8024
-0.25	0.9214	0.8034
-0.2400000000000002	0.9218	0.8046
-0.23000000000000043	0.9222	0.8058
-0.22000000000000064	0.9225	0.8068
-0.20999999999999908	0.9229	0.8079
-0.1999999999999993	0.9233	0.8089
-0.1899999999999995	0.9236	0.81
-0.17999999999999972	0.9239	0.8112
-0.16999999999999993	0.9242	0.8121
-0.16000000000000014	0.9245	0.8132
-0.15000000000000036	0.9248	0.8143
-0.14000000000000057	0.9251	0.8152
-0.129999999999999	0.9256	0.8162
-0.11999999999999922	0.9259	0.8171
-0.10999999999999943	0.9262	0.8181
-0.09999999999999964	0.9266	0.819
-0.08999999999999986	0.9268	0.8201
-0.08000000000000007	0.9271	0.8211
-0.07000000000000028	0.9274	0.8222
-0.0600000000000005	0.9277	0.8231
-0.049999999999998934	0.9281	0.8241
-0.03999999999999915	0.9284	0.8251
-0.02999999999999936	0.9287	0.8261
-0.019999999999999574	0.9291	0.8268
-0.009999999999999787	0.9294	0.8278
0.0	0.9298	0.8286
0.009999999999999787	0.9301	0.8296
0.019999999999999574	0.9304	0.8306
0.030000000000001137	0.9307	0.8315
0.040000000000000924	0.9309	0.8322
0.05000000000000071	0.9313	0.8333
0.0600000000000005	0.9317	0.8341
0.07000000000000028	0.932	0.835
0.08000000000000007	0.9324	0.8359
0.08999999999999986	0.9326	0.8366
0.09999999999999964	0.933	0.8374
0.10999999999999943	0.9333	0.8382
0.120000000000001	0.9336	0.8391
0.13000000000000078	0.9339	0.84
0.14000000000000057	0.9341	0.8408
0.15000000000000036	0.9344	0.8416
0.16000000000000014	0.9348	0.8427
0.16999999999999993	0.9351	0.8436
0.17999999999999972	0.9355	0.8443
0.1899999999999995	0.9358	0.8451
0.20000000000000107	0.9361	0.8461
0.21000000000000085	0.9364	0.8471
0.22000000000000064	0.9368	0.8479
0.23000000000000043	0.9371	0.8488
0.2400000000000002	0.9374	0.8496
0.25	0.9378	0.8503
0.2599999999999998	0.9381	0.8511
0.2699999999999996	0.9385	0.8518
0.28000000000000114	0.9389	0.8526
0.2900000000000009	0.9391	0.8534
0.3000000000000007	0.9394	0.8541
0.3100000000000005	0.9397	0.8549
0.3200000000000003	0.94	0.8558
0.33000000000000007	0.9403	0.8567
0.33999999999999986	0.9406	0.8574
0.34999999999999964	0.9409	0.8581
0.35999999999999943	0.9412	0.8587
0.370000000000001	0.9415	0.8596
0.3800000000000008	0.9418	0.8605
0.39000000000000057	0.9421	0.8613
0.40000000000000036	0.9423	0.8622
0.41000000000000014	0.9426	0.8629
0.41999999999999993	0.9429	0.8637
0.4299999999999997	0.9433	0.8643
0.4399999999999995	0.9434	0.8651
0.45000000000000107	0.9438	0.8659
0.46000000000000085	0.9441	0.8667
0.47000000000000064	0.9443	0.8673
0.4800000000000004	0.9446	0.8681
0.4900000000000002	0.9451	0.8688
0.5	0.9453	0.8696
0.5099999999999998	0.9456	0.8702
0.5199999999999996	0.9459	0.871
0.5300000000000011	0.9461	0.8718
0.5400000000000009	0.9463	0.8725
0.5500000000000007	0.9465	0.8733
0.5600000000000005	0.9467	0.8739
0.5700000000000003	0.947	0.8745
0.5800000000000001	0.9473	0.8755
0.5899999999999999	0.9475	0.8762
0.5999999999999996	0.9477	0.877
0.6099999999999994	0.9479	0.8777
0.620000000000001	0.9482	0.8784
0.6300000000000008	0.9485	0.8791
0.6400000000000006	0.9487	0.8799
0.6500000000000004	0.949	0.8805
0.6600000000000001	0.9492	0.8811
0.6699999999999999	0.9494	0.8819
0.6799999999999997	0.9496	0.8825
0.6899999999999995	0.95	0.8832
0.7000000000000011	0.9502	0.8839
0.7100000000000009	0.9505	0.8846
0.7200000000000006	0.9507	0.8852
0.7300000000000004	0.9509	0.8859
0.7400000000000002	0.9512	0.8865
0.75	0.9514	0.8871
0.7599999999999998	0.9516	0.8877
0.7699999999999996	0.9518	0.8882
0.7800000000000011	0.9519	0.8888
0.7900000000000009	0.9521	0.8894
0.8000000000000007	0.9523	0.89
0.8100000000000005	0.9525	0.8905
0.8200000000000003	0.9528	0.8912
0.8300000000000001	0.953	0.8918
0.8399999999999999	0.9533	0.8924
0.8499999999999996	0.9535	0.8932
0.8599999999999994	0.9537	0.8937
0.870000000000001	0.9539	0.8943
0.8800000000000008	0.9541	0.895
0.8900000000000006	0.9543	0.8958
0.9000000000000004	0.9545	0.8964
0.9100000000000001	0.9549	0.8971
0.9199999999999999	0.9551	0.8977
0.9299999999999997	0.9553	0.8983
0.9399999999999995	0.9555	0.8989
0.9500000000000011	0.9558	0.8995
0.9600000000000009	0.9559	0.9
0.9700000000000006	0.9562	0.9007
0.9800000000000004	0.9563	0.9013
0.9900000000000002	0.9565	0.9018
1.0	0.9567	0.9025
1.0099999999999998	0.9569	0.9031
1.0199999999999996	0.9571	0.9037
1.0300000000000011	0.9573	0.9043
1.040000000000001	0.9575	0.9047
1.0500000000000007	0.9577	0.9053
1.0600000000000005	0.9579	0.9059
1.0700000000000003	0.9581	0.9066
1.08	0.9583	0.9072
1.0899999999999999	0.9585	0.9077
1.0999999999999996	0.9587	0.9083
1.1099999999999994	0.9589	0.9088
1.120000000000001	0.9591	0.9093
1.1300000000000008	0.9593	0.9097
1.1400000000000006	0.9595	0.9105
1.1500000000000004	0.9596	0.9109
1.1600000000000001	0.9599	0.9114
1.17	0.96	0.9118
1.1799999999999997	0.9601	0.9124
1.1899999999999995	0.9603	0.9128
1.200000000000001	0.9605	0.9132
1.2100000000000009	0.9608	0.9136
1.2200000000000006	0.961	0.9142
1.2300000000000004	0.9611	0.9147
1.2400000000000002	0.9614	0.9151
1.25	0.9616	0.9155
1.2599999999999998	0.9618	0.916
1.2699999999999996	0.962	0.9164
1.2800000000000011	0.9621	0.9169
1.290000000000001	0.9623	0.9176
1.3000000000000007	0.9625	0.9181
1.3100000000000005	0.9627	0.9187
1.3200000000000003	0.9628	0.9191
1.33	0.963	0.9195
1.3399999999999999	0.9632	0.92
1.3499999999999996	0.9633	0.9204
1.3599999999999994	0.9635	0.9208
1.370000000000001	0.9636	0.9214
1.3800000000000008	0.9639	0.922
1.3900000000000006	0.9641	0.9224
1.4000000000000004	0.9641	0.9229
1.4100000000000001	0.9643	0.9233
1.42	0.9645	0.9237
1.4299999999999997	0.9646	0.9243
1.4399999999999995	0.9648	0.9247
1.450000000000001	0.9649	0.9252
1.4600000000000009	0.965	0.9258
1.4700000000000006	0.9651	0.9263
1.4800000000000004	0.9652	0.9268
1.4900000000000002	0.9655	0.9274
1.5	0.9655	0.9278
1.5099999999999998	0.9658	0.9283
1.5199999999999996	0.9659	0.9287
1.5300000000000011	0.966	0.9292
1.540000000000001	0.9661	0.9296
1.5500000000000007	0.9663	0.93
1.5600000000000005	0.9664	0.9306
1.5700000000000003	0.9666	0.931
1.58	0.9667	0.9313
1.5899999999999999	0.9668	0.9318
1.5999999999999996	0.9669	0.9324
1.6099999999999994	0.9671	0.9327
1.620000000000001	0.9673	0.9331
1.6300000000000008	0.9675	0.9336
1.6400000000000006	0.9677	0.934
1.6500000000000004	0.9678	0.9343
1.6600000000000001	0.9679	0.9347
1.67	0.9681	0.9352
1.6799999999999997	0.9683	0.9356
1.6899999999999995	0.9684	0.9361
1.700000000000001	0.9685	0.9365
1.7100000000000009	0.9686	0.9369
1.7200000000000006	0.9688	0.9374
1.7300000000000004	0.9688	0.9378
1.7400000000000002	0.9689	0.9381
1.75	0.969	0.9385
1.7599999999999998	0.9691	0.939
1.7699999999999996	0.9693	0.9394
1.7800000000000011	0.9694	0.9398
1.790000000000001	0.9695	0.9401
1.8000000000000007	0.9696	0.9405
1.8100000000000005	0.9697	0.9409
1.8200000000000003	0.9698	0.9412
1.83	0.9698	0.9415
1.8399999999999999	0.9699	0.9419
1.8499999999999996	0.97	0.9422
1.8599999999999994	0.9701	0.9425
1.870000000000001	0.9703	0.9428
1.8800000000000008	0.9705	0.9434
1.8900000000000006	0.9705	0.9438
1.9000000000000004	0.9706	0.9443
1.9100000000000001	0.9707	0.9446
1.92	0.9708	0.945
1.9299999999999997	0.9709	0.9453
1.9399999999999995	0.971	0.9459
1.950000000000001	0.9711	0.9464
1.9600000000000009	0.9712	0.9468
1.9700000000000006	0.9713	0.947
1.9800000000000004	0.9713	0.9473
1.9900000000000002	0.9714	0.9476
2.0	0.9714	0.948
2.01	0.9715	0.9482
2.0199999999999996	0.9716	0.9485
2.030000000000001	0.9717	0.9489
2.040000000000001	0.9718	0.9493
2.0500000000000007	0.9719	0.9497
2.0600000000000005	0.972	0.9499
2.0700000000000003	0.9721	0.9502
2.08	0.9722	0.9506
2.09	0.9723	0.9509
2.0999999999999996	0.9725	0.9513
2.1099999999999994	0.9726	0.9516
2.120000000000001	0.9726	0.952
2.130000000000001	0.9727	0.9525
2.1400000000000006	0.9728	0.9528
2.1500000000000004	0.9729	0.9531
2.16	0.973	0.9533
2.17	0.973	0.9536
2.1799999999999997	0.9732	0.9538
2.1899999999999995	0.9732	0.9541
2.200000000000001	0.9733	0.9544
2.210000000000001	0.9734	0.9548
2.2200000000000006	0.9735	0.9551
2.2300000000000004	0.9736	0.9555
2.24	0.9737	0.9557
2.25	0.9737	0.9561
2.26	0.9738	0.9564
2.2699999999999996	0.974	0.9567
2.280000000000001	0.974	0.9569
2.290000000000001	0.9741	0.9571
2.3000000000000007	0.9742	0.9574
2.3100000000000005	0.9744	0.9577
2.3200000000000003	0.9744	0.958
2.33	0.9745	0.9582
2.34	0.9746	0.9585
2.3499999999999996	0.9747	0.9587
2.3599999999999994	0.9748	0.9591
2.370000000000001	0.9748	0.9594
2.380000000000001	0.975	0.9597
2.3900000000000006	0.9751	0.96
2.4000000000000004	0.9752	0.9603
2.41	0.9753	0.9606
2.42	0.9759	0.961
2.4299999999999997	0.9759	0.9612
2.4399999999999995	0.976	0.9615
2.450000000000001	0.976	0.9618
2.460000000000001	0.9761	0.9621
2.4700000000000006	0.9761	0.9624
2.4800000000000004	0.9762	0.9626
2.49	0.9763	0.9628
2.5	0.9763	0.9632
2.51	0.9763	0.9636
2.5199999999999996	0.9764	0.9639
2.530000000000001	0.9764	0.9641
2.540000000000001	0.9767	0.9645
2.5500000000000007	0.9772	0.9647
2.5600000000000005	0.9782	0.965
2.5700000000000003	0.9784	0.9652
2.58	0.9785	0.9655
2.59	0.9787	0.9658
2.5999999999999996	0.9789	0.9661
2.6099999999999994	0.9789	0.9663
2.620000000000001	0.979	0.9665
2.630000000000001	0.9791	0.9668
2.6400000000000006	0.9792	0.9672
2.6500000000000004	0.9793	0.9675
2.66	0.9793	0.9677
2.67	0.9793	0.968
2.6799999999999997	0.9794	0.9682
2.6899999999999995	0.9796	0.9685
2.700000000000001	0.9797	0.9688
2.710000000000001	0.9835	0.9691
2.7200000000000006	0.9836	0.9694
2.7300000000000004	0.9837	0.9696
2.74	0.9838	0.9697
2.75	0.9841	0.97
2.76	0.9841	0.9702
2.7699999999999996	0.9842	0.9704
2.780000000000001	0.9842	0.9706
2.790000000000001	0.9842	0.9709
2.8000000000000007	0.9843	0.9711
2.8100000000000005	0.9844	0.9713
2.8200000000000003	0.9844	0.9715
2.83	0.9845	0.9718
2.84	0.9847	0.9721
2.8499999999999996	0.9851	0.9723
2.8599999999999994	0.9851	0.9726
2.870000000000001	0.9852	0.9728
2.880000000000001	0.9852	0.973
2.8900000000000006	0.9852	0.9732
2.9000000000000004	0.9859	0.9735
2.91	0.9864	0.9737
2.92	0.9865	0.9738
2.9299999999999997	0.9866	0.974
2.9399999999999995	0.9868	0.9743
2.950000000000001	0.987	0.9745
2.960000000000001	0.987	0.9747
2.9700000000000006	0.9871	0.9749
2.9800000000000004	0.9873	0.9752
2.99	0.9873	0.9753
3.0	0.9875	0.9755
3.01	0.9877	0.9756
3.0199999999999996	0.9878	0.9758
3.030000000000001	0.9878	0.976
3.040000000000001	0.9879	0.9762
3.0500000000000007	0.9879	0.9766
3.0600000000000005	0.988	0.9768
3.0700000000000003	0.9888	0.977
3.08	0.9896	0.9771
3.09	0.9897	0.9773
3.0999999999999996	0.9906	0.9775
3.1099999999999994	0.9907	0.9777
3.120000000000001	0.9916	0.9779
3.130000000000001	0.9916	0.9781
3.1400000000000006	0.9916	0.9783
3.1500000000000004	0.9917	0.9785
3.16	0.9917	0.9786
3.17	0.9917	0.9787
3.1799999999999997	0.9917	0.9789
3.1899999999999995	0.9918	0.9791
3.200000000000001	0.9918	0.9793
3.210000000000001	0.9919	0.9794
3.2200000000000006	0.9921	0.9796
3.2300000000000004	0.9921	0.9798
3.24	0.9923	0.98
3.25	0.9924	0.9801
3.26	0.9931	0.9803
3.2699999999999996	0.9931	0.9805
3.280000000000001	0.9931	0.9806
3.290000000000001	0.9932	0.9808
3.3000000000000007	0.9935	0.9809
3.3100000000000005	0.9935	0.9811
3.3200000000000003	0.9936	0.9812
3.33	0.9937	0.9814
3.34	0.9937	0.9816
3.3499999999999996	0.9938	0.9818
3.3599999999999994	0.9938	0.9819
3.370000000000001	0.9938	0.9821
3.380000000000001	0.9938	0.9822
3.3900000000000006	0.9939	0.9824
3.4000000000000004	0.9939	0.9825
3.41	0.994	0.9826
3.42	0.994	0.9828
3.4299999999999997	0.994	0.9829
3.4399999999999995	0.9941	0.983
3.450000000000001	0.9943	0.983
3.460000000000001	0.9943	0.9832
3.4700000000000006	0.9943	0.9833
3.4800000000000004	0.9943	0.9835
3.49	0.9947	0.9836
3.5	0.9947	0.9839
3.51	0.9949	0.984
3.5199999999999996	0.9949	0.9841
3.530000000000001	0.995	0.9841
3.540000000000001	0.9951	0.9842
3.5500000000000007	0.9952	0.9844
3.5600000000000005	0.9952	0.9845
3.5700000000000003	0.9953	0.9846
3.58	0.9953	0.9847
3.59	0.9954	0.9849
3.5999999999999996	0.9954	0.985
3.610000000000001	0.9954	0.9852
3.620000000000001	0.9954	0.9854
3.630000000000001	0.9957	0.9854
3.6400000000000006	0.9957	0.9855
3.6500000000000004	0.9958	0.9856
3.66	0.9958	0.9858
3.67	0.9959	0.9859
3.6799999999999997	0.9959	0.986
3.6899999999999995	0.9959	0.9862
3.700000000000001	0.996	0.9864
3.710000000000001	0.996	0.9864
3.7200000000000006	0.996	0.9865
3.7300000000000004	0.996	0.9867
3.74	0.9961	0.9868
3.75	0.9961	0.9869
3.76	0.9962	0.987
3.7699999999999996	0.9962	0.9871
3.780000000000001	0.9962	0.9872
3.790000000000001	0.9963	0.9873
3.8000000000000007	0.9963	0.9874
3.8100000000000005	0.9964	0.9875
3.8200000000000003	0.9964	0.9876
3.83	0.9964	0.9877
3.84	0.9964	0.9879
3.8499999999999996	0.9965	0.9881
3.860000000000001	0.9965	0.9881
3.870000000000001	0.9965	0.9882
3.880000000000001	0.9966	0.9883
3.8900000000000006	0.9966	0.9884
3.9000000000000004	0.9966	0.9886
3.91	0.9966	0.9887
3.92	0.9967	0.9888
3.9299999999999997	0.9967	0.9889
3.9399999999999995	0.9968	0.989
3.950000000000001	0.9968	0.9891
3.960000000000001	0.9968	0.9892
3.9700000000000006	0.9968	0.9893
3.9800000000000004	0.9969	0.9894
3.99	0.9969	0.9895
4.0	0.9969	0.9896
4.010000000000002	0.9969	0.9897
4.02	0.997	0.9898
4.030000000000001	0.9971	0.9899
4.039999999999999	0.9971	0.99
4.050000000000001	0.9971	0.9901
4.059999999999999	0.9972	0.9902
4.07	0.9972	0.9903
4.080000000000002	0.9972	0.9904
4.09	0.9973	0.9904
4.100000000000001	0.9973	0.9906
4.109999999999999	0.9974	0.9907
4.120000000000001	0.9974	0.9908
4.129999999999999	0.9974	0.9909
4.140000000000001	0.9975	0.991
4.149999999999999	0.9975	0.9911
4.16	0.9975	0.9912
4.170000000000002	0.9975	0.9913
4.18	0.9975	0.9915
4.190000000000001	0.9976	0.9916
4.199999999999999	0.9976	0.9917
4.210000000000001	0.9976	0.9917
4.219999999999999	0.9976	0.9918
4.23	0.9976	0.9919
4.240000000000002	0.9977	0.9919
4.25	0.9977	0.9921
4.260000000000002	0.9977	0.9921
4.27	0.9977	0.9922
4.280000000000001	0.9978	0.9923
4.289999999999999	0.9978	0.9923
4.300000000000001	0.9978	0.9924
4.309999999999999	0.9978	0.9925
4.32	0.9979	0.9926
4.330000000000002	0.9979	0.9926
4.34	0.998	0.9927
4.350000000000001	0.998	0.9928
4.359999999999999	0.998	0.9928
4.370000000000001	0.998	0.9928
4.379999999999999	0.998	0.9929
4.390000000000001	0.998	0.9929
4.399999999999999	0.998	0.993
4.41	0.9981	0.9931
4.420000000000002	0.9981	0.9932
4.43	0.9981	0.9932
4.440000000000001	0.9981	0.9933
4.449999999999999	0.9981	0.9934
4.460000000000001	0.9981	0.9935
4.469999999999999	0.9981	0.9936
4.48	0.9982	0.9937
4.490000000000002	0.9982	0.9938
4.5	0.9982	0.9939
4.510000000000002	0.9983	0.9939
4.52	0.9983	0.994
4.530000000000001	0.9983	0.9941
4.539999999999999	0.9983	0.9941
4.550000000000001	0.9983	0.9942
4.559999999999999	0.9984	0.9943
4.57	0.9984	0.9943
4.580000000000002	0.9984	0.9944
4.59	0.9984	0.9945
4.600000000000001	0.9984	0.9946
4.609999999999999	0.9985	0.9946
4.620000000000001	0.9985	0.9947
4.629999999999999	0.9985	0.9947
4.640000000000001	0.9985	0.9948
4.649999999999999	0.9985	0.9948
4.66	0.9985	0.9949
4.670000000000002	0.9985	0.9949
4.68	0.9986	0.9949
4.690000000000001	0.9986	0.995
4.699999999999999	0.9986	0.9951
4.710000000000001	0.9986	0.9951
4.719999999999999	0.9986	0.9952
4.73	0.9986	0.9953
4.740000000000002	0.9986	0.9954
4.75	0.9986	0.9954
4.760000000000002	0.9987	0.9954
4.77	0.9987	0.9955
4.780000000000001	0.9987	0.9956
4.789999999999999	0.9987	0.9956
4.800000000000001	0.9987	0.9957
4.809999999999999	0.9987	0.9957
4.82	0.9987	0.9958
4.830000000000002	0.9987	0.9958
4.84	0.9987	0.9959
4.850000000000001	0.9987	0.9959
4.859999999999999	0.9988	0.996
4.870000000000001	0.9988	0.996
4.879999999999999	0.9988	0.9961
4.890000000000001	0.9988	0.9961
4.899999999999999	0.9988	0.9962
4.91	0.9988	0.9962
4.920000000000002	0.9988	0.9963
4.93	0.9988	0.9963
4.940000000000001	0.9988	0.9964
4.949999999999999	0.9989	0.9964
4.960000000000001	0.9989	0.9964
4.969999999999999	0.9989	0.9964
4.98	0.9989	0.9965
4.990000000000002	0.9989	0.9965
5.0	0.9989	0.9966
5.010000000000002	0.9989	0.9966
5.02	0.999	0.9967
5.030000000000001	0.999	0.9967
5.039999999999999	0.999	0.9968
5.050000000000001	0.999	0.9968
5.059999999999999	0.999	0.9968
5.07	0.999	0.9969
5.080000000000002	0.999	0.9969
5.09	0.999	0.997
5.100000000000001	0.9991	0.9971
5.109999999999999	0.9991	0.9971
5.120000000000001	0.9991	0.9972
5.129999999999999	0.9991	0.9972
5.140000000000001	0.9991	0.9972
5.150000000000002	0.9991	0.9972
5.16	0.9991	0.9972
5.170000000000002	0.9991	0.9973
5.18	0.9991	0.9973
5.190000000000001	0.9992	0.9973
5.199999999999999	0.9992	0.9973
5.210000000000001	0.9992	0.9974
5.219999999999999	0.9992	0.9974
5.23	0.9992	0.9975
5.240000000000002	0.9992	0.9975
5.25	0.9992	0.9975
5.260000000000002	0.9992	0.9976
5.27	0.9992	0.9976
5.280000000000001	0.9993	0.9976
5.289999999999999	0.9993	0.9977
5.300000000000001	0.9993	0.9977
5.309999999999999	0.9993	0.9977
5.32	0.9993	0.9978
5.330000000000002	0.9993	0.9978
5.34	0.9993	0.9978
5.350000000000001	0.9994	0.9979
5.359999999999999	0.9994	0.9979
5.370000000000001	0.9994	0.9979
5.379999999999999	0.9994	0.9979
5.390000000000001	0.9994	0.998
5.400000000000002	0.9994	0.998
5.41	0.9994	0.998
5.420000000000002	0.9994	0.9981
5.43	0.9994	0.9981
5.440000000000001	0.9995	0.9981
5.449999999999999	0.9995	0.9981
5.460000000000001	0.9995	0.9982
5.469999999999999	0.9995	0.9982
5.48	0.9995	0.9982
5.490000000000002	0.9995	0.9983
5.5	0.9995	0.9983
5.510000000000002	0.9995	0.9983
5.52	0.9995	0.9984
5.530000000000001	0.9995	0.9984
5.539999999999999	0.9995	0.9985
5.550000000000001	0.9995	0.9985
5.559999999999999	0.9995	0.9985
5.57	0.9995	0.9985
5.580000000000002	0.9995	0.9985
5.59	0.9995	0.9986
5.600000000000001	0.9995	0.9986
5.609999999999999	0.9995	0.9986
5.620000000000001	0.9995	0.9986
5.629999999999999	0.9995	0.9986
5.640000000000001	0.9996	0.9987
5.650000000000002	0.9996	0.9987
5.66	0.9996	0.9987
5.670000000000002	0.9996	0.9988
5.68	0.9996	0.9988
5.690000000000001	0.9996	0.9988
5.699999999999999	0.9996	0.9989
5.710000000000001	0.9996	0.9989
5.719999999999999	0.9996	0.9989
5.73	0.9996	0.9989
5.740000000000002	0.9996	0.9989
5.75	0.9996	0.9989
5.760000000000002	0.9996	0.999
5.77	0.9996	0.999
5.780000000000001	0.9996	0.999
5.789999999999999	0.9996	0.9991
5.800000000000001	0.9997	0.9991
5.809999999999999	0.9997	0.9991
5.82	0.9997	0.9991
5.830000000000002	0.9997	0.9991
5.84	0.9997	0.9991
5.850000000000001	0.9997	0.9991
5.859999999999999	0.9997	0.9992
5.870000000000001	0.9997	0.9992
5.879999999999999	0.9997	0.9992
5.890000000000001	0.9997	0.9992
5.900000000000002	0.9997	0.9992
5.91	0.9997	0.9992
5.920000000000002	0.9997	0.9992
5.93	0.9997	0.9992
5.940000000000001	0.9997	0.9992
5.949999999999999	0.9997	0.9992
5.960000000000001	0.9997	0.9992
5.969999999999999	0.9998	0.9992
5.98	0.9998	0.9992
5.990000000000002	0.9998	0.9992
6.0	0.9998	0.9992
6.010000000000002	0.9998	0.9992
6.02	0.9998	0.9992
6.030000000000001	0.9998	0.9993
6.039999999999999	0.9998	0.9993
6.050000000000001	0.9998	0.9993
6.059999999999999	0.9998	0.9993
6.07	0.9998	0.9993
6.080000000000002	0.9998	0.9993
6.09	0.9998	0.9993
6.100000000000001	0.9998	0.9994
6.109999999999999	0.9998	0.9994
6.120000000000001	0.9998	0.9994
6.129999999999999	0.9998	0.9994
6.140000000000001	0.9998	0.9994
6.150000000000002	0.9999	0.9994
6.16	0.9999	0.9995
6.170000000000002	0.9999	0.9995
6.18	0.9999	0.9995
6.190000000000001	0.9999	0.9995
6.199999999999999	0.9999	0.9995
6.210000000000001	0.9999	0.9995
6.219999999999999	0.9999	0.9995
6.23	0.9999	0.9995
6.240000000000002	0.9999	0.9996
6.25	0.9999	0.9996
6.260000000000002	0.9999	0.9996
6.27	0.9999	0.9996
6.280000000000001	0.9999	0.9996
6.289999999999999	0.9999	0.9996
6.300000000000001	0.9999	0.9996
6.309999999999999	0.9999	0.9996
6.32	0.9999	0.9996
6.330000000000002	0.9999	0.9996
6.34	0.9999	0.9996
6.350000000000001	0.9999	0.9997
6.359999999999999	0.9999	0.9997
6.370000000000001	0.9999	0.9997
6.379999999999999	0.9999	0.9997
6.390000000000001	0.9999	0.9997
6.400000000000002	0.9999	0.9997
6.41	0.9999	0.9997
6.420000000000002	0.9999	0.9997
6.43	0.9999	0.9998
6.440000000000001	0.9999	0.9998
6.449999999999999	0.9999	0.9998
6.460000000000001	0.9999	0.9998
6.469999999999999	0.9999	0.9998
6.48	0.9999	0.9998
6.490000000000002	0.9999	0.9998
6.5	0.9999	0.9998
6.510000000000002	0.9999	0.9998
6.52	0.9999	0.9998
6.530000000000001	0.9999	0.9998
6.539999999999999	0.9999	0.9998
6.550000000000001	0.9999	0.9998
6.559999999999999	0.9999	0.9998
6.57	0.9999	0.9998
6.580000000000002	0.9999	0.9998
6.59	0.9999	0.9998
6.600000000000001	0.9999	0.9998
6.609999999999999	0.9999	0.9998
6.620000000000001	0.9999	0.9998
6.629999999999999	0.9999	0.9998
6.640000000000001	0.9999	0.9998
6.650000000000002	0.9999	0.9998
6.66	0.9999	0.9999
6.670000000000002	0.9999	0.9999
6.68	0.9999	0.9999
6.690000000000001	0.9999	0.9999
6.699999999999999	0.9999	0.9999
6.710000000000001	0.9999	0.9999
6.719999999999999	0.9999	0.9999
6.73	0.9999	0.9999
6.740000000000002	0.9999	0.9999
6.75	0.9999	0.9999
6.760000000000002	0.9999	0.9999
6.77	0.9999	0.9999
6.780000000000001	0.9999	0.9999
6.789999999999999	0.9999	0.9999
6.800000000000001	0.9999	0.9999
6.809999999999999	0.9999	0.9999
6.82	0.9999	0.9999
6.830000000000002	0.9999	0.9999
6.84	0.9999	0.9999
6.850000000000001	0.9999	0.9999
6.859999999999999	0.9999	0.9999
6.870000000000001	0.9999	0.9999
6.879999999999999	0.9999	0.9999
6.890000000000001	0.9999	1.0
6.900000000000002	0.9999	1.0
6.91	0.9999	1.0
6.920000000000002	0.9999	1.0
6.93	0.9999	1.0
6.940000000000001	0.9999	1.0
6.949999999999999	1.0	1.0
6.960000000000001	1.0	1.0
6.969999999999999	1.0	1.0
6.98	1.0	1.0
6.990000000000002	1.0	1.0
7.0	1.0	1.0
7.010000000000002	1.0	1.0
7.02	1.0	1.0
7.030000000000001	1.0	1.0
7.039999999999999	1.0	1.0
7.050000000000001	1.0	1.0
7.059999999999999	1.0	1.0
7.07	1.0	1.0
7.080000000000002	1.0	1.0
7.09	1.0	1.0
7.100000000000001	1.0	1.0
7.109999999999999	1.0	1.0
7.120000000000001	1.0	1.0
7.129999999999999	1.0	1.0
7.140000000000001	1.0	1.0
7.150000000000002	1.0	1.0
7.16	1.0	1.0
7.170000000000002	1.0	1.0
7.18	1.0	1.0
7.190000000000001	1.0	1.0
7.199999999999999	1.0	1.0
7.210000000000001	1.0	1.0
7.219999999999999	1.0	1.0
7.23	1.0	1.0
7.240000000000002	1.0	1.0
7.25	1.0	1.0
7.260000000000002	1.0	1.0
7.27	1.0	1.0
7.280000000000001	1.0	1.0
7.289999999999999	1.0	1.0
7.300000000000001	1.0	1.0
7.309999999999999	1.0	1.0
7.32	1.0	1.0
7.330000000000002	1.0	1.0
7.34	1.0	1.0
7.350000000000001	1.0	1.0
7.359999999999999	1.0	1.0
7.370000000000001	1.0	1.0
7.379999999999999	1.0	1.0
7.390000000000001	1.0	1.0
7.400000000000002	1.0	1.0
7.41	1.0	1.0
7.420000000000002	1.0	1.0
7.43	1.0	1.0
7.440000000000001	1.0	1.0
7.449999999999999	1.0	1.0
7.460000000000001	1.0	1.0
7.469999999999999	1.0	1.0
7.48	1.0	1.0
7.490000000000002	1.0	1.0
7.5	1.0	1.0
7.510000000000002	1.0	1.0
7.52	1.0	1.0
7.530000000000001	1.0	1.0
7.539999999999999	1.0	1.0
7.550000000000001	1.0	1.0
7.559999999999999	1.0	1.0
7.57	1.0	1.0
7.580000000000002	1.0	1.0
7.59	1.0	1.0
7.600000000000001	1.0	1.0
7.609999999999999	1.0	1.0
7.620000000000001	1.0	1.0
7.629999999999999	1.0	1.0
7.640000000000001	1.0	1.0
7.650000000000002	1.0	1.0
7.66	1.0	1.0
7.670000000000002	1.0	1.0
7.68	1.0	1.0
7.690000000000001	1.0	1.0
7.699999999999999	1.0	1.0
7.710000000000001	1.0	1.0
7.719999999999999	1.0	1.0
7.73	1.0	1.0
7.740000000000002	1.0	1.0
7.75	1.0	1.0
7.760000000000002	1.0	1.0
7.77	1.0	1.0
7.780000000000001	1.0	1.0
7.789999999999999	1.0	1.0
7.800000000000001	1.0	1.0
7.809999999999999	1.0	1.0
7.82	1.0	1.0
7.830000000000002	1.0	1.0
7.84	1.0	1.0
7.850000000000001	1.0	1.0
7.859999999999999	1.0	1.0
7.870000000000001	1.0	1.0
7.879999999999999	1.0	1.0
7.890000000000001	1.0	1.0
7.900000000000002	1.0	1.0
7.91	1.0	1.0
7.920000000000002	1.0	1.0
7.93	1.0	1.0
7.940000000000001	1.0	1.0
7.949999999999999	1.0	1.0
7.960000000000001	1.0	1.0
7.969999999999999	1.0	1.0
7.98	1.0	1.0
7.990000000000002	1.0	1.0
8.0	1.0	1.0
8.010000000000002	1.0	1.0
8.02	1.0	1.0
8.030000000000001	1.0	1.0
8.04	1.0	1.0
8.05	1.0	1.0
8.059999999999999	1.0	1.0
8.07	1.0	1.0
8.080000000000002	1.0	1.0
8.09	1.0	1.0
8.100000000000001	1.0	1.0
8.11	1.0	1.0
8.120000000000001	1.0	1.0
8.129999999999999	1.0	1.0
8.14	1.0	1.0
8.150000000000002	1.0	1.0
8.16	1.0	1.0
8.170000000000002	1.0	1.0
8.18	1.0	1.0
8.190000000000001	1.0	1.0
8.2	1.0	1.0
8.21	1.0	1.0
8.219999999999999	1.0	1.0
8.23	1.0	1.0
8.240000000000002	1.0	1.0
8.25	1.0	1.0
8.260000000000002	1.0	1.0
8.27	1.0	1.0
8.280000000000001	1.0	1.0
8.29	1.0	1.0
8.3	1.0	1.0
8.309999999999999	1.0	1.0
8.32	1.0	1.0
8.330000000000002	1.0	1.0
8.34	1.0	1.0
8.350000000000001	1.0	1.0
8.36	1.0	1.0
8.370000000000001	1.0	1.0
8.379999999999999	1.0	1.0
8.39	1.0	1.0
8.400000000000002	1.0	1.0
8.41	1.0	1.0
8.420000000000002	1.0	1.0
8.43	1.0	1.0
8.440000000000001	1.0	1.0
8.45	1.0	1.0
8.46	1.0	1.0
8.469999999999999	1.0	1.0
8.48	1.0	1.0
8.490000000000002	1.0	1.0
8.5	1.0	1.0
8.510000000000002	1.0	1.0
8.52	1.0	1.0
8.530000000000001	1.0	1.0
8.54	1.0	1.0
8.55	1.0	1.0
8.559999999999999	1.0	1.0
8.57	1.0	1.0
8.580000000000002	1.0	1.0
8.59	1.0	1.0
8.600000000000001	1.0	1.0
8.61	1.0	1.0
8.620000000000001	1.0	1.0
8.629999999999999	1.0	1.0
8.64	1.0	1.0
8.650000000000002	1.0	1.0
8.66	1.0	1.0
8.670000000000002	1.0	1.0
8.68	1.0	1.0
8.690000000000001	1.0	1.0
8.7	1.0	1.0
8.71	1.0	1.0
8.719999999999999	1.0	1.0
8.73	1.0	1.0
8.740000000000002	1.0	1.0
8.75	1.0	1.0
8.760000000000002	1.0	1.0
8.77	1.0	1.0
8.780000000000001	1.0	1.0
8.79	1.0	1.0
8.8	1.0	1.0
8.809999999999999	1.0	1.0
8.82	1.0	1.0
8.830000000000002	1.0	1.0
8.84	1.0	1.0
8.850000000000001	1.0	1.0
8.86	1.0	1.0
8.870000000000001	1.0	1.0
8.879999999999999	1.0	1.0
8.89	1.0	1.0
8.900000000000002	1.0	1.0
8.91	1.0	1.0
8.920000000000002	1.0	1.0
8.93	1.0	1.0
8.940000000000001	1.0	1.0
8.95	1.0	1.0
8.96	1.0	1.0
8.969999999999999	1.0	1.0
8.98	1.0	1.0
8.990000000000002	1.0	1.0
9.0	1.0	1.0
9.010000000000002	1.0	1.0
9.02	1.0	1.0
9.030000000000001	1.0	1.0
9.04	1.0	1.0
9.05	1.0	1.0
9.059999999999999	1.0	1.0
9.07	1.0	1.0
9.080000000000002	1.0	1.0
9.09	1.0	1.0
9.100000000000001	1.0	1.0
9.11	1.0	1.0
9.120000000000001	1.0	1.0
9.129999999999999	1.0	1.0
9.14	1.0	1.0
9.150000000000002	1.0	1.0
9.16	1.0	1.0
9.170000000000002	1.0	1.0
9.18	1.0	1.0
9.190000000000001	1.0	1.0
9.2	1.0	1.0
9.21	1.0	1.0
9.219999999999999	1.0	1.0
9.23	1.0	1.0
9.240000000000002	1.0	1.0
9.25	1.0	1.0
9.260000000000002	1.0	1.0
9.27	1.0	1.0
9.280000000000001	1.0	1.0
9.29	1.0	1.0
9.3	1.0	1.0
9.309999999999999	1.0	1.0
9.32	1.0	1.0
9.330000000000002	1.0	1.0
9.34	1.0	1.0
9.350000000000001	1.0	1.0
9.36	1.0	1.0
9.370000000000001	1.0	1.0
9.379999999999999	1.0	1.0
9.39	1.0	1.0
9.400000000000002	1.0	1.0
9.41	1.0	1.0
9.420000000000002	1.0	1.0
9.43	1.0	1.0
9.440000000000001	1.0	1.0
9.45	1.0	1.0
9.46	1.0	1.0
9.469999999999999	1.0	1.0
9.48	1.0	1.0
9.490000000000002	1.0	1.0
9.5	1.0	1.0
9.510000000000002	1.0	1.0
9.52	1.0	1.0
9.530000000000001	1.0	1.0
9.54	1.0	1.0
9.55	1.0	1.0
9.559999999999999	1.0	1.0
9.57	1.0	1.0
9.580000000000002	1.0	1.0
9.59	1.0	1.0
9.600000000000001	1.0	1.0
9.61	1.0	1.0
9.620000000000001	1.0	1.0
9.629999999999999	1.0	1.0
9.64	1.0	1.0
9.650000000000002	1.0	1.0
9.66	1.0	1.0
9.670000000000002	1.0	1.0
9.68	1.0	1.0
9.690000000000001	1.0	1.0
9.7	1.0	1.0
9.71	1.0	1.0
9.719999999999999	1.0	1.0
9.73	1.0	1.0
9.740000000000002	1.0	1.0
9.75	1.0	1.0
9.760000000000002	1.0	1.0
9.77	1.0	1.0
9.780000000000001	1.0	1.0
9.79	1.0	1.0
9.8	1.0	1.0
9.809999999999999	1.0	1.0
9.82	1.0	1.0
9.830000000000002	1.0	1.0
9.84	1.0	1.0
9.850000000000001	1.0	1.0
9.86	1.0	1.0
9.870000000000001	1.0	1.0
9.879999999999999	1.0	1.0
9.89	1.0	1.0
9.900000000000002	1.0	1.0
9.91	1.0	1.0
9.920000000000002	1.0	1.0
9.93	1.0	1.0
9.940000000000001	1.0	1.0
9.95	1.0	1.0
9.96	1.0	1.0
9.969999999999999	1.0	1.0
9.98	1.0	1.0
9.990000000000002	1.0	1.0
10.0	1.0	1.0
10.010000000000002	1.0	1.0
10.02	1.0	1.0
10.030000000000001	1.0	1.0
10.04	1.0	1.0
10.05	1.0	1.0
10.059999999999999	1.0	1.0
10.07	1.0	1.0
10.080000000000002	1.0	1.0
10.09	1.0	1.0
10.100000000000001	1.0	1.0
10.11	1.0	1.0
10.120000000000001	1.0	1.0
10.129999999999999	1.0	1.0
10.14	1.0	1.0
10.150000000000002	1.0	1.0
10.16	1.0	1.0
10.170000000000002	1.0	1.0
10.18	1.0	1.0
10.190000000000001	1.0	1.0
10.2	1.0	1.0
10.21	1.0	1.0
10.219999999999999	1.0	1.0
10.23	1.0	1.0
10.240000000000002	1.0	1.0
10.25	1.0	1.0
10.260000000000002	1.0	1.0
10.27	1.0	1.0
10.280000000000001	1.0	1.0
10.29	1.0	1.0
10.3	1.0	1.0
10.309999999999999	1.0	1.0
10.32	1.0	1.0
10.330000000000002	1.0	1.0
10.34	1.0	1.0
10.350000000000001	1.0	1.0
10.36	1.0	1.0
10.370000000000001	1.0	1.0
10.379999999999999	1.0	1.0
10.39	1.0	1.0
10.400000000000002	1.0	1.0
10.41	1.0	1.0
10.420000000000002	1.0	1.0
10.43	1.0	1.0
10.440000000000001	1.0	1.0
10.45	1.0	1.0
10.46	1.0	1.0
10.469999999999999	1.0	1.0
10.48	1.0	1.0
10.490000000000002	1.0	1.0
10.5	1.0	1.0
10.510000000000002	1.0	1.0
10.52	1.0	1.0
10.530000000000001	1.0	1.0
10.54	1.0	1.0
10.55	1.0	1.0
10.559999999999999	1.0	1.0
10.57	1.0	1.0
10.580000000000002	1.0	1.0
10.59	1.0	1.0
10.600000000000001	1.0	1.0
10.61	1.0	1.0
10.620000000000001	1.0	1.0
10.629999999999999	1.0	1.0
10.64	1.0	1.0
10.650000000000002	1.0	1.0
10.66	1.0	1.0
10.670000000000002	1.0	1.0
10.68	1.0	1.0
10.690000000000001	1.0	1.0
10.7	1.0	1.0
10.71	1.0	1.0
10.719999999999999	1.0	1.0
10.73	1.0	1.0
10.740000000000002	1.0	1.0
10.75	1.0	1.0
10.760000000000002	1.0	1.0
10.77	1.0	1.0
10.780000000000001	1.0	1.0
10.79	1.0	1.0
10.8	1.0	1.0
10.809999999999999	1.0	1.0
10.82	1.0	1.0
10.830000000000002	1.0	1.0
10.84	1.0	1.0
10.850000000000001	1.0	1.0
10.86	1.0	1.0
10.870000000000001	1.0	1.0
10.879999999999999	1.0	1.0
10.89	1.0	1.0
10.900000000000002	1.0	1.0
10.91	1.0	1.0
10.920000000000002	1.0	1.0
10.93	1.0	1.0
10.940000000000001	1.0	1.0
10.95	1.0	1.0
10.96	1.0	1.0
10.969999999999999	1.0	1.0
10.98	1.0	1.0
10.990000000000002	1.0	1.0
11.0	1.0	1.0
inf	1.0	1.0"""
