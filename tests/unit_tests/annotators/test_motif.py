import io
import re

import pandas as pd
import pytest

from psite_annotation.annotators.annotator_base import MissingColumnsError
from psite_annotation.annotators.motif import MotifAnnotator

# Define the mock input file as a string
mock_input_file = """Identifier	Regex
MOD_CDK_SPK_2	...([st])P[RK]
MOD_CDK_SPxK_1	...([st])P.[KR]
MOD_CDK_SPxxK_3	...([st])P..[RK]
MOD_CK1_1	S..([st])...
MOD_CK2_1	...([st])..E
"""


@pytest.fixture
def annotator() -> MotifAnnotator:
    """Fixture to create a MotifAnnotator object using a mock input file.

    Returns:
        MotifAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return MotifAnnotator(file)


@pytest.fixture
def input_df() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Site sequence context": [
                "AAAsPRAAA",
                "AAAsPRKAA",
                "AAAsPKAAA",
                "AAAsPAAAA",
                "AAASPRAAA",
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
                "AAAsPRAAA",
                "AAAsPRKAA",
                "AAAsPKAAA",
                "AAAsPAAAA",
                "AAASPRAAA",
            ],
            "Motifs": [
                "MOD_CDK_SPK_2",
                "MOD_CDK_SPK_2;MOD_CDK_SPxK_1",
                "MOD_CDK_SPK_2",
                "",
                "",
            ],
        }
    )


class TestMotifAnnotator:
    """Test the MotifAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file into `motif_dict`.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        assert set(annotator.motif_dict.keys()) == {
            "MOD_CDK_SPK_2",
            "MOD_CDK_SPxK_1",
            "MOD_CDK_SPxxK_3",
            "MOD_CK1_1",
            "MOD_CK2_1",
        }

        assert re.match(annotator.motif_dict["MOD_CDK_SPK_2"], "AAAsPRAAA")
        assert re.match(annotator.motif_dict["MOD_CDK_SPK_2"], "AAAsPKAAA")
        assert not re.match(annotator.motif_dict["MOD_CDK_SPK_2"], "AAAsPAAAA")

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
            "Motifs",
        }

        # Assert that the output dataframe has the expected values

        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)

    def test_raise_error_on_missing_input_columns(self, annotator):
        """Test that a value error is raised if the "Matching proteins", "Start positions" or "End positions" column is missing.

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
