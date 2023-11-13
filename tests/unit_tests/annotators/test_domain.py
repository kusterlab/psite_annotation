import io

import pandas as pd
import pytest
from psite_annotation.annotators.annotator_base import MissingColumnsError

from psite_annotation.annotators.domain import DomainAnnotator

# Define the mock input file as a string
mock_input_file = """Uniprot_ACC,Domain,START_POSITION,END_POSITION
Q86VB7-4,internal_repeat_2,476,705
Q8WZ42-10,low_complexity_region,29208,29219
Q8WZ42-10,Pfam:PPAK,10295,10321
O60566-3,low_complexity_region,455,471
Q6ZU52-2,low_complexity_region,542,561
"""  # noqa: E501,B950


@pytest.fixture
def annotator() -> DomainAnnotator:
    """Fixture to create a DomainAnnotator object using a mock input file.

    Returns:
        DomainAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return DomainAnnotator(file)


@pytest.fixture
def input_df() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Proteins": [
                "Q86VB7-4",
                "Q8WZ42-10;O60566-3",
                "Q6ZU52-2",
                "P05067",
            ],
            "Matched proteins": [
                "Q86VB7-4",
                "Q8WZ42-10;O60566-3",
                "Q6ZU52-2",
                "P05067",
            ],
            "Start positions": [
                "480",
                "29218;460",
                "400",
                "100",
            ],
            "End positions": [
                "500",
                "29225;480",
                "420",
                "120",
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
            "Proteins": [
                "Q86VB7-4",
                "Q8WZ42-10;O60566-3",
                "Q6ZU52-2",
                "P05067",
            ],
            "Matched proteins": [
                "Q86VB7-4",
                "Q8WZ42-10;O60566-3",
                "Q6ZU52-2",
                "P05067",
            ],
            "Start positions": [
                "480",
                "29218;460",
                "400",
                "100",
            ],
            "End positions": [
                "500",
                "29225;480",
                "420",
                "120",
            ],
            "Domains": [
                "internal_repeat_2",
                "low_complexity_region",
                "",
                "",
            ],
        }
    )


class TestDomainAnnotator:
    """Test the DomainAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file into `psp_dict`.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        expected_dict = {
            "Q86VB7-4": [(476, 705, "internal_repeat_2")],
            "Q8WZ42-10": [
                (29208, 29219, "low_complexity_region"),
                (10295, 10321, "Pfam:PPAK"),
            ],
            "O60566-3": [(455, 471, "low_complexity_region")],
            "Q6ZU52-2": [(542, 561, "low_complexity_region")],
        }
        assert annotator.domain_dict == expected_dict

    def test_annotate(self, annotator, input_df, expected_output_df):
        """Test that the annotate method correctly adds the domain annotations to the input dataframe as new columns.

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
            "Proteins",
            "Matched proteins",
            "Start positions",
            "End positions",
            "Domains",
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
