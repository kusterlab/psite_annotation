import io
from unittest.mock import patch

import pandas as pd
import pytest

from psite_annotation.annotators.annotator_base import MissingColumnsError
from psite_annotation.annotators.clinical_basket import ClinicalBasketAnnotator

# Define the mock input file as a string
mock_input_file = """BASKET	GENE NAME
PI3K-AKT-mTOR	PIK3CA
PI3K-AKT-mTOR	PIK3R1
EGFR	GRB2
VEGFR	SOS1
VEGFR	SOS1
FGFR	SOS1
"""


@pytest.fixture
def mock_df() -> pd.DataFrame:
    """Fixture to mock the return value of pandas.read_excel.

    Returns:
        ClinicalBasketAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return pd.read_csv(file, sep="\t")


@pytest.fixture
def input_df() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Gene names": [
                "PIK3CA;PIK3R1",
                "PIK3R1",
                "GRB2",
                "SOS1",
                "SOS2",
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
            "Gene names": [
                "PIK3CA;PIK3R1",
                "PIK3R1",
                "GRB2",
                "SOS1",
                "SOS2",
            ],
            "Clinical baskets": [
                "PI3K-AKT-mTOR",
                "PI3K-AKT-mTOR",
                "EGFR",
                "FGFR;VEGFR",
                "",
            ],
        }
    )


class TestClinicalBasketAnnotator:
    """Test the ClinicalBasketAnnotator class."""

    @patch("pandas.read_excel")
    def test_load_annotations(self, read_excel_method, mock_df: pd.DataFrame):
        """Test that the load_annotations method correctly loads the annotations from the input file.

        Args:
            read_excel_method: patched pandas.read_excel method
            mock_df: mocked annotation dataframe
        """
        read_excel_method.return_value = mock_df

        annotator = ClinicalBasketAnnotator("")
        annotator.load_annotations()

        expected_output_df = pd.DataFrame(
            {
                "Gene names": [
                    "PIK3CA",
                    "PIK3R1",
                    "GRB2",
                    "SOS1",
                ],
                "Clinical baskets": [
                    "PI3K-AKT-mTOR",
                    "PI3K-AKT-mTOR",
                    "EGFR",
                    "FGFR;VEGFR",
                ],
            }
        )
        pd.testing.assert_frame_equal(
            annotator.basket_df, expected_output_df, check_like=True
        )

    @patch("pandas.read_excel")
    def test_annotate(
        self,
        read_excel_method,
        mock_df: pd.DataFrame,
        input_df: pd.DataFrame,
        expected_output_df: pd.DataFrame,
    ):
        """Test that the annotate method correctly adds the basket annotations to the input dataframe as new column.

        Args:
            read_excel_method: patched pandas.read_excel method
            mock_df: mocked annotation dataframe
            input_df: Example input dataframe
            expected_output_df: Expected output dataframe

        """
        read_excel_method.return_value = mock_df

        annotator = ClinicalBasketAnnotator("")
        annotator.load_annotations()

        # Annotate the input dataframe
        output_df = annotator.annotate(input_df)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Gene names",
            "Clinical baskets",
        }

        # Assert that the output dataframe has the expected values

        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)

    @patch("pandas.read_excel")
    def test_raise_error_on_missing_input_columns(
        self, read_excel_method, mock_df: pd.DataFrame
    ):
        """Test that a value error is raised if the "Gene names" column is missing.

        Args:
            annotator: Annotator object with mock file loaded
        """
        read_excel_method.return_value = mock_df

        annotator = ClinicalBasketAnnotator("")
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

    @patch("pandas.read_excel")
    def test_annotate_input_df_unchanged(
        self, read_excel_method, mock_df: pd.DataFrame, input_df: pd.DataFrame
    ):
        """Test that the annotate method does not alter the input dataframe.

        Args:
            read_excel_method: patched pandas.read_excel method
            mock_df: mocked annotation dataframe
            input_df: Example input dataframe

        """
        input_df_copy = input_df.copy()

        read_excel_method.return_value = mock_df

        annotator = ClinicalBasketAnnotator("")
        annotator.load_annotations()

        # Annotate the input dataframe
        _ = annotator.annotate(input_df)

        pd.testing.assert_frame_equal(input_df, input_df_copy, check_like=True)

    @patch("pandas.read_excel")
    def test_annotate_inplace(
        self,
        read_excel_method,
        mock_df: pd.DataFrame,
        input_df: pd.DataFrame,
        expected_output_df: pd.DataFrame,
    ):
        """Test that the annotate method can annotate the input dataframe inplace.

        Args:
            read_excel_method: patched pandas.read_excel method
            mock_df: mocked annotation dataframe
            input_df: Example input dataframe
            expected_output_df: Expected output dataframe

        """
        read_excel_method.return_value = mock_df

        annotator = ClinicalBasketAnnotator("")
        annotator.load_annotations()

        # Annotate the input dataframe
        annotator.annotate(input_df, inplace=True)

        pd.testing.assert_frame_equal(input_df, expected_output_df, check_like=True)
