import io

import pandas as pd
import pytest

from psite_annotation.annotators.in_vitro_kinases import InVitroKinasesAnnotator

# Define the mock input file as a string
mock_input_file = """Kinase	UniprotId_Site
HRI	P05198_S52
PKCD	Q9UQL6_S259
PKCD	P18433-2_S204
PKCD	P10415_S70
"""  # noqa: E501,B950


@pytest.fixture
def annotator() -> InVitroKinasesAnnotator:
    """Fixture to create a InVitroKinasesAnnotator object using a mock input file.

    Returns:
        InVitroKinasesAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return InVitroKinasesAnnotator(file)


@pytest.fixture
def annotator_output_gene_names() -> InVitroKinasesAnnotator:
    """Fixture to create a InVitroKinasesAnnotator object using a mock input file.

    Returns:
        InVitroKinasesAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return InVitroKinasesAnnotator(file, output_gene_names=True)


class TestInVitroKinasesAnnotator:
    """Test the InVitroKinasesAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "P05198_S52",
                    "Q9UQL6_S259",
                    "P18433-2_S204",
                    "P10415_S70",
                ],
                "In Vitro Kinases": [
                    "HRI",
                    "PKCD",
                    "PKCD",
                    "PKCD",
                ],
            }
        )
        pd.testing.assert_frame_equal(
            annotator.kinase_df, expected_output_df, check_like=True
        )

    def test_annotate(self, annotator):
        """Test that the annotate method correctly adds the InVitro annotations to the input dataframe as new columns.

        Args:
            annotator: Annotator object with mock file loaded
        """
        # Create a sample input dataframe
        input_df = pd.DataFrame(
            {
                "Site positions": [
                    "P05198_S52;Q9UQL6_S259",
                    "P18433-2_S204",
                    "P10415_S70",
                    "P05067_S100",
                ]
            }
        )

        annotator.load_annotations()

        # Annotate the input dataframe
        output_df = annotator.annotate(input_df)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Site positions",
            "In Vitro Kinases",
        }

        # Assert that the output dataframe has the expected values
        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "P05198_S52;Q9UQL6_S259",
                    "P18433-2_S204",
                    "P10415_S70",
                    "P05067_S100",
                ],
                "In Vitro Kinases": [
                    "HRI;PKCD",
                    "PKCD",
                    "PKCD",
                    "",
                ],
            }
        )
        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)
