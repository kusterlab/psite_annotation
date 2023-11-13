import io

import pandas as pd
import pytest
from psite_annotation.annotators.annotator_base import MissingColumnsError

from psite_annotation.annotators.psp_kinases import PSPKinasesAnnotator

# Define the mock input file as a string
mock_input_file = """121420
Data extracted from PhosphoSitePlus(R)

GENE	KINASE	KIN_ORGANISM	SUB_ACC_ID	SUB_ORGANISM	SUB_MOD_RSD
Pak2	PAK2	rat	P01236	human	S207
EIF2AK1	HRI	human	P05198	human	S52
PRKCD	PKCD	human	P34901	rat	S183
PRKCD	PKCD	human	Q9UQL6	human	S259
PRKCD	PKCD	human	P18433-2	human	S204
PRKCD	PKCD	human	P10415	human	S70
"""  # noqa: E501,B950


@pytest.fixture
def annotator() -> PSPKinasesAnnotator:
    """Fixture to create a PSPKinasesAnnotator object using a mock input file.

    Returns:
        PSPKinasesAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return PSPKinasesAnnotator(file)


@pytest.fixture
def annotator_output_gene_names() -> PSPKinasesAnnotator:
    """Fixture to create a PSPKinasesAnnotator object using a mock input file.

    Returns:
        PSPKinasesAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return PSPKinasesAnnotator(file, output_gene_names=True)


class TestPSPKinasesAnnotator:
    """Test the PSPKinasesAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file into `psp_dict`.

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
                "PSP Kinases": [
                    "HRI",
                    "PKCD",
                    "PKCD",
                    "PKCD",
                ],
            }
        )
        pd.testing.assert_frame_equal(
            annotator.psp_df, expected_output_df, check_like=True
        )

    def test_load_annotations_output_gene_names(self, annotator_output_gene_names):
        """Test that the load_annotations method correctly loads the annotations from the input file into `psp_dict`.

        Args:
            annotator_output_gene_names: Annotator object with mock file loaded with output_gene_names set to True
        """
        annotator_output_gene_names.load_annotations()

        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "P05198_S52",
                    "Q9UQL6_S259",
                    "P18433-2_S204",
                    "P10415_S70",
                ],
                "PSP Kinases": [
                    "EIF2AK1",
                    "PRKCD",
                    "PRKCD",
                    "PRKCD",
                ],
            }
        )
        pd.testing.assert_frame_equal(
            annotator_output_gene_names.psp_df, expected_output_df, check_like=True
        )

    def test_annotate(self, annotator):
        """Test that the annotate method correctly adds the PSP annotations to the input dataframe as new columns.

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
            "PSP Kinases",
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
                "PSP Kinases": [
                    "HRI;PKCD",
                    "PKCD",
                    "PKCD",
                    "",
                ],
            }
        )
        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)
    
    def test_raise_error_on_missing_input_columns(self, annotator):
        """Test that a value error is raised if the "Site positions" column is missing.

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
