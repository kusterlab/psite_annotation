import io

import pandas as pd
import pytest

from psite_annotation.annotators.psp_studies import PSPStudiesAnnotator

# Define the mock input file as a string
mock_input_file = """110220
PhosphoSitePlus(R) (PSP)

ACC_ID	MOD_RSD	GENE	ORGANISM	LT_LIT	MS_LIT	MS_CST
P05067	S9	A2M	human	1	2	3
P05067	S10	A2M	human	4	5	6
P05067	S11	A2M	human	7	8	9
P05067	S12	A2M	human	10	11	12
P05067	S12	A2M	mouse	13	14	15
"""


@pytest.fixture
def annotator() -> PSPStudiesAnnotator:
    """Fixture to create a PSPStudiesAnnotator object using a mock input file.

    Returns:
        PSPStudiesAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return PSPStudiesAnnotator(file)


class TestPSPStudiesAnnotator:
    """Test the PSPStudiesAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file into `psp_dict`.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "P05067_S9",
                    "P05067_S10",
                    "P05067_S11",
                    "P05067_S12",
                ],
                "PSP_LT_LIT": [1, 4, 7, 10],
                "PSP_MS_LIT": [2, 5, 8, 11],
                "PSP_MS_CST": [3, 6, 9, 12],
            }
        )
        pd.testing.assert_frame_equal(
            annotator.psp_df, expected_output_df, check_like=True
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
                    "P05067_S9;P05067_S10",
                    "P05067_S11",
                    "P05067_S12",
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
            "PSP_LT_LIT",
            "PSP_MS_LIT",
            "PSP_MS_CST",
        }

        # Assert that the output dataframe has the expected values
        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "P05067_S9;P05067_S10",
                    "P05067_S11",
                    "P05067_S12",
                    "P05067_S100",
                ],
                "PSP_LT_LIT": ["1;4", "7", "10", "0"],
                "PSP_MS_LIT": ["2;5", "8", "11", "0"],
                "PSP_MS_CST": ["3;6", "9", "12", "0"],
            }
        )
        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)
