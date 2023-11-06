import io

import pandas as pd
import pytest

from psite_annotation.annotators.psp_regulatory import PSPRegulatoryAnnotator

# Define the mock input file as a string
mock_input_file = """110220
PhosphoSitePlus(R) (PSP)

ACC_ID	MOD_RSD	ORGANISM	ON_FUNCTION	ON_PROCESS	ON_PROT_INTERACT	ON_OTHER_INTERACT	NOTES
O75822	S11	human	regulation of transcription; protein degradation; intracellular localization; phosphorylation	cell cycle	Protein Kinase C	regulation of gene expression
O75823	S20	human	regulation of cell growth; phosphorylation	cell cycle	Protein Kinase A	regulation of cell division
O75824	T100	human	regulation of cell signaling; phosphorylation	cell signaling	Phosphatidylinositol 3-kinase	regulation of intracellular signaling cascades
O75823	S20	rat	regulation of cell signaling; phosphorylation	cell signaling	Phosphatidylinositol 3-kinase	regulation of intracellular signaling cascades
"""  # noqa: E501,B950


@pytest.fixture
def annotator() -> PSPRegulatoryAnnotator:
    """Fixture to create a PSPRegulatoryAnnotator object using a mock input file.

    Returns:
        PSPRegulatoryAnnotator: Annotator object with mock file loaded
    """
    file = io.StringIO(mock_input_file)
    return PSPRegulatoryAnnotator(file)


class TestPSPRegulatorysAnnotator:
    """Test the PSPRegulatoryAnnotator class."""

    def test_load_annotations(self, annotator):
        """Test that the load_annotations method correctly loads the annotations from the input file into `psp_dict`.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "O75822_S11",
                    "O75823_S20",
                    "O75824_T100",
                ],
                "PSP_ON_FUNCTION": [
                    "regulation of transcription; protein degradation; intracellular localization; phosphorylation",
                    "regulation of cell growth; phosphorylation",
                    "regulation of cell signaling; phosphorylation",
                ],
                "PSP_ON_PROCESS": ["cell cycle", "cell cycle", "cell signaling"],
                "PSP_ON_PROT_INTERACT": [
                    "Protein Kinase C",
                    "Protein Kinase A",
                    "Phosphatidylinositol 3-kinase",
                ],
                "PSP_ON_OTHER_INTERACT": [
                    "regulation of gene expression",
                    "regulation of cell division",
                    "regulation of intracellular signaling cascades",
                ],
                "PSP_NOTES": ["", "", ""],
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
            {"Site positions": ["O75822_S11;O75823_S20", "O75824_T100", "P05067_S100"]}
        )

        annotator.load_annotations()

        # Annotate the input dataframe
        output_df = annotator.annotate(input_df)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Site positions",
            "PSP_ON_FUNCTION",
            "PSP_ON_PROCESS",
            "PSP_ON_PROT_INTERACT",
            "PSP_ON_OTHER_INTERACT",
            "PSP_NOTES",
        }

        # Assert that the output dataframe has the expected values
        expected_output_df = pd.DataFrame(
            {
                "Site positions": [
                    "O75822_S11;O75823_S20",
                    "O75824_T100",
                    "P05067_S100",
                ],
                "PSP_ON_FUNCTION": [
                    "intracellular localization; phosphorylation; protein degradation; "
                    + "regulation of cell growth; regulation of transcription",
                    "phosphorylation; regulation of cell signaling",
                    "",
                ],
                "PSP_ON_PROCESS": ["cell cycle", "cell signaling", ""],
                "PSP_ON_PROT_INTERACT": [
                    "Protein Kinase A; Protein Kinase C",
                    "Phosphatidylinositol 3-kinase",
                    "",
                ],
                "PSP_ON_OTHER_INTERACT": [
                    "regulation of cell division; regulation of gene expression",
                    "regulation of intracellular signaling cascades",
                    "",
                ],
                "PSP_NOTES": ["", "", ""],
            }
        )
        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)
