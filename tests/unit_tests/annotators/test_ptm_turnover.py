import io

import pandas as pd
import pytest
from psite_annotation.annotators.annotator_base import MissingColumnsError

from psite_annotation.annotators.ptm_turnover import PTMTurnoverAnnotator


@pytest.fixture
def annotator() -> PTMTurnoverAnnotator:
    """Fixture to create a PTMTurnoverAnnotator object using a mock input file.

    Returns:
        PTMTurnoverAnnotator: Annotator object with mock file loaded
    """
    mock_file_content = """Significant global,Peptidoforms
FALSE,PEPTIDE1;PEPTIDE2;PEPTIDE3
slower,PEPTIDE4;PEPTIDE5
faster,PEPTIDE6;PEPTIDE7;PEPTIDE8
,PEPTIDE9;PEPTIDE10
"""
    file = io.StringIO(mock_file_content)
    return PTMTurnoverAnnotator(file)


class TestPSPStudiesAnnotator:
    def test_load_annotations(self, annotator: PTMTurnoverAnnotator):
        """Test that the annotations are loaded correctly in the annotator object.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()
        assert annotator.turnover_df.shape == (10, 1)
        assert annotator.turnover_df.loc["PEPTIDE1", "PTM_Turnover"] == "no difference"
        assert annotator.turnover_df.loc["PEPTIDE4", "PTM_Turnover"] == "slower"
        assert annotator.turnover_df.loc["PEPTIDE6", "PTM_Turnover"] == "faster"
        assert (
            annotator.turnover_df.loc["PEPTIDE9", "PTM_Turnover"] == "not enough replicates"
        )


    def test_annotate(self, annotator: PTMTurnoverAnnotator):
        """Test that the annotate method correctly annotates the input dataframe.

        Args:
            annotator: Annotator object with mock file loaded
        """
        df = pd.DataFrame(
            {"Modified sequence": ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3", "PEPTIDE4"]}
        )
        annotator.load_annotations()
        annotated_df = annotator.annotate(df)
        assert annotated_df.shape == (4, 2)
        assert (
            annotated_df.loc[
                annotated_df["Modified sequence"] == "PEPTIDE1", "PTM_Turnover"
            ].values[0]
            == "no difference"
        )
        assert (
            annotated_df.loc[
                annotated_df["Modified sequence"] == "PEPTIDE4", "PTM_Turnover"
            ].values[0]
            == "slower"
        )
    
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
