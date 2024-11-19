import collections
import unittest
import unittest.mock

import pandas as pd
import pytest

import psite_annotation.annotators.site_sequence_context as pa
from psite_annotation.annotators.annotator_base import MissingColumnsError

# Define the mock input file as a string
mock_input_file = """>tr|E9PL77|E9PL77_HUMAN Uncharacterized protein
MKGEGCGHRAGSTKDAGARPGLCRMCGRGSESIKIPRSKQSINNQRQQKNLPKNNYKK
>tr|F2Z2P1|F2Z2P1_HUMAN Uncharacterized protein
MNVLVWEDCIAEQAEVLHNDSYGVIIDCSPKGMFSLNCTSQSACHGHTMFSWSEQNGQMVEMIRSMARVPIIWKHGGIVAPQPQMIWP
AVGAKHKDLWKLLMALNKIKIWERIKKHLEGHSRNLDIAKLKEQIFKASQAHLTLMPGTGVLEGAADGLAAINPLK
"""

mock_input_file_psp = """;110220
;Data extracted from PhosphoSitePlus(R), created by Cell Signaling Technology Inc. PhosphoSitePlus is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. Attribution must be given in written, oral and digital presentations to PhosphoSitePlus, www.phosphosite.org. Written documents should additionally cite Hornbeck PV, Kornhauser JM, Tkachev S, Zhang B, Skrzypek E, Murray B, Latham V, Sullivan M (2012) PhosphoSitePlus: a comprehensive resource for investigating the structure and function of experimentally determined post-translational modifications in man and mouse. Nucleic Acids Res. 40, D261<D0>70.; www.phosphosite.org.

>GN:Cbln1|CBLN1|mouse|Q9R171
MLGVVELLLLGTAWLAGPARGQNETEPIVLEGKCLVVCDSNPTSDPTGTALGISVRSGSA
KVAFSAIRSTNHEPSEMSNRTMIIYFDQVLVNIGNNFDSERSTFIAPRKGIYSFNFHVVK
VYNRQTIQVSLMLNGWPVISAFAGDQDVTREAASNGVLIQMEKGDRAYLKLERGNLMGGW
KYSTFSGFLVFPL
>GN:Cox7a2|COX7A2|mouse|P48771
MLRNLLALRQIAQRTISTTSRRHFENKVPEKQKLFQEDNGMPVHLKGGASDALLYRATMA
LTLGGTAYAIYLLAMAAFPKKQN
>GN:RPS20|RPS20|human|P60866
MAFKDTGKTPVEPEVAIHRIRITLTSRNVKSLEKVCADLIRGAKEKNLKVKGPVRMPTKT
LRITTRKTPCGEGSKTWDRFQMRIHKRLIDLHSPSEIVKQITSISIEPGVEVEVTIADA
>GN:UFC1|UFC1|human|Q9Y3C8
MADEATRRVVSEIPVLKTNAGPRDRELWVQRLKEEYQSLIRYVENNKNADNDWFRLESNK
EGTRWFGKCWYIHDLLKYEFDIEFDIPITYPTTAPEIAVPELDGKTAKMYRGGKICLTDH
FKPLWARNVPKFGLAHLMALGLGPWLAVEIPDLIQKGVIQHKEKCNQ
"""  # noqa: E501, B950


@pytest.fixture
def annotator() -> pa.SiteSequenceContextAnnotator:
    """Fixture to create a MotifAnnotator object using a mock input file.

    Returns:
        MotifAnnotator: Annotator object with mock file loaded
    """
    return pa.SiteSequenceContextAnnotator("")


@pytest.fixture
def annotator_psp() -> pa.SiteSequenceContextAnnotator:
    """Fixture to create a MotifAnnotator object using a mock input file.

    Returns:
        MotifAnnotator: Annotator object with mock file loaded
    """
    return pa.SiteSequenceContextAnnotator("", pspInput=True)


@pytest.fixture
def annotator_custom_context() -> pa.SiteSequenceContextAnnotator:
    """Fixture to create a MotifAnnotator object using a mock input file.

    Returns:
        MotifAnnotator: Annotator object with mock file loaded
    """
    return pa.SiteSequenceContextAnnotator("", context_left=5, context_right=4)


@pytest.fixture
def input_df() -> pd.DataFrame:
    """Fixture to create an example input dataframe to the annotate() method.

    Returns:
        pd.DataFrame: Example input dataframe

    """
    return pd.DataFrame(
        {
            "Site positions": [
                "E9PL77_T13",
                "E9PL77_S30;E9PL77_S32",
                "",
                "F2Z2P1_S53",
                "",
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
            "Site positions": [
                "E9PL77_T13",
                "E9PL77_S30;E9PL77_S32",
                "",
                "F2Z2P1_S53",
                "",
            ],
            "Site sequence context": [
                "___MKGEGCGHRAGStKDAGARPGLCRMCGR",
                "DAGARPGLCRMCGRGsESIKIPRSKQSINNQ;GARPGLCRMCGRGSEsIKIPRSKQSINNQRQ",
                "",
                "CTSQSACHGHTMFSWsEQNGQMVEMIRSMAR",
                "",
            ],
        }
    )


@pytest.fixture
def expected_output_df_custom_context() -> pd.DataFrame:
    """Fixture to create an example ouput dataframe to the annotate() method given the input_df fixture.

    Returns:
        pd.DataFrame: Output dataframe corresponding to input_df fixture

    """
    return pd.DataFrame(
        {
            "Site positions": [
                "E9PL77_T13",
                "E9PL77_S30;E9PL77_S32",
                "",
                "F2Z2P1_S53",
                "",
            ],
            "Site sequence context": [
                "HRAGStKDAG",
                "MCGRGsESIK;GRGSEsIKIP",
                "",
                "TMFSWsEQNG",
                "",
            ],
        }
    )


class TestSiteSequenceContextAnnotator:
    """Test the SiteSequenceContextAnnotator class."""

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file),
        create=True,
    )
    def test_load_annotations(self, annotator):
        """Test that load_annotations method loads annotations from the input file into `protein_sequences`.

        Args:
            annotator: Annotator object with mock file loaded
        """
        annotator.load_annotations()

        assert set(annotator.protein_sequences.keys()) == {
            "E9PL77",
            "F2Z2P1",
        }

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file_psp),
        create=True,
    )
    def test_load_annotations_psp(self, annotator_psp):
        """Test that load_annotations method loads annotations from the input file into `protein_sequences`.

        The PSP fasta reader should ignore non-human sequences.

        Args:
            annotator_psp: Annotator object with mock file loaded
        """
        annotator_psp.load_annotations()

        assert set(annotator_psp.protein_sequences.keys()) == {
            "P60866",
            "Q9Y3C8",
        }

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file),
        create=True,
    )
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
            "Site positions",
            "Site sequence context",
        }

        # Assert that the output dataframe has the expected values
        pd.testing.assert_frame_equal(output_df, expected_output_df, check_like=True)

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file),
        create=True,
    )
    def test_annotate_custom_context(
        self, annotator_custom_context, input_df, expected_output_df_custom_context
    ):
        """Test that the annotate method correctly adds the motif annotations to the input dataframe as new column.

        Args:
            annotator_custom_context: Annotator object with mock file loaded
            input_df: Example input dataframe
            expected_output_df_custom_context: Expected output dataframe

        """
        annotator_custom_context.load_annotations()

        # Annotate the input dataframe
        output_df = annotator_custom_context.annotate(input_df)

        # Assert that the output dataframe has the expected columns
        assert set(output_df.columns) == {
            "Site positions",
            "Site sequence context",
        }

        # Assert that the output dataframe has the expected values
        pd.testing.assert_frame_equal(
            output_df, expected_output_df_custom_context, check_like=True
        )

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file),
        create=True,
    )
    def test_annotate_missing_site_positions_column(self, annotator):
        """Test that the annotate method correctly raises an error if Site positions column is missing.

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

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file),
        create=True,
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

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file),
        create=True,
    )
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


class TestGetSiteSequenceContexts:
    """Test the _get_site_sequence_contexts function."""

    def test_get_site_sequence_contexts(self, proteinSequences):
        """Test the _get_site_sequence_contexts function with two isoforms with the identical site sequence context.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        site_position_string = "Q86U42-2_S19;Q86U42_S19"

        assert (
            pa._get_site_sequence_contexts(site_position_string, proteinSequences)
            == "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE;AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE"
        )

    def test_get_site_sequence_contexts_unique(self, proteinSequences):
        """Test the _get_site_sequence_contexts function with duplicate elimination.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        site_position_string = "Q86U42-2_S19;Q86U42_S19"

        assert (
            pa._get_site_sequence_contexts(site_position_string, proteinSequences, return_unique=True)
            == "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE"
        )

    def test_get_site_sequence_contexts_sorted(self, proteinSequences):
        """Test the _get_site_sequence_contexts function with sorting of output

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        site_position_string = "Q86U42_Y46;Q86U42_S19"

        assert (
            pa._get_site_sequence_contexts(site_position_string, proteinSequences, return_sorted=True)
            == "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE;AGGEAGEGAPGGAGDyGNGLESEELEPEELL"
        )

    def test_get_site_sequence_contexts_custom_context(self, proteinSequences):
        """Test the _get_site_sequence_contexts function with two isoforms with the identical site sequence context.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        site_position_string = "Q86U42-2_S19;Q86U42_S19"

        assert (
            pa._get_site_sequence_contexts(
                site_position_string, proteinSequences, context_left=5, context_right=4
            )
            == "AGGRGsGPGR;AGGRGsGPGR"
        )

    def test_get_site_sequence_contexts_all_potential_sites(
        self, proteinSequencesExtraPhospho
    ):
        """Test the _get_site_sequence_contexts function with multiple sites within one peptide.

        Args:
            proteinSequencesExtraPhospho: dictionary of UniProt identifiers to protein sequences

        """
        site_position_string = (
            "Q86U42_S19;Q86U42_T20;Q86U42_Y21;Q86U42-2_S19;Q86U42-2_T20;Q86U42-2_Y21"
        )

        assert (
            pa._get_site_sequence_contexts(
                site_position_string, proteinSequencesExtraPhospho
            )
            == "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG;AAAAAAAGAAGGRGStYGPGRRRHLVPGAGG;AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE;" +
            "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG;AAAAAAAGAAGGRGStYGPGRRRHLVPGAGG;AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE"
        )

    def test_get_site_sequence_context_missing_protein(self):
        """Test the _get_site_sequence_contexts function when the protein cannot be found in the fasta file."""
        proteinSequences = collections.defaultdict(str)
        sitePosString = "Q86U42_S19"

        assert pa._get_site_sequence_contexts(sitePosString, proteinSequences) == ""

    def test_get_site_sequence_context_out_of_bounds(self, proteinSequences):
        """Test the _get_site_sequence_contexts function when the position is larger than the length of the sequence"""
        sitePosString = "Q86U42_S307" # (Q86U42 has length 306)
        assert pa._get_site_sequence_contexts(sitePosString, proteinSequences) == ""



class TestAddModificationToSequenceContext:
    def test_same_position(self):
        result = pa._add_modification_to_sequence_context(
            "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG",
            "Q86U42-2_S19",
            "Q86U42-2_T20",
            15,
        )
        assert result == "AAAAAAAAGAAGGRGstYGPGRRRHLVPGAG"

    def test_different_protein(self):
        result = pa._add_modification_to_sequence_context(
            "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG",
            "Q86U42-2_S19",
            "Q86U42_T20",
            15,
        )
        assert result == "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG"

    def test_before_context(self):
        result = pa._add_modification_to_sequence_context(
            "AGGRGsTYGPG",
            "Q86U42-2_S19",
            "Q86U42-2_T13",
            5,
        )
        assert result == "AGGRGsTYGPG"

    def test_after_context(self):
        result = pa._add_modification_to_sequence_context(
            "AGGRGsTYGPG",
            "Q86U42-2_S19",
            "Q86U42-2_T25",
            5,
        )
        assert result == "AGGRGsTYGPG"

    def test_wrong_amino_acid(self):
        with pytest.raises(ValueError, match="Incorrect modified amino acid at position"):
            pa._add_modification_to_sequence_context(
                "AGGRGsTYGPG",
                "Q86U42-2_S19",
                "Q86U42-2_T24",
                5,
            )


class TestUnpackSitePositionString:
    def test_valid_format(self):
        result = pa._unpack_site_position_string("protein1_S5")
        assert result == ("protein1", 4, "s")

    def test_invalid_format(self):
        with pytest.raises(ValueError, match="Invalid format for site_position_string"):
            pa._unpack_site_position_string("invalid_string")



# You may need to adjust the imports and module names based on your actual module structure.
