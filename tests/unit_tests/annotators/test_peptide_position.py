import psite_annotation.annotators.peptide_position as pa


class TestGetPeptidePositions:
    """Test the _get_peptide_positions function."""

    def test_get_peptide_positions(self, proteinSequences):
        """Test the _get_peptide_positions function with two isoforms.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        proteinIds = "Q86U42;Q86U42-2"
        modPeptideSequence = "(ac)AAAAAAAAAAGAAGGRGS(ph)GPGR"

        assert pa._get_peptide_positions(
            proteinIds, proteinSequences, modPeptideSequence
        ) == (
            "1;1",
            "23;23",
            "Q86U42-2_S19;Q86U42_S19",
        )

    def test_get_peptide_positions_empty_sequence(self, proteinSequences):
        """Test the _get_peptide_positions function with empty modified sequence.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        proteinIds = "Q86U42;Q86U42-2"
        modPeptideSequence = ""

        assert pa._get_peptide_positions(
            proteinIds, proteinSequences, modPeptideSequence
        ) == (
            "0;0",
            "0;0",
            "",
        )

    def test_get_peptide_positions_empty_proteins(self, proteinSequences):
        """Test the _get_peptide_positions function with empty modified sequence.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        proteinIds = ""
        modPeptideSequence = "(ac)AAAAAAAAAAGAAGGRGS(ph)GPGR"

        assert pa._get_peptide_positions(
            proteinIds, proteinSequences, modPeptideSequence
        ) == (
            "",
            "",
            "",
        )

    def test_get_peptide_positions_peptide_not_in_proteins(self, proteinSequences):
        """Test the _get_peptide_positions function with empty modified sequence.

        Args:
            proteinSequences: dictionary of UniProt identifiers to protein sequences

        """
        proteinIds = "Q86U42;Q86U42-2"
        modPeptideSequence = "(ac)AXXXAAAAAAAAAGAAGGRGS(ph)GPGR"

        assert pa._get_peptide_positions(
            proteinIds, proteinSequences, modPeptideSequence
        ) == (
            "-1;-1",
            "-1;-1",
            "",
        )

    def test_get_peptide_positions_all_potential_sites(
        self, proteinSequencesExtraPhospho
    ):
        """Test the _get_peptide_positions function with two isoforms and the returnAllPotentialSites option.

        Args:
            proteinSequencesExtraPhospho: dictionary of UniProt identifiers to protein sequences

        """
        proteinIds = "Q86U42;Q86U42-2"
        modPeptideSequence = "(ac)AAAAAAAAAAGAAGGRGS(ph)TYGPGR"

        assert pa._get_peptide_positions(
            proteinIds,
            proteinSequencesExtraPhospho,
            modPeptideSequence,
            returnAllPotentialSites=True,
        ) == (
            "1;1",
            "25;25",
            "Q86U42-2_S19;Q86U42-2_T20;Q86U42-2_Y21;Q86U42_S19;Q86U42_T20;Q86U42_Y21",
        )
