import collections

import psite_annotation.annotators.site_sequence_context as pa


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
            == "AAAAAAAAGAAGGRGsGPGRRRHLVPGAGGE"
        )

    def test_get_site_sequence_contexts_all_potential_sites(
        self, proteinSequencesExtraPhospho
    ):
        """Test the _get_site_sequence_contexts function with multiple sites within one peptide.

        Args:
            proteinSequencesExtraPhospho: dictionary of UniProt identifiers to protein sequences

        """
        site_position_string = (
            "Q86U42-2_S19;Q86U42-2_T20;Q86U42-2_Y21;Q86U42_S19;Q86U42_T20;Q86U42_Y21"
        )

        assert (
            pa._get_site_sequence_contexts(
                site_position_string, proteinSequencesExtraPhospho
            )
            == "AAAAAAAAGAAGGRGsTYGPGRRRHLVPGAG;AAAAAAAGAAGGRGStYGPGRRRHLVPGAGG;AAAAAAGAAGGRGSTyGPGRRRHLVPGAGGE"
        )

    def test_get_site_sequence_context_missing_protein(self):
        """Test the _get_site_sequence_contexts function when the protein cannot be found in the fasta file."""
        proteinSequences = collections.defaultdict(str)
        sitePosString = "Q86U42_S19"

        assert pa._get_site_sequence_contexts(sitePosString, proteinSequences) == ""
