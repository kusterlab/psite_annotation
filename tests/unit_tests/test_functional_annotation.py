import unittest
import unittest.mock
import pandas as pd

import psite_annotation.functional_annotation as pa


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


class TestAddSiteSequenceContext:
    """Test the addSiteSequenceContext function."""

    @unittest.mock.patch(
        "builtins.open",
        new=unittest.mock.mock_open(read_data=mock_input_file_psp),
        create=True,
    )
    def test_add_site_sequence_context(self):
        df = pd.DataFrame({"Site positions": ["Q9Y3C8_T6", "Q9Y3C8_Q167", "Q9Y3C8_Q168"]})
        result_df = pa.addSiteSequenceContext(df, "mock_input_file_psp", pspInput=True)
        
        assert result_df["Site sequence context"].iloc[0] == "__________MADEAtRRVVSEIPVLKTNAG"
        assert result_df["Site sequence context"].iloc[1] == "DLIQKGVIQHKEKCNq_______________"
        assert result_df["Site sequence context"].iloc[2] == ""
