# flake8: noqa

import os

import pandas as pd
import pytest

import psite_annotation.functional_annotation as pa

basePath = "tests/data"
curveFile = "curves.txt"

fastaFile = os.path.join(basePath, "uniprot_proteins.fasta")


@pytest.fixture
def curves_df():
    return pd.read_csv(os.path.join(basePath, curveFile), sep="\t")


class TestAddPsitePositions:
    def test_addPeptideAndPsitePositions(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(curves_df, fastaFile)
        # print(curves_df.head())

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site positions"
            ].values[0]
            == "P62877_T9"
        )

    def test_addPeptideAndPsitePositionsPSP(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, pa.pspFastaFile, pspInput=True
        )
        # print(curves_df.head())

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site positions"
            ].values[0]
            == "P62877_T9"
        )


class TestAddPSPAnnotations:
    def test_addPSPAnnotations(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, pa.pspFastaFile, pspInput=True
        )
        curves_df = pa.addPSPAnnotations(curves_df, pa.pspAnnotationFile)
        print(list(curves_df.columns))

        assert (
            curves_df[
                curves_df["Modified sequence"] == "(ac)AAAAAAAGDS(ph)DSWDADAFSVEDPVRK"
            ]["PSP_MS_LIT"].values[0]
            == "22"
        )

    def test_addPSPRegulatoryAnnotations(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, pa.pspFastaFile, pspInput=True
        )
        curves_df = pa.addPSPRegulatoryAnnotations(curves_df, pa.pspRegulatoryFile)
        # print(curves_df[curves_df['PSP_ON_FUNCTION'] != ""].head(n = 100))

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)ADFDDRVS(ph)DEEKVR"][
                "PSP_ON_FUNCTION"
            ].values[0]
            == "molecular association, regulation; phosphorylation"
        )

    def test_addPSPKinaseSubstrateAnnotations(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, pa.pspFastaFile, pspInput=True
        )
        curves_df = pa.addPSPKinaseSubstrateAnnotations(
            curves_df, pa.pspKinaseSubstrateFile
        )
        # print(curves_df[curves_df['PSP Kinases'] != ""].head(n = 100))

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AESSESFTMASS(ph)PAQR"][
                "PSP Kinases"
            ].values[0]
            == "CDK2;CDK7;CK2A1"
        )

    def test_addPSPKinaseSubstrateAnnotationsReturnAllPotentialSites(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, pa.pspFastaFile, pspInput=True, returnAllPotentialSites=True
        )
        curves_df = pa.addPSPKinaseSubstrateAnnotations(
            curves_df, pa.pspKinaseSubstrateFile
        )
        # print(curves_df[curves_df['PSP Kinases'] != ""].head(n = 100))

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AESSESFTMASS(ph)PAQR"][
                "PSP Kinases"
            ].values[0]
            == "CDC7;CDK2;CDK7;CK2A1"
        )


class TestAddDomains:
    def test_addDomains(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(curves_df, fastaFile)
        curves_df = pa.addDomains(curves_df, pa.domainMappingFile)
        # print(curves_df.head())

        assert (
            curves_df[
                curves_df["Modified sequence"] == "(ac)AAAAAAAGDS(ph)DSWDADAFSVEDPVRK"
            ]["Domains"].values[0]
            == "Pfam:eIF3_subunit;low_complexity_region"
        )


class TestAddInVitroKinases:
    def test_addInVitroKinases(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(curves_df, fastaFile)
        curves_df = pa.addInVitroKinases(
            curves_df, pa.inVitroKinaseSubstrateMappingFile
        )
        # print(curves_df[curves_df['In Vitro Kinases'] != ""].head(n = 100))

        assert (
            curves_df[
                curves_df["Modified sequence"]
                == "(ac)ANQVNGNAVQLKEEEEPMDTSS(ph)VTHTEHYK"
            ]["In Vitro Kinases"].values[0]
            == "CK1g1"
        )


class TestAddMotifs:
    def test_addMotifs(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(curves_df, fastaFile)
        curves_df = pa.addMotifs(curves_df, pa.motifsFile)
        # print(curves_df[curves_df['Motifs'] != ""].head(n = 100))

        assert (
            curves_df[
                curves_df["Modified sequence"] == "(ac)AAAAAAAGDS(ph)DSWDADAFSVEDPVRK"
            ]["Motifs"].values[0]
            == "HPRD_CDK2_2;HPRD_bARK_1"
        )


class TestAddTurnoverRates:
    def test_addTurnoverRates(self, curves_df):
        curves_df = pa.addTurnoverRates(curves_df, pa.turnoverFile)
        # print(curves_df[curves_df['PTM_Turnover'] != ""])

        assert (
            curves_df[
                curves_df["Modified sequence"] == "(ac)AAAAAAAGDS(ph)DSWDADAFSVEDPVRK"
            ]["PTM_Turnover"].values[0]
            == "slower"
        )


if __name__ == "__main__":
    unittest.main()
