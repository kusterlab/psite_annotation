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

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site sequence context"
            ].values[0]
            == "_______MAAAMDVDtPSGTNSGAGKKRFEV"
        )

    def test_addPeptideAndPsitePositions_custom_context(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, fastaFile, context_left=5, context_right=7
        )
        # print(curves_df.head())

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site positions"
            ].values[0]
            == "P62877_T9"
        )

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site sequence context"
            ].values[0]
            == "AMDVDtPSGTNSG"
        )

    def test_addPeptideAndPsitePositions_retain_other_mods(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, fastaFile, retain_other_mods=True
        )
        # print(curves_df.head())

        assert (
            curves_df[
                curves_df["Modified sequence"] == "(ac)AAS(ph)DTERDGLAPEKT(ph)SPDRDK"
            ]["Site positions"].values[0]
            == "Q7L4I2_S4;Q7L4I2_T16;E9PI52_S4;E9PI52_T16"
        )

        assert (
            curves_df[
                curves_df["Modified sequence"] == "(ac)AAS(ph)DTERDGLAPEKT(ph)SPDRDK"
            ]["Site sequence context"].values[0]
            == "____________MAAsDTERDGLAPEKtSPD;MAAsDTERDGLAPEKtSPDRDKKKEQSEVSV;____________MAAsDTERDGLAPEKtSPD;MAAsDTERDGLAPEKtSPDRDKKKEQSEVSV"
        )

    def test_addPeptideAndPsitePositions_localization_uncertainty(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, fastaFile, localization_uncertainty=3
        )
        # print(curves_df.head())

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site positions"
            ].values[0]
            == "P62877_T9;P62877_S11"
        )

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site sequence context"
            ].values[0]
            == "_______MAAAMDVDtPSGTNSGAGKKRFEV;_____MAAAMDVDTPsGTNSGAGKKRFEVKK"
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

    def test_addPeptideAndPsitePositions_custom_mod_dict(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, fastaFile, mod_dict={"(ac)A": "a"}
        )
        # print(curves_df.head())

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site positions"
            ].values[0]
            == "P62877_A2"
        )

        assert (
            curves_df[curves_df["Modified sequence"] == "(ac)AAAMDVDT(ph)PSGTNSGAGK"][
                "Site sequence context"
            ].values[0]
            == "______________MaAAMDVDTPSGTNSGA"
        )


class TestAddPSPAnnotations:
    def test_addPSPAnnotations(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(
            curves_df, pa.pspFastaFile, pspInput=True
        )
        curves_df = pa.addPSPAnnotations(curves_df, pa.pspAnnotationFile)
        # print(list(curves_df.columns))

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
        # print(curves_df[
        #         curves_df["Modified sequence"] == "(ac)AAAAAAAGDS(ph)DSWDADAFSVEDPVRK"
        #     ])

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


class TestAddKinaseLibraryAnnotations:
    def test_addKinaseLibraryAnnotations(self, curves_df):
        curves_df = pa.addPeptideAndPsitePositions(curves_df, fastaFile)
        curves_df = pa.addKinaseLibraryAnnotations(
            curves_df,
            pa.kinaseLibraryMotifsFile,
            pa.kinaseLibraryQuantilesFile,
            split_sequences=True,
        )
        # print(curves_df[curves_df["Motif Kinases"] != ""].head(n=100))

        assert (
            curves_df[
                curves_df["Modified sequence"]
                == "(ac)ANQVNGNAVQLKEEEEPMDTSS(ph)VTHTEHYK"
            ]["Motif Kinases"].values[0]
            == "CSNK1A1"
        )
