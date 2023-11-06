# flake8: noqa

"""Import all annotators so they can be imported from the psite_annotation.annotators subpackage.

Typical usage example:
    from psite_annotation import annotators

    annotator = annotators.PSPKinasesAnnotator(<path_to_annotation_file>)
    annotator.load_annotations()
    df = annotator.annotate(df)
"""

from .clinical_basket import ClinicalBasketAnnotator
from .domain import DomainAnnotator
from .in_vitro_kinases import InVitroKinasesAnnotator
from .motif import MotifAnnotator
from .peptide_position import PeptidePositionAnnotator
from .psp_kinases import PSPKinasesAnnotator
from .psp_regulatory import PSPRegulatoryAnnotator
from .psp_studies import PSPStudiesAnnotator
from .ptm_turnover import PTMTurnoverAnnotator
from .site_sequence_context import SiteSequenceContextAnnotator
