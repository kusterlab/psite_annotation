Basic usage
===========

The psite-annotation package consists of 
`annotator functions <./_autosummary/psite_annotation.html>`_ and 
`annotator classes <./_autosummary/psite_annotation.annotators.html>`_.
Using the `annotator functions <./_autosummary/psite_annotation.html>`_ is 
generally easier, whereas using the 
`annotator classes <./_autosummary/psite_annotation.annotators.html>`_ offers 
more flexibility.

The first step is always to import the package:

.. code-block:: python

    import psite_annotation as pa

After loading your dataframe with pandas, you can then add annotations using
one of the `annotator functions <./_autosummary/psite_annotation.html>`_. These
functions always have the same structure, following general pandas principles:

.. code-block:: python

    df = pa.addSomeAnnotation(df, other_arguments, optional_argument=optional_argument)

This adds one or more new column(s) to the dataframe with the annotation(s).


Annotator functions
-------------------

.. autosummary::
    :template: custom-module-template.rst

    psite_annotation.addPeptideAndPsitePositions
    psite_annotation.addSiteSequenceContext
    psite_annotation.addPSPAnnotations
    psite_annotation.addPSPKinaseSubstrateAnnotations
    psite_annotation.addPSPRegulatoryAnnotations
    psite_annotation.addDomains
    psite_annotation.addTurnoverRates
    psite_annotation.addInVitroKinases
    psite_annotation.addMotifs
    psite_annotation.addKinaseLibraryAnnotations

Please note the following:

- Each annotator function has one or more `required columns`, which are listed in the documentation of the corresponding function. 
- Multiple annotator functions can (and some times have to) be applied to the dataframe in succession. 
- The documentation of each annotator function also includes one or more examples.


Example: Add upstream kinases
-----------------------------

To add upstream kinases to a pandas dataframe :code:`df` with columns :code:`Proteins` 
(UniProt identifiers separated by semicolons, e.g. :code:`Q86U42-2;Q86U42`) and 
:code:`Modified sequence` (standard MaxQuant notation, e.g. :code:`(ac)AAAAAAAAAAGAAGGRGS(ph)GPGR`):

.. code-block:: python

    import psite_annotation as pa

    df = pa.addPeptideAndPsitePositions(df, pa.pspFastaFile, pspInput = True)
    df = pa.addPSPKinaseSubstrateAnnotations(df, pa.pspKinaseSubstrateFile)

We first annotate the peptide and modification positions within the protein 
using :code:`addPeptideAndPsitePositions()`. This adds a column with an identifier 
for the phosphosite, which can then be mapped to the phosphorylating kinase 
using :code:`addPSPKinaseSubstrateAnnotations()`.