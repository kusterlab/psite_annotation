"""Convenience functions for annotating a pandas dataframe with a variety of annotations."""

from .functional_annotation import *  # noqa: F401,F403

"""Get version from distribution and set copyright."""
__version__ = "0.0.0"
try:
    from importlib.metadata import (
        PackageNotFoundError as _PackageNotFoundError,
        version as _version,
    ) # making import aliases private prevents sphinx-autodoc from listing them

    try:
        __version__ = _version(__name__)
    except _PackageNotFoundError:
        pass
except ImportError:
    from pkg_resources import DistributionNotFound, get_distribution

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        pass

__copyright__ = """Copyright (c) 2020-2022 Matthew The. All rights reserved.
Written by Matthew The (matthew.the@tum.de) at the
Chair of Proteomics and Bioanalytics at the Technical University of Munich."""
