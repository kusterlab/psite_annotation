"""Get version from distribution and set copyright."""

from .functional_annotation import *  # noqa: F401,F403

__version__ = "0.0.0"
try:
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
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
