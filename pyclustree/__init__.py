from importlib.metadata import version

from ._clustree import clustree

pyclustree = clustree

__all__ = ["clustree", "pyclustree"]

__version__ = version("pyclustree")
