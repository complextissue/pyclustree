import matplotlib.pyplot as plt
import scanpy as sc

import pyclustree
from pyclustree import clustree


def test_package_has_version():
    assert pyclustree.__version__ is not None


def test_clustree():
    adata = sc.datasets.pbmc3k_processed()

    # Run leiden clustering for different resolutions
    for resolution in [0.2, 0.4, 0.6, 0.8, 1.0]:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=f"leiden_{str(resolution).replace('.', '_')}",
        )

    # Create a clustree visualization
    fig = clustree(
        adata,
        [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4, 0.6, 0.8, 1.0]],
    )

    assert isinstance(fig, plt.Figure), "pyclustree should return a matplotlib Figure object."
