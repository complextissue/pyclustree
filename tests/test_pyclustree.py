import matplotlib.pyplot as plt
import pytest
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

    fig = clustree(
        adata,
        [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4, 0.6, 0.8, 1.0]],
        node_color_gene="CD8A",
    )

    assert isinstance(fig, plt.Figure), "pyclustree should return a matplotlib Figure object."


def test_scatter_reference():
    adata = sc.datasets.pbmc3k_processed()

    # Run leiden clustering for different resolutions
    for resolution in [0.2, 1.0]:
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
        [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 1.0]],
        title="Clusters projected on UMAP",
        scatter_reference="X_umap",
    )

    assert isinstance(fig, plt.Figure), "pyclustree should return a matplotlib Figure object."

    # Testing errorhandling

    # Testing none existing gene fo node color
    with pytest.raises(AssertionError):
        clustree(
            adata,
            [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 1.0]],
            node_color_gene="Non-existing gene",
        )

    with pytest.raises(AssertionError):
        clustree(
            adata,
            [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 1.0]],
            node_color_gene="Non-existing gene",
            node_color_gene_use_raw=False,
        )

    # Testing node colormap argument
    with pytest.raises(AssertionError):
        cluster_keys = [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 1.0]]
        clustree(adata, cluster_keys, node_colormap=["#FF0000"] * (len(cluster_keys) + 1))

    # Testing node_color_gene when node colormap argument is provided
    with pytest.raises(AssertionError):
        cluster_keys = [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 1.0]]
        clustree(
            adata,
            cluster_keys,
            node_colormap=["#FF0000"] * len(cluster_keys),
            node_color_gene="CD8A",
        )


def test_sankey_basic():
    """Test basic sankey plot functionality."""
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

    # Create a sankey visualization
    fig = clustree(
        adata,
        [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4, 0.6, 0.8, 1.0]],
        transition_plot="sankey",
    )

    assert isinstance(fig, plt.Figure), "clustree with sankey plot should return a matplotlib Figure object."


def test_sankey_with_gene_coloring():
    """Test sankey plot with gene expression coloring."""
    adata = sc.datasets.pbmc3k_processed()

    # Run leiden clustering for different resolutions
    for resolution in [0.2, 0.4, 0.6]:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=f"leiden_{str(resolution).replace('.', '_')}",
        )

    # Create a sankey visualization with gene coloring
    fig = clustree(
        adata,
        [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4, 0.6]],
        transition_plot="sankey",
        node_color_gene="CD8A",
    )

    assert isinstance(fig, plt.Figure), "clustree with sankey and gene coloring should return a matplotlib Figure."


def test_sankey_parameters():
    """Test sankey plot with various parameter combinations."""
    adata = sc.datasets.pbmc3k_processed()

    # Run leiden clustering for different resolutions
    for resolution in [0.2, 0.4, 0.6]:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=f"leiden_{str(resolution).replace('.', '_')}",
        )

    cluster_keys = [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4, 0.6]]

    # Test different curve types
    for curve_type in ["curve3", "curve4", "line"]:
        fig = clustree(
            adata,
            cluster_keys,
            transition_plot="sankey",
            sankey_kwargs={"curve_type": curve_type},
        )
        assert isinstance(fig, plt.Figure), f"Sankey plot with curve_type={curve_type} should work."

    # Test different spacing values
    fig = clustree(
        adata,
        cluster_keys,
        transition_plot="sankey",
        sankey_kwargs={"spacing": 0.05},
    )
    assert isinstance(fig, plt.Figure), "Sankey plot with custom spacing should work."

    # Test different ribbon parameters
    fig = clustree(
        adata,
        cluster_keys,
        transition_plot="sankey",
        sankey_kwargs={"ribbon_alpha": 0.5, "ribbon_color": "blue"},
    )
    assert isinstance(fig, plt.Figure), "Sankey plot with custom ribbon parameters should work."

    # Test different annotation types
    for annotate_type in ["index", "weight", "weight_percent", None]:
        fig = clustree(
            adata,
            cluster_keys,
            transition_plot="sankey",
            sankey_kwargs={"annotate_columns": annotate_type},
        )
        assert isinstance(fig, plt.Figure), f"Sankey plot with annotate_columns={annotate_type} should work."

    # Test with column width parameter
    fig = clustree(
        adata,
        cluster_keys,
        transition_plot="sankey",
        sankey_kwargs={"rel_column_width": 0.25},
    )
    assert isinstance(fig, plt.Figure), "Sankey plot with custom column width should work."

    # Test with legend
    fig = clustree(
        adata,
        cluster_keys,
        transition_plot="sankey",
        sankey_kwargs={"show_legend": True},
    )
    assert isinstance(fig, plt.Figure), "Sankey plot with legend should work."

    # Test with colorbar
    fig = clustree(
        adata,
        cluster_keys,
        transition_plot="sankey",
        sankey_kwargs={"show_colorbar": True},
    )
    assert isinstance(fig, plt.Figure), "Sankey plot with colorbar should work."


def test_sankey_with_scatter_reference_warning():
    """Test that sankey plot ignores scatter_reference with a warning."""
    adata = sc.datasets.pbmc3k_processed()

    # Run leiden clustering for different resolutions
    for resolution in [0.2, 0.4]:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=f"leiden_{str(resolution).replace('.', '_')}",
        )

    # This should produce a warning but still work
    with pytest.warns(UserWarning, match="scatter_reference.*not supported.*sankey"):
        fig = clustree(
            adata,
            [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4]],
            transition_plot="sankey",
            scatter_reference="X_umap",
        )
        assert isinstance(fig, plt.Figure), "Sankey plot should work despite scatter_reference warning."


def test_sankey_with_title():
    """Test sankey plot with title."""
    adata = sc.datasets.pbmc3k_processed()

    # Run leiden clustering for different resolutions
    for resolution in [0.2, 0.4]:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=2,
            key_added=f"leiden_{str(resolution).replace('.', '_')}",
        )

    fig = clustree(
        adata,
        [f"leiden_{str(resolution).replace('.', '_')}" for resolution in [0.2, 0.4]],
        transition_plot="sankey",
        title="Sankey Clustree Visualization",
    )

    assert isinstance(fig, plt.Figure), "Sankey plot with title should work."
