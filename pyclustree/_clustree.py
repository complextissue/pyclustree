from typing import Optional, Union

import networkx as nx
import numpy as np
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.colors import Colormap

from ._utils import calculate_transition_matrix, get_centered_positions


def clustree(
    adata: AnnData,
    cluster_keys: list[str],
    title: Optional[str] = None,
    node_colormap: Union[Colormap, str] = "tab20",
    node_color_gene: Optional[str] = None,
    node_color_gene_use_raw: bool = True,
    node_color_gene_transformer: Optional[callable] = None,
    node_size_range: tuple[float, float] = (100, 1000),
    edge_width_range: tuple[float, float] = (0.5, 5.0),
    edge_weight_threshold: float = 0.0,
    x_spacing: float = 2.5,
    y_spacing: float = 1.25,
    order_clusters: bool = True,
    show_colorbar: bool = False,
    show_cluster_keys: bool = True,
    graph_plot_kwargs: Optional[dict] = None,
) -> plt.Figure:
    """Create a hierarchical clustering tree visualization to compare different clustering resolutions.

    .. code-block:: python

        import scanpy as sc
        from pyclustree import clustree

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

    Args:
        adata (AnnData): The AnnData object from `scanpy` or any other library.
        cluster_keys (list[str]): The list of cluster keys to visualize, in hierarchical order. Keys should be present
            in `adata.obs`.
        title (str, optional): The title of the plot. Defaults to None.
        node_colormap (Union[Colormap, str], optional): The colormap to use for coloring the nodes. Defaults to "tab20".
        node_color_gene (str, optional): The gene to use for coloring the nodes. If provided, node colors will be based
            on the expression of this gene. If None, node colors will be based on the cluster key/level.
            Defaults to None.
        node_color_gene_use_raw (bool, optional): Whether to use the raw data for the gene expression if available.
            Defaults to True.
        node_color_gene_transformer (Optional[callable], optional): A function to transform the gene expression values
            to a single value for coloring the nodes. If None, the mean expression of the gene will be used.
            Defaults to None.
        node_size_range (tuple[float, float], optional): The range of node sizes to use. Defaults to (100, 1000).
        edge_width_range (tuple[float, float], optional): The range of edge widths to use. Defaults to (0.5, 5.0).
        edge_weight_threshold (float, optional): The threshold for edge weights to include in the visualization.
            Defaults to 0.0.
        x_spacing (float, optional): The horizontal spacing between nodes. Defaults to 2.5.
        y_spacing (float, optional): The vertical spacing between nodes. Defaults to 1.25.
        order_clusters (bool, optional): Whether to order the clusters based on the transition matrix. Defaults to True.
        show_colorbar (bool, optional): Whether to show the colorbar. Defaults to False.
        show_cluster_keys (bool, optional): Whether to show the cluster keys on the left side of the plot.
            Defaults to True.
        graph_plot_kwargs (Optional[dict], optional): Additional keyword arguments to pass to `nx.draw`. Will override
            the default arguments. Defaults to None.

    Returns:
        plt.Figure: _description_
    """
    # Ensure all cluster keys are present in adata.obs
    assert all(key in adata.obs for key in cluster_keys), "All cluster keys should be present in adata.obs."

    if node_color_gene is not None:
        if node_color_gene_use_raw:
            assert (
                node_color_gene in adata.raw.var_names
            ), "The provided gene should be present in the adata.raw.var_names."
        else:
            assert node_color_gene in adata.var_names, "The provided gene should be present in the adata.var_names."

    if isinstance(node_colormap, str):
        node_colormap = plt.get_cmap(node_colormap)

    df_cluster_assignments = adata.obs[cluster_keys]

    transition_matrices = [
        calculate_transition_matrix(
            df_cluster_assignments[cluster_keys[i]], df_cluster_assignments[cluster_keys[i + 1]]
        )
        for i in range(len(cluster_keys) - 1)
    ]

    unique_clusters = [np.unique(df_cluster_assignments[key]).tolist() for key in cluster_keys]

    if order_clusters:
        ordered_clusters = [unique_clusters[0]]

        for i, transition_matrix in enumerate(transition_matrices):
            ordered_child_clusters = []

            for parent_cluster in ordered_clusters[i]:
                # Find the child clusters where the argmax of the transition matrix is the parent cluster
                child_clusters = np.array(unique_clusters[i + 1])[
                    (transition_matrix.idxmax(axis=0).values == parent_cluster)
                ]
                for child_cluster in child_clusters:
                    ordered_child_clusters.append(child_cluster)

            ordered_clusters.append(ordered_child_clusters)

        unique_clusters = ordered_clusters

    # Create the Graph
    G = nx.DiGraph()

    # Add the nodes for each cluster key (unique clusters)
    node_names = []

    for i, key in enumerate(cluster_keys):
        level_nodes = [f"{key}_{cluster_name}" for cluster_name in unique_clusters[i]]
        node_names.append(level_nodes)

    G.add_nodes_from([node for nodes in node_names for node in nodes])

    # Add edges between each level and the next level
    for i, transition_matrix in enumerate(transition_matrices):
        for parent_cluster in transition_matrix.index:
            for child_cluster in transition_matrix.columns:
                weight = transition_matrix.loc[parent_cluster, child_cluster]

                if weight > edge_weight_threshold:
                    G.add_edge(
                        f"{cluster_keys[i]}_{parent_cluster}",
                        f"{cluster_keys[i + 1]}_{child_cluster}",
                        weight=weight,
                    )

    # Calculate the positions of the nodes
    node_positions = {}
    for i in range(len(cluster_keys)):
        node_positions.update(
            get_centered_positions(
                nodes=node_names[i],
                level=i,
                y_spacing=y_spacing,
                x_spacing=x_spacing,
            )
        )

    # Use the provided colormap to color the nodes
    node_levels = [[i] * len(unique_clusters[i]) for i in range(len(unique_clusters))]
    node_levels_flat = [item for sublist in node_levels for item in sublist]
    norm = plt.Normalize(vmin=0, vmax=len(unique_clusters[0]))
    node_colors = [node_colormap(norm(value)) for value in node_levels_flat]

    # If a node color gene is provided, use the expression of the gene to color the nodes
    if node_color_gene is not None:
        # Get the expression values for the specified gene
        if node_color_gene_use_raw:
            gene_counts = adata.raw.X[:, adata.raw.var_names == node_color_gene]
        else:
            gene_counts = adata.X[:, adata.var_names == node_color_gene]

        # Flatten gene_counts if it's 2D (i.e., n_cells x 1)
        gene_counts = gene_counts.toarray().flatten() if hasattr(gene_counts, "toarray") else gene_counts.flatten()

        # Compute median expression of the gene for each cluster
        node_color_gene_transformer = np.mean if node_color_gene_transformer is None else node_color_gene_transformer
        gene_cluster_means = []
        for i, key in enumerate(cluster_keys):
            gene_cluster_means.append(
                [
                    np.nan_to_num(
                        node_color_gene_transformer(gene_counts[df_cluster_assignments[key] == cluster]),
                        nan=0,
                    )
                    for cluster in unique_clusters[i]
                ]
            )

        gene_cluster_means_flat = [item for sublist in gene_cluster_means for item in sublist]
        print(gene_cluster_means_flat, len(gene_cluster_means_flat))

        gene_min, gene_max = np.min(gene_cluster_means_flat), np.max(gene_cluster_means_flat)
        norm_gene = plt.Normalize(vmin=gene_min, vmax=gene_max)

        node_colors = [node_colormap(norm_gene(gene_cluster_mean)) for gene_cluster_mean in gene_cluster_means_flat]

    # Scale the node size based on the number of cells in each cluster
    node_sizes = [adata.obs[key].value_counts().sort_index().values for key in cluster_keys]
    node_sizes_flat = [item for sublist in node_sizes for item in sublist]
    node_sizes_flat = np.clip(
        (np.array(node_sizes_flat) * node_size_range[1]) / np.max(node_sizes_flat), *node_size_range
    )

    # Create the plot
    figsize = x_spacing * len(cluster_keys), y_spacing * len(unique_clusters[0])
    fig, ax = plt.subplots(figsize=figsize, dpi=300)

    # Scale the edge width based on the edge weight
    edge_widths = [G.edges[edge]["weight"] for edge in G.edges]

    # Set the maximum edge width to be the maximum and clip the values to be within the range
    edge_widths = np.clip((np.array(edge_widths) * edge_width_range[1]) / np.max(edge_widths), *edge_width_range)

    # Draw the graph
    graph_plot_kwargs = {
        "G": G,
        "with_labels": True,
        "labels": {node: node.split("_")[-1] for node in G.nodes},
        "pos": node_positions,
        "node_color": node_colors,
        "edge_color": "black",
        "width": edge_widths,
        "font_size": 10,
        "font_color": "white",
        "font_weight": "bold",
        "node_size": node_sizes_flat,
        "edge_cmap": plt.cm.Greys,
        "edge_vmin": 0,
        "edge_vmax": 1,
        "arrowsize": 20,
        "ax": ax,
        "connectionstyle": "arc3,rad=0.05",
    }

    if graph_plot_kwargs is not None:
        graph_plot_kwargs.update(graph_plot_kwargs)

    nx.draw(
        **graph_plot_kwargs,
    )

    # Plot the colorbar
    if show_colorbar:
        if node_color_gene is not None:
            sm = plt.cm.ScalarMappable(cmap=node_colormap, norm=norm_gene)
        else:
            sm = plt.cm.ScalarMappable(cmap=node_colormap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.02)

    # Show the name of the cluster key on the left side of the plot
    if show_cluster_keys:
        x_min, _x_max = ax.get_xlim()
        y_positions_levels = [node_positions[node_names[i][0]][1] for i in range(len(node_names))]

        for i, key in enumerate(cluster_keys):
            x = x_min - 0.5
            y = y_positions_levels[i]

            # plot on top of a rounded rectangle in the color of the node
            ax.text(
                x,
                y,
                key,
                fontsize=12,
                color="black",
                ha="center",
                va="center",
                bbox={"boxstyle": "round", "facecolor": "white", "edgecolor": "black"},
            )

    # Set the title of the plot
    if title is not None:
        ax.set_title(title, fontsize=16, fontweight="bold", pad=20)

    return fig
