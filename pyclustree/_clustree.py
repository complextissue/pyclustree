from logging import warning
from typing import Optional, Union

import networkx as nx
import numpy as np
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.colors import Colormap

from ._utils import calculate_transition_matrix, get_centered_positions, order_unique_clusters


def clustree(
    adata: AnnData,
    cluster_keys: list[str],
    title: Optional[str] = None,
    scatter_reference: Optional[str] = None,
    node_colormap: Union[list[Colormap], Colormap, str] = "tab20",
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
    show_fraction: bool = False,
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
        scatter_reference (str, optional): The key in `adata.obsm` to use as a reference for the scatter plot. If None,
            the nodes will be placed in a hierarchical tree. Defaults to None.
        node_colormap (Union[Colormap, str], optional): The colormap to use for coloring the nodes. If a list is
            provided, the first colormap will be used for the first clustering, the second colormap for the second
            clustering, and so on. For each clustering, the colors will be scaled based on the number of clusters.
            Defaults to "tab20".
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
        show_fraction (bool, optional): Whether to show the fraction of cells from the parent cluster that transitioned
            to the child cluster. Defaults to False.
        show_cluster_keys (bool, optional): Whether to show the cluster keys on the left side of the plot.
            Defaults to True.
        graph_plot_kwargs (Optional[dict], optional): Additional keyword arguments to pass to `nx.draw`. Will override
            the default arguments. Defaults to None.

    Returns:
        plt.Figure: The matplotlib figure object of the clustree visualization.
    """
    # Ensure all cluster keys are present in adata.obs
    assert all(key in adata.obs for key in cluster_keys), "All cluster keys should be present in adata.obs."

    assert (
        scatter_reference is None or scatter_reference in adata.obsm
    ), "The provided scatter reference is not valid. It should be present in adata.obsm."

    if node_color_gene is not None:
        if scatter_reference is not None:
            raise ValueError("Currently, you cannot provide both a scatter reference and a node color gene.")

        if node_color_gene_use_raw and node_color_gene not in adata.obs.columns:
            assert (
                node_color_gene in adata.raw.var_names
            ), "The provided gene should be present in the adata.raw.var_names."
        else:
            assert (
                node_color_gene in adata.obs.columns or node_color_gene in adata.var_names
            ), "The provided gene should be present in the adata.var_names/adata.raw.var_names or adata.obs."

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
    unique_clusters_sorted = [sorted(unique_clusters_level) for unique_clusters_level in unique_clusters]

    if order_clusters:
        unique_clusters = order_unique_clusters(unique_clusters, transition_matrices)

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
        if scatter_reference is not None:
            # Get the median x and y positions of the scatter reference for each cluster
            x_positions = adata.obsm[scatter_reference][:, 0]
            y_positions = adata.obsm[scatter_reference][:, 1]

            cluster_positions = {
                cluster: (
                    np.median(x_positions[df_cluster_assignments[cluster_keys[i]] == cluster]),
                    np.median(y_positions[df_cluster_assignments[cluster_keys[i]] == cluster]),
                )
                for cluster in unique_clusters[i]
            }

            node_positions.update(
                {
                    f"{cluster_keys[i]}_{cluster}": (
                        cluster_positions[cluster][0],
                        cluster_positions[cluster][1],
                    )
                    for j, cluster in enumerate(unique_clusters[i])
                }
            )
        else:
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
    norm = plt.Normalize(vmin=0, vmax=len(np.unique(node_levels_flat)) - 1)

    if not isinstance(node_colormap, list):
        node_colors = [node_colormap(norm(value)) for value in node_levels_flat]
    else:
        # Ensure that the length of the list corresponds to the number of clusterings
        assert len(node_colormap) == len(
            cluster_keys
        ), "The length of the colormap list should match the number of cluster keys."
        assert (
            node_color_gene is None
        ), "The node_color_gene argument is not supported when providing a list of colormaps."

        node_colors = []
        for i in range(len(cluster_keys)):
            node_colormap_level = (
                plt.cm.get_cmap(node_colormap[i]) if isinstance(node_colormap[i], str) else node_colormap[i]
            )
            norm_level = plt.Normalize(vmin=0, vmax=len(unique_clusters[i]) - 1)

            node_colors.extend(
                [
                    node_colormap_level(norm_level(unique_clusters_sorted[i].index(unique_cluster)))
                    for unique_cluster in unique_clusters[i]
                ]
            )

    # If a node color gene is provided, use the expression of the gene to color the nodes
    if node_color_gene is not None:
        # Get the expression values for the specified gene
        if node_color_gene in adata.obs.columns:
            gene_counts = adata.obs[node_color_gene].values
        elif node_color_gene_use_raw:
            if adata.raw is None:
                raise ValueError(
                    "The raw data is not available. Please set `node_color_gene_use_raw` to False or provide raw data."
                )

            gene_counts = adata.raw.X[:, adata.raw.var_names == node_color_gene]
        else:
            gene_counts = adata.X[:, adata.var_names == node_color_gene]

        # Flatten gene_counts if it's 2D (i.e., n_cells x 1)
        gene_counts = gene_counts.toarray().flatten() if hasattr(gene_counts, "toarray") else gene_counts.flatten()

        # Compute median expression of the gene for each cluster
        node_color_gene_transformer = np.mean if node_color_gene_transformer is None else node_color_gene_transformer
        gene_cluster_means = []
        for i, key in enumerate(cluster_keys):
            for cluster in unique_clusters[i]:
                gene_cluster_means.append(
                    np.nan_to_num(
                        node_color_gene_transformer(gene_counts[df_cluster_assignments[key] == cluster]),
                        nan=0,
                    )
                )

        gene_min, gene_max = np.min(gene_cluster_means), np.max(gene_cluster_means)
        norm_gene = plt.Normalize(vmin=gene_min, vmax=gene_max)

        node_colors = [node_colormap(norm_gene(gene_cluster_mean)) for gene_cluster_mean in gene_cluster_means]

    # Scale the node size based on the number of cells in each cluster
    node_sizes = []
    for i in range(len(cluster_keys)):
        for cluster in unique_clusters[i]:
            node_sizes.append(np.sum(df_cluster_assignments[cluster_keys[i]] == cluster))

    node_sizes = np.clip((np.array(node_sizes) * node_size_range[1]) / np.max(node_sizes), *node_size_range)

    # Create the plot
    figsize = x_spacing * len(cluster_keys), y_spacing * len(unique_clusters[0])
    fig, ax = plt.subplots(figsize=figsize, dpi=300)

    # Plot the scatter reference if provided
    if scatter_reference is not None:
        ax.scatter(
            adata.obsm[scatter_reference][:, 0],
            adata.obsm[scatter_reference][:, 1],
            c="lightgrey",
            alpha=0.5,
            s=10,
        )

    # Scale the edge width based on the edge weight
    edge_widths = [G.edges[edge]["weight"] for edge in G.edges]

    # Set the maximum edge width to be the maximum and clip the values to be within the range
    edge_widths = np.clip((np.array(edge_widths) * edge_width_range[1]) / np.max(edge_widths), *edge_width_range)

    # Draw the graph
    graph_plot_kwargs_base = {
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
        "node_size": node_sizes,
        "edge_cmap": plt.cm.Greys,
        "edge_vmin": 0,
        "edge_vmax": 1,
        "arrowsize": 10,
        "ax": ax,
        "connectionstyle": "arc3,rad=0.05",
    }

    if graph_plot_kwargs is not None:
        graph_plot_kwargs_base.update(graph_plot_kwargs)

    nx.draw(
        **graph_plot_kwargs_base,
    )

    # Add the edge labels
    if show_fraction:
        edge_labels = nx.get_edge_attributes(G, "weight")
        formatted_edge_labels = {}
        for edge_label in zip(edge_labels.keys(), edge_labels.values()):
            formatted_edge_labels[edge_label[0]] = f"{np.round(edge_label[1], 2)}"

        nx.draw_networkx_edge_labels(G, node_positions, edge_labels=formatted_edge_labels, font_size=6)

    # Plot the colorbar
    if show_colorbar and not isinstance(node_colormap, list):
        if node_color_gene is not None:
            sm = plt.cm.ScalarMappable(cmap=node_colormap, norm=norm_gene)
        else:
            sm = plt.cm.ScalarMappable(cmap=node_colormap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", fraction=0.02, pad=0.02)
    elif show_colorbar and isinstance(node_colormap, list):
        warning("Colorbars are not supported when providing a list of colormaps. Ignoring the argument.")

    # Show the name of the cluster key on the left side of the plot
    if show_cluster_keys:
        x_min = ax.get_xlim()[0] if scatter_reference is None else ax.get_xlim()[1] + 2
        if scatter_reference is not None:
            # Position them in equal intervals along the y-axis
            y_positions_levels = np.linspace(
                ax.get_ylim()[1], ax.get_ylim()[1] - len(cluster_keys) * 1.0, len(cluster_keys)
            )
            # Use the level colors for the facecolor
            if isinstance(node_colormap, list):
                warning("Cannot show colored cluster keys when providing a list of colormaps. Showing white keys.")
                facecolor = ["white"] * len(cluster_keys)
            else:
                facecolor = [node_colormap(norm(i)) for i in range(len(cluster_keys))]
        else:
            y_positions_levels = [node_positions[node_names[i][0]][1] for i in range(len(node_names))]
            facecolor = ["white"] * len(cluster_keys)

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
                bbox={"boxstyle": "round", "facecolor": facecolor[i], "edgecolor": "black"},
            )

    # Set the title of the plot
    if title is not None:
        ax.set_title(title, fontsize=16, fontweight="bold", pad=20)

    return fig
