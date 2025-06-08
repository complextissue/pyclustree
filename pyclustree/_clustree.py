from collections.abc import Callable
from logging import warning

import networkx as nx
import numpy as np
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.colors import Colormap
from numpy.typing import NDArray

from ._utils import order_unique_clusters, transition_matrix


def clustree(
    adata: AnnData,
    cluster_keys: list[str],
    title: str | None = None,
    scatter_reference: str | None = None,
    node_colormap: list[Colormap] | Colormap | str = "tab20",
    node_color_gene: str | None = None,
    node_color_gene_use_raw: bool = True,
    node_color_gene_transformer: Callable | None = None,
    node_size_range: tuple[float, float] = (100, 1000),
    edge_width_range: tuple[float, float] = (0.5, 5.0),
    edge_weight_threshold: float = 0.0,
    x_spacing: float = 2.5,
    y_spacing: float = 1.25,
    order_clusters: bool = True,
    show_colorbar: bool = False,
    show_fraction: bool = False,
    show_cluster_keys: bool = True,
    graph_plot_kwargs: dict | None = None,
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
        transition_matrix(
            df_cluster_assignments[cluster_keys[i]],
            df_cluster_assignments[cluster_keys[i + 1]],
        )
        for i in range(len(cluster_keys) - 1)
    ]

    def cluster_sort_key(cluster):
        try:
            return int(cluster)
        except (ValueError, TypeError):
            return str(cluster)

    unique_clusters = [np.unique(df_cluster_assignments[key]).tolist() for key in cluster_keys]
    unique_clusters = [sorted(unique_clusters_level, key=cluster_sort_key) for unique_clusters_level in unique_clusters]

    if order_clusters:
        unique_clusters = order_unique_clusters(unique_clusters, transition_matrices)

    # Create the Graph
    G: nx.Graph = nx.DiGraph()

    layers: dict[int, list[str]] = {}

    # Add the nodes and store cluster info directly in the graph
    for i, key in enumerate(cluster_keys):
        for cluster in unique_clusters[i]:
            node_name = f"{key}_{cluster}"
            cluster_cells = df_cluster_assignments[key] == cluster
            node_size = np.sum(cluster_cells)

            # Store info directly in the node
            G.add_node(
                node_name,
                level=i,
                key=key,
                cluster=cluster,
                size=node_size,
                cells=cluster_cells,
            )

            layer_key = int(len(cluster_keys) - i)
            if layer_key not in layers:
                layers[layer_key] = []
            layers[layer_key].append(node_name)

    # Compute node sizes scaled to the desired range
    max_size = max(G.nodes[node]["size"] for node in G.nodes)
    for node in G.nodes:
        G.nodes[node]["display_size"] = np.clip(
            (G.nodes[node]["size"] * node_size_range[1]) / max_size,
            *node_size_range,
        )

    # Add edges between each level and the next level
    for i, transition in enumerate(transition_matrices):
        for parent_cluster in transition.index:
            parent_node = f"{cluster_keys[i]}_{parent_cluster}"

            for child_cluster in transition.columns:
                child_node = f"{cluster_keys[i + 1]}_{child_cluster}"
                weight = transition.loc[parent_cluster, child_cluster]

                if weight > edge_weight_threshold:
                    G.add_edge(parent_node, child_node, weight=weight)

    # Prepare node colors
    if isinstance(node_colormap, list):
        # Multiple colormaps for different levels
        assert len(node_colormap) == len(
            cluster_keys
        ), "The length of the colormap list should match the number of cluster keys."
        assert (
            node_color_gene is None
        ), "The node_color_gene argument is not supported when providing a list of colormaps."

        # Apply appropriate colormap to each node based on its level
        for node in G.nodes:
            level = G.nodes[node]["level"]
            cluster = G.nodes[node]["cluster"]
            cmap_level = (
                plt.cm.get_cmap(node_colormap[level]) if isinstance(node_colormap[level], str) else node_colormap[level]
            )
            norm_level = plt.Normalize(vmin=0, vmax=len(unique_clusters[level]) - 1)
            color_idx = unique_clusters[level].index(cluster)
            G.nodes[node]["color"] = cmap_level(norm_level(color_idx))

    elif node_color_gene is not None:
        # Color nodes by gene expression
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

        # Flatten gene_counts if it's 2D
        gene_counts = gene_counts.toarray().flatten() if hasattr(gene_counts, "toarray") else gene_counts.flatten()

        # Apply transformer to get expression value per cluster
        transformer = np.mean if node_color_gene_transformer is None else node_color_gene_transformer

        # Calculate expression for each node
        gene_values = {}
        for node in G.nodes:
            cells = G.nodes[node]["cells"]
            gene_values[node] = np.nan_to_num(transformer(gene_counts[cells]), nan=0)

        # Create color normalization
        gene_min, gene_max = min(gene_values.values()), max(gene_values.values())
        norm_gene = plt.Normalize(vmin=gene_min, vmax=gene_max)

        # Assign colors
        for node in G.nodes:
            G.nodes[node]["color"] = node_colormap(norm_gene(gene_values[node]))

    else:
        # Single colormap for all nodes based on level
        norm = plt.Normalize(vmin=0, vmax=len(cluster_keys) - 1)
        for node in G.nodes:
            G.nodes[node]["color"] = node_colormap(norm(G.nodes[node]["level"]))

    # Calculate the positions of the nodes
    if scatter_reference is not None:
        # Get the median x and y positions of the scatter reference for each cluster
        x_positions = adata.obsm[scatter_reference][:, 0]
        y_positions = adata.obsm[scatter_reference][:, 1]

        node_positions = {}
        for node in G.nodes:
            cells = G.nodes[node]["cells"]
            node_positions[node] = (np.median(x_positions[cells]), np.median(y_positions[cells]))
    else:
        node_positions = nx.multipartite_layout(
            G,
            align="horizontal",
            subset_key=layers,  # type: ignore
            scale=1.0,
        )

    # Prepare edge widths scaled to desired range
    edge_weights = [G.edges[edge]["weight"] for edge in G.edges]
    if edge_weights:  # Only scale if there are edges
        max_weight = max(edge_weights)
        edge_widths = np.clip(
            (np.array(edge_weights) * edge_width_range[1]) / max_weight,
            *edge_width_range,
        )
    else:
        edge_widths = []

    # Create the plot
    figsize = (x_spacing * len(cluster_keys), y_spacing * len(unique_clusters[0]))
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

    # Prepare graph drawing parameters
    graph_plot_kwargs_base = {
        "G": G,
        "with_labels": True,
        "labels": {node: G.nodes[node]["cluster"] for node in G.nodes},
        "pos": node_positions,
        "node_size": [G.nodes[node]["display_size"] for node in G.nodes],
        "node_color": [G.nodes[node]["color"] for node in G.nodes],
        "edge_color": "black",
        "width": edge_widths,
        "font_size": 8,
        "font_color": "white",
        "font_weight": "bold",
        "edge_vmin": 0,
        "edge_vmax": 1,
        "arrowsize": 10,
        "ax": ax,
    }

    if graph_plot_kwargs is not None:
        graph_plot_kwargs_base.update(graph_plot_kwargs)

    # Draw the graph
    nx.draw(**graph_plot_kwargs_base)

    # Add edge labels if requested
    if show_fraction:
        edge_labels = {edge: f"{weight:.2f}" for edge, weight in nx.get_edge_attributes(G, "weight").items()}
        nx.draw_networkx_edge_labels(
            G,
            node_positions,
            label_pos=0.4,
            edge_labels=edge_labels,
            font_size=6,
            rotate=False,
            alpha=0.8,
            bbox={
                "alpha": 0.8,
                "ec": [1.0, 1.0, 1.0],
                "fc": [1.0, 1.0, 1.0],
                "boxstyle": "round,pad=0.3",
            },
        )

    # Plot the colorbar
    if show_colorbar and not isinstance(node_colormap, list):
        if node_color_gene is not None:
            # For gene expression, use the gene value min/max
            gene_min, gene_max = min(gene_values.values()), max(gene_values.values())
            norm_for_colorbar = plt.Normalize(vmin=gene_min, vmax=gene_max)
        else:
            # For clustering levels
            norm_for_colorbar = plt.Normalize(vmin=0, vmax=len(cluster_keys) - 1)

        sm = plt.cm.ScalarMappable(cmap=node_colormap, norm=norm_for_colorbar)
        sm.set_array([])
        label = f"{node_color_gene} expression" if node_color_gene is not None else "Cluster color"
        colorbar = fig.colorbar(
            sm,
            ax=ax,
            orientation="vertical",
            fraction=0.02,
            pad=0.02,
            label=label,
            aspect=10,
        )
        colorbar.ax.yaxis.set_label_position("left")
    elif show_colorbar and isinstance(node_colormap, list):
        warning("Colorbars are not supported when providing a list of colormaps. Ignoring the argument.")

    # Calculate positions for cluster keys
    need_level_positions = show_cluster_keys and scatter_reference is not None
    y_positions_levels: NDArray[np.floating] | list[float]
    if need_level_positions:
        y_positions_levels = np.linspace(
            ax.get_ylim()[1],
            ax.get_ylim()[1] - len(cluster_keys) * 1.0,
            len(cluster_keys),
        )

    if show_cluster_keys:
        x_min = ax.get_xlim()[0] if scatter_reference is None else ax.get_xlim()[1] + 2

        # Determine y positions and colors
        facecolor: list[str] | list[tuple[float, float, float, float]]
        if scatter_reference is None:
            # Use node positions for y-coordinates in hierarchical layout
            level_nodes: dict[int, list[str]] = {}
            for node in G.nodes:
                level = G.nodes[node]["level"]
                if level not in level_nodes:
                    level_nodes[level] = []
                level_nodes[level].append(node)

            y_positions_levels = [
                max(node_positions[node][1] for node in level_nodes.get(i, [])) for i in range(len(cluster_keys))
            ]
            facecolor = ["white"] * len(cluster_keys)
        else:
            # Use previously calculated y positions for scatter layout
            if isinstance(node_colormap, list):
                warning("Cannot show colored cluster keys when providing a list of colormaps. Showing white keys.")
                facecolor = ["white"] * len(cluster_keys)
            else:
                norm = plt.Normalize(vmin=0, vmax=len(cluster_keys) - 1)
                facecolor = [node_colormap(norm(i)) for i in range(len(cluster_keys))]

        # Add text labels for cluster keys
        for i, key in enumerate(cluster_keys):
            ax.text(
                x_min - 0.5,
                y_positions_levels[i],
                key,
                fontsize=12,
                color="black",
                ha="center",
                va="center",
                bbox={
                    "boxstyle": "round",
                    "facecolor": facecolor[i],
                    "edgecolor": "black",
                },
            )

    # Set the title of the plot
    if title is not None:
        ax.set_title(title, fontsize=16, fontweight="bold", pad=10)

    return fig
