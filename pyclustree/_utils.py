from collections.abc import Callable
from typing import Literal

import pandas as pd
from anndata import AnnData
from numpy.typing import ArrayLike, NDArray


def transition_matrix(
    cluster_a: pd.Series,
    cluster_b: pd.Series,
) -> pd.DataFrame:
    """Calculate the transition matrix between two clusterings.

    Args:
        cluster_a (pd.Series): The first clustering.
        cluster_b (pd.Series): The second clustering.

    Returns:
        pd.DataFrame: The transition matrix between the two clusterings.
    """
    # Crosstab for counting transitions
    transition_matrix = pd.crosstab(cluster_a, cluster_b, normalize="index")
    return transition_matrix


def get_centered_positions(
    nodes: list[str],
    level: int,
    y_spacing: int = 1,
    x_spacing: int = 1,
) -> dict[str, tuple[float, int]]:
    """Center nodes horizontally and assign y positions based on level.

    Args:
        nodes (list[str]): List of node identifiers.
        level (int): The vertical level at which to place the nodes.
        y_spacing (int, optional): The vertical spacing between levels. Defaults to 1.
        x_spacing (int, optional): The horizontal spacing between nodes. Defaults to 1.

    Returns:
        dict[str, tuple[float, int]]: A dictionary mapping each node to its (x, y) position.
    """
    n_nodes = len(nodes)
    x_start = -(n_nodes - 1) / 2 * x_spacing  # Center the nodes on the x-axis
    pos = {}
    for i, node in enumerate(nodes):
        pos[node] = (
            x_start + i * x_spacing,
            -level * y_spacing,
        )  # Calculate and assign position
    return pos


def order_unique_clusters(
    unique_clusters: list[list[str]],
    transition_matrices: list[pd.DataFrame],
) -> list[list[str]]:
    """Order unique clusters by parent-child relationships.

    Args:
        unique_clusters (list[list[str]]): List of unique clusters.
        transition_matrices (list[pd.DataFrame]): List of transition matrices.

    Returns:
        list[list[str]]: List of unique clusters ordered by parent-child relationships.
    """
    # Initialize with the first hierarchy level, assuming it's already ordered
    ordered_clusters = [unique_clusters[0]]

    # Iterate over each transition matrix (between hierarchy levels)
    for i, transition_matrix in enumerate(transition_matrices):
        ordered_child_clusters = []

        # Get the current level of parent clusters (already ordered)
        parent_clusters = ordered_clusters[i]

        # For each parent cluster, find and order its child clusters
        for parent_cluster in parent_clusters:
            # Find the child clusters where this parent is the maximum contributor
            # idxmax() gives the row (parent cluster) contributing most to each column (child cluster)
            child_max_contributor = transition_matrix.idxmax(axis=0) == parent_cluster

            # Filter relevant child clusters
            relevant_child_clusters = transition_matrix.columns[child_max_contributor]

            if len(relevant_child_clusters) == 0:
                continue  # Skip if no relevant child clusters are found

            # Get the transition fractions for these child clusters
            transition_fractions = transition_matrix.loc[parent_cluster, relevant_child_clusters]

            # Sort relevant child clusters based on the transition fractions in descending order
            sorted_child_clusters = transition_fractions.sort_values(ascending=False).index

            # Add the ordered child clusters for this parent
            ordered_child_clusters.extend(sorted_child_clusters)

        # Add the ordered child clusters to the hierarchy
        ordered_clusters.append(list(ordered_child_clusters))

    return ordered_clusters


def calculate_clustering_score(
    adata: AnnData,
    cluster_key: str,
    score_method: Literal["silhouette", "davies_bouldin", "calinski_harabasz"]
    | Callable[[ArrayLike, ArrayLike], float],
    score_basis: Literal["X", "raw", "pca"] = "pca",
) -> float:
    """Calculate clustering score using specified method and data basis.

    Args:
        adata (AnnData): Annotated data matrix.
        cluster_key (str): Key in `adata.obs` for the clustering labels.
        score_method (Union[Literal["silhouette", "davies_bouldin", "calinski_harabasz"], Callable]): Method to calculate the score.
            Can be a string or a callable function.
        score_basis (Literal["X", "raw", "pca"], optional): Basis for scoring. Defaults to "pca".
            Can be "X" (raw data), "raw" (raw data if available), or "pca" (PCA representation).

    Returns:
        float: The calculated clustering score.

    Raises:
        ValueError: If an invalid score method or basis is provided.
    """
    # Assign basis for scoring
    basis: NDArray | None = None

    if score_basis == "X":
        basis = adata.X.copy()
    elif score_basis == "pca":
        basis = adata.obsm["X_pca"].copy()
    elif score_basis == "raw":
        basis = adata.raw.X.copy()

    assert basis is not None, f"Can't score clustering on basis '{score_basis}'"

    # Calculate score using the appropriate method
    if isinstance(score_method, str):
        if score_method == "calinski_harabasz":
            from sklearn.metrics import calinski_harabasz_score

            return float(calinski_harabasz_score(X=basis, labels=adata.obs[cluster_key]))
        elif score_method == "davies_bouldin":
            from sklearn.metrics import davies_bouldin_score

            return float(davies_bouldin_score(X=basis, labels=adata.obs[cluster_key]))
        elif score_method == "silhouette":
            from sklearn.metrics import silhouette_score

            return float(silhouette_score(X=basis, labels=adata.obs[cluster_key]))
        else:
            raise ValueError(f"Score '{score_method}' not a valid scoring method")
    elif callable(score_method):
        return float(score_method(basis, adata.obs[cluster_key]))

    raise ValueError("Invalid score_method provided")
