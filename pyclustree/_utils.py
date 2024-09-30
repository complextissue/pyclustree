import pandas as pd


def calculate_transition_matrix(
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
