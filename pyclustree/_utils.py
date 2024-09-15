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
        pos[node] = (x_start + i * x_spacing, -level * y_spacing)  # Calculate and assign position
    return pos
