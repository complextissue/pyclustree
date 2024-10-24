# pyclustree

[![Version](https://img.shields.io/pypi/v/pyclustree)](https://pypi.org/project/pyclustree/)
[![License](https://img.shields.io/pypi/l/pyclustree)](https://img.shields.io/pypi/l/pyclustree)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/complextissue/pyclustree/test.yaml)
[![Documentation Status](https://readthedocs.org/projects/pyclustree/badge/?version=stable)](https://pyclustree.readthedocs.io/stable/?badge=stable)
[![Codecov](https://codecov.io/gh/complextissue/pyclustree/graph/badge.svg?token=45BNU20CBP)](https://codecov.io/gh/complextissue/pyclustree)
[![Python Version Required](https://img.shields.io/pypi/pyversions/pyclustree)](https://pypi.org/project/pyclustree/)
[![DOI](https://zenodo.org/badge/857752929.svg)](https://doi.org/10.5281/zenodo.13987570)

Visualize cluster assignments at different resolutions. Possbile applications include finding the optimal resolution for
single-cell RNA-sequencing clusterings.

`pyclustree` is inspired by the R package `clustree`, however, while we aim to provide the same functionality, the API
will differ between the implementations.

## Getting started

Please refer to the [documentation][link-docs].

## Installation

You need to have Python 3.9 or newer installed on your system. If you don't have
Python installed, we recommend installing [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

There are several alternative options to install pyclustree:

1. Install the latest release of `pyclustree` from [PyPI][link-pypi]:

```bash
pip install pyclustree
```

2. Install the latest development version:

```bash
pip install git+https://github.com/complextissue/pyclustree.git@dev
```

## Contact

If you found a bug, please use the [issue tracker][issue-tracker].

## Authors

@maltekuehl
@harryhaller001

Unaffiliated with the creators of the R package `clustree`.

## License

Please refer to the [LICENSE][license] file.

## Citation

Please cite both the original R package as well as this implementation when using `pyclustree`. For example: Cluster resolution was determined based on visualization with pyclustree (Kuehl et al., 2024), a Python implementation of clustree (Zappia et al., 2018).

-   pyclustree: Kuehl, M., Hellmig, M., & Puelles, V. G. (2024). pyclustree: Visualizing cluster resolution optimization for biomedical data (0.3.1). Zenodo. https://doi.org/10.5281/zenodo.13987570
-   clustree: Zappia, L., & Oshlack, A. (2018). Clustering trees: a visualization for evaluating clusterings at multiple resolutions. Gigascience, 7(7), giy083.

[license]: https://github.com/complextissue/pyclustree/blob/main/LICENSE
[issue-tracker]: https://github.com/complextissue/pyclustree/issues
[changelog]: https://pyclustree.readthedocs.io/latest/changelog.html
[link-docs]: https://pyclustree.readthedocs.io
[link-pypi]: https://pypi.org/project/pyclustree
