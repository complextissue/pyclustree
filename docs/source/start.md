# pyclustree

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/complextissue/pyclustree/test.yaml?branch=main
[link-tests]: https://github.com/complextissue/pyclustree/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/pyclustree

[![Version](https://img.shields.io/pypi/v/pyclustree)](https://pypi.org/project/pyclustree/)
[![License](https://img.shields.io/pypi/l/pyclustree)](https://github.com/complextissue/pyclustree)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/complextissue/pyclustree/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/pyclustree/badge/?version=latest)](https://pyclustree.readthedocs.io/en/latest/?badge=latest)
[![Codecov](https://codecov.io/gh/complextissue/pyclustree/graph/badge.svg)](https://codecov.io/gh/complextissue/pyclustree)
[![Python Version Required](https://img.shields.io/pypi/pyversions/pyclustree)](https://pypi.org/project/pyclustree/)

Visualize cluster assignments at different resolutions. Possbile applications include finding the optimal resolution for
single-cell RNA-sequencing clusterings.

`pyclustree` is inspired by the R package `clustree`, however, while we aim to provide the same functionality, the API
will differ between the implementations.

## Getting started

Please refer to the [documentation][link-docs]. In particular, the

-   [API documentation][link-api].

## Installation

You need to have Python 3.10 or newer installed on your system. If you don't have
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

> t.b.a

[license]: https://github.com/complextissue/pyclustree/blob/main/LICENSE
[issue-tracker]: https://github.com/complextissue/pyclustree/issues
[changelog]: https://pyclustree.readthedocs.io/latest/changelog.html
[link-docs]: https://pyclustree.readthedocs.io
[link-api]: https://pyclustree.readthedocs.io/latest/api.html
[link-pypi]: https://pypi.org/project/pyclustree
