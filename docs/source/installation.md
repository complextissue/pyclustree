# Installation

:::{card}
:class-card: sd-bg-warning
:class-body: sd-bg-text-warning
**pyclustree** only supports Python versions greater than or equal to **3.10**.
:::

## Installation Options

Choose an option to install this package.

::::{tab-set}

:::{tab-item} PyPi
Install `pyclustree` package using `pip`:

```bash
pip install pyclustree
```

For scoring of clusters functions from the `scikit-learn` package are used. To install `pyclustree` with this optional
dependency use the following `pip` command:

```bash
pip install "pyclustree[sklearn]"
```

:::

:::{tab-item} Source
Install `pyclustree` from source:

```bash
# Clone repo
git clone --depth 1 https://github.com/complextissue/pyclustree.git
cd pyclustree
pyenv global 3.10
make create-venv
source .venv/bin/activate
make install
```

:::

::::
