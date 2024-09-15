# Installation

:::{card}
:class-card: sd-bg-warning
:class-body: sd-bg-text-warning
**pyclustree** only supports Python versions greater than or equal to **3.9**.
:::

## Installation Options

Choose an option to install this package.

::::{tab-set}

:::{tab-item} PyPi
Install `pyclustree` package using `pip`:

```bash
pip install pyclustree
```

:::

:::{tab-item} Source
Install `pyclustree` from source:

```bash
# Clone repo
git clone --depth 1 https://github.com/complextissue/pyclustree.git
cd pyclustree
pyenv global 3.9
make create-venv
source .venv/bin/activate
make install
```

:::

::::
