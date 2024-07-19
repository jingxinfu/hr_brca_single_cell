# The resulting directory structure
```
|-- INSTALL.R <- The required R packages for reproducing the analysis environment
|-- LICENSE
|-- Makefile <- Makefile with commands like `make create_environment`
|-- README.md <- The top-level README for developers using this project.
|-- conda-env.yml <- The requirements file for reproducing the analysis environment
|-- data
│   ├── external       <- Data from third party sources.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump.
|-- notebooks <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
|
|-- report <- cleaned figure and table to share an present
|   |-- figure
|   `-- table
|-- requirements.txt <- The required python packages for reproducing the analysis environment
|-- setup.py <- makes project pip installable (pip install -e .) so src can be imported
`-- src
    |-- Rscripts <- R Scripts to process the data
    |-- __init__.py
    |-- data <- Scripts to download or generate data
    |   `-- __init__.py
    |-- external
    |   `-- __init__.py
    |-- models.py  <- Functions to test hypothesis
    |-- utils.py <- Miscellaneous functions
    `-- visualization.py  <- Functions to create exploratory and results oriented visualizations
```