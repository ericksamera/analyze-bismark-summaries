## analyze-bismark-summaries
**A Python-wrapped toolkit for summarizing methylation data from `process-bismark-coverage.py`.**

## Table of contents
* [Background](#background)
* [Requirements](#requirements)
* [Usage](#usage)
    * [Per-CpG boxplots](#per-cpg-boxplots)
    * [PCoA across all CpG](#pcoa-across-all-cpg)

## Background
The script, `process-bismark-coverage.py`, produces technically human-readable 
`.csv` files. This set of scripts generates visualizations of some relevant
methylation statistics.

## Requirements
* `python`>=3.11
* `plotly`>=5.14.0
* `pandas`>=1.5.3
* `scikit-bio`>0.5.8

Alternatively use the provided `environment.yml` file to create a conda
environment containing all the prequisites.

```
conda env create --name methylation-analysis --file=environment.yml
```

## Usage
TODO: add some actual usage descriptions.

1. If you'd like the plug-and-play version, clone this repo and copy over the 
`individual_primer/` and `summaries/` to this repo's root directory.

2. Then activate your preferred conda environment which fulfills the requirements.

    ```
    conda activate methylation-analysis
    ```

3. Simply run the plug-and-play bash script.
    ```
    bash src/i_pipeline.sh
    ```
4. It should produce 2 `.html` summary documents in the current working 
    directory.

### Per-CpG boxplots
TODO: add some pictures

### PCoA across all CpG
TODO: add some pictures
