# Analysis of individual records

Use these scripts if you have the prepared FitRec records and want to refit the session, daily or weekly models, run the subgroup and sensitivity analyses, or resample runners. The [main analysis](../README.md) reproduces the main contrasts from the aggregate data included in the repository.

To run the individual-record models, you need the licensed panels and the comparison files listed in the [input inventory](input_checksums.csv). The [input guide](INPUTS.md) describes each group of files. Install R packages `data.table`, `DBI`, `duckdb`, `jsonlite` and `sf`. From the repository root, run:

```sh
Rscript --vanilla record_level/run.R
```

The script looks for inputs under `data/record_level/`. Set `FITREC_INPUT_DIR` to use a different directory. If your files have other names or are stored in several places, set `FITREC_INPUT_MAP` to a local CSV as described in the input guide. The command writes tables, figures and numerical checks under `outputs/record_level/`. Set `FITREC_OUTPUT` to write the files to a single directory of your choice.

By default, R draws the runners for each bootstrap sample. To reproduce the exact draws used for the paper, set `FITREC_PYTHON` to the path of a Python executable with NumPy installed. The models are fitted in R in both cases.
