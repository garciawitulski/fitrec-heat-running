# Analysis of individual records

This workflow fits the session, daily and weekly models from prepared FitRec records. It also runs the subgroup analyses, bootstrap and sensitivity checks. Readers can reproduce the main results using the [repository's standard analysis](../README.md) and the aggregate inputs included with it.

To run the individual-record models, you need the licensed panels and the comparison files listed in the [input inventory](input_checksums.csv). The [input guide](INPUTS.md) describes each group of files. Install R packages `data.table`, `DBI`, `duckdb`, `jsonlite` and `sf`. From the repository root, run:

```sh
Rscript --vanilla record_level/run.R
```

The script looks for inputs under `data/record_level/`. Set `FITREC_INPUT_DIR` to use a different directory. If your files have other names or are stored in several places, set `FITREC_INPUT_MAP` to a local CSV as described in the input guide. The command writes tables, figures and numerical checks under `outputs/record_level/`. Set `FITREC_OUTPUT` to write the files to a single directory of your choice.

The bootstrap draws runners using R's random-number generator. Setting `FITREC_PYTHON` to a Python executable with NumPy adds a comparison using the draw sequence used in the paper.
