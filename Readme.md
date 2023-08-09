# PALAGRAM Monograph

This repository contains the data and codes associated with the PALAGRAM consortium.

## Repository organization

```
  palagram_monograph
  │
  └───data: data are stored here
  │   └───input_data: input data as sent by everyone
  │       └─── ... : NETCDF files
  │   └───output_data: processed data output by analysis.py (also contains input_data)
  │       └─── ... : NETCDF files
  └───paper: contains source files for article
  │   └───figures: contains source figures
  │       └─── ... : PDF files
  │       └─── figure_scripts: contains figure scripts that reads data in data/output_data and writes figures in paper/figures
  │            └─── *.py : python scripts for figures
  │   └─── ... : various files (.tex, .bib, ...)
  │   └─── main.pdf : article preprint
  └───analysis:
      └───analysis.py: analysis code, that reads input_data and writes output_data

  ```

## Workflow

One can work separately on the analysis (`analysis.py`), or the figures. While changing and running analysis.py, `output_data` are automatically updated. Figure scripts then needs to be re-run to update them (Use `figure_scripts/run_all.py` to update them all).

### Working on the repo

If you have git installed, clone the repository: `git clone git@github.com:Cgadal/palagram_monograph.git`.

> **Warning**: When working locally, please create and work on a separate branch, different from `main`. 

Then:

- commit your work
- pull to update the `main` branch
- merge the `main` branch with your work branch and solve potential conflicts
- only then, push on the distant repo (GitHub)

If you don't have git installed, you won't be able to work on the repo. You can still download the whole content by clicking on the green <> code button, or simply clicking [here](https://github.com/Cgadal/palagram_monograph/archive/refs/heads/main.zip).