# CDC_transcriptomics_ms
This repository contains R scripts to reproduce the main and supplementary figures presented in [Gargiuli et al.](https://www.mdpi.com/2072-6694/13/12/2903).

## Instructions
To run the scripts create the CDC folder in your user home directory. Download the repo and save the scripts in the newly created folder.
For example, for Unix-based systems:

```
$ mkdir -p CDC
```

Open R/RStudio to run the first script 0-data_retrieval.R. This script generates all necessary datasets for downstream analyses.
You can simply tipe in R/RStudio :

```
source("~/CDC/scripts/0-data_retrieval.R")
```

After completion of the execution you can run other scripts in the same manner (just change the name of the script in the command line above).
