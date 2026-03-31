Data and code to reproduce analyses in Lawlor et al., "Global Patterns Predict Local Biodiversity Shifts in a Climate Change Hotspot," https://doi.org/10.32942/X20H2J

Code tasks are annotated and organized by task in the `scripts/` folder. The provided workflow creates all objects in the `data-processed/` folder, and subsequent scripts use these processed data to conduct analyses that are presented in the manuscript. 

One data product referenced within this repository is not provided here, because the object is too largetostore on GitHub: `data-processed/global-temp-data/global_temps_1982_2023.rds`. This object is an array of global sea surface temperatures, referenced within `scripts/03-calculate-STIs/03-match-occurrences-to-global-temp.R` to match OBIS records of species to the temperatures in the time and location of their occurrences. Therefore that script cannot be run from a forked repository without first running `scripts/02-get-sst-data/02-extract-global-sst-1982-2023.R` to recreate the temperature object locally. However, intermediate data objects produced in the temperature matching script and all subsequent scripts are stored in the `data-processed/` folder, so all downstream analyses should be reproducible. 

All outputs, including figures in the manuscript are contained within the `outputs/` folder, some of which were further edited (to add species' silhouettes, etc) outside of R before addition to the final manuscript. 

The script, `master.R` lists all necessary scripts to replicate our analyses, with annotations describing the goals and outputs of each sourced script from the `scripts/` folder. 


