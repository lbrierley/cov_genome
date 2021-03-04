# Predicting hosts of coronaviruses from genome composition biases

Supporting data and code associated with Brierley, L., & Fowler, A. (2020). Predicting the animal hosts of coronaviruses from compositional biases of spike protein and whole genome sequences through machine learning. bioRxiv preprint, 2020.11.02.350439. https://doi.org/10.1101/2020.11.02.350439

## Scripts

`startup.R` initialises options and runs all analytical scripts. Note that an NCBI API key can be filled in here in order to use `process_cov_seq.R` to extract new data at a faster rate. See https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#rate-limiting-and-api-keys for more information.

`functions.R` defines the custom functions required by `process_cov_seq.R` and `figures_tables.R`

`process_cov_seq.R` extracts and processes sequences and their associated metadata from NCBI's GenBank (or loads previously extracted results), before identifying corresponding proteins for each coding sequence and calculating genome composition bias features (nucleotide, dinucleotide, and codon biases).

`build_random_forests.R` trains random forest models to predict host category of coronavirus sequences from genomic bias features using an outer loop of hold-one-species-out cross-validation and an inner loop of 10-fold cross-validation. Models are trained separately for spike protein sequences and whole genome sequences.

`partial_dependence.R` calculates partial dependence profiles for the four most important features that can then be plotted to interpret and visualise how predicted probability of host category changes as a function of specific genomic bias measures.

`build_random_forests_samplers.R` acts as `build_random_forests.R` but tests five additional data resampling methods (downsampling, upsampling, SMOTE, class weighting and no resampling).

These three scripts are not run by default as they require significant computing time. For `build_random_forest.R` and `partial_dependence.R` saved output .RData files are provided in `/outputs`.

`figures_tables.R` reconstructs figures and tables associated with Brierley & Fowler 2020.
