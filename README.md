## Repository to hold data and code for ICESJMS cod allee effects

Data and R code for the analysis of depletion, depensation and recovery in 24 Atlantic cod stocks (ICES Journal of Marine Science, manuscript ICESJMS-2026-338).

### Contents

- `data/cod_clean.RDS`: annual catch and biomass by stock, compiled from DFO stock assessments and the RAM Legacy Stock Assessment Database  
- `R/analysis.R`: surplus production, segmented and linear regressions, Allee effect thresholds and strength, depletion and recovery status, main figures  
- `R/supp.R`: GAM production curves and Bmsy proxy  

### How to run

1. Open the R project file at the repository root (paths use the `here` package).
2. Run `R/analysis.R`.
3. Run `R/supp.R` in the same session; it uses objects created by `analysis.R`.

Missing packages are installed automatically.

### Citation

[Authors] ([Year]). [Title]. ICES Journal of Marine Science. [DOI]

### Contact

François Turcotte, Fisheries and Oceans Canada (francois.turcotte@dfo-mpo.gc.ca)
