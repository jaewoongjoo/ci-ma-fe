# ci-ma-fe

This repository is for the paper **"Valid Causal Inference in Linear Regression From Several Observational Studies None of Which Account for All Confounders"** by Hani Doss and Jaewoong Joo (2025). Here, you'll find all the essential resources, including the manuscript, supplementary materials, and reproducible code.

### Contents

* `ci-ma-fe.pdf`: The complete manuscript of the paper.
* `ci-ma-fe-supp.pdf`: The supplementary material for the paper.
* `ci-ma-fe.R`: The R script for the simulation studies. Running this script generates an `.Rdata` file named `gendata-ci-ma-fe.Rdata`, which is necessary for the analyses performed in `ci-ma-fe-figure1.R` and `ci-ma-fe-supp-figure-s1-s2.R`.
* `ci-ma-fe-figure1.R`: An R script that generates Figure 1 from the paper.
* `ci-ma-fe-figure2.R`: An R script that generates Figure 2 from the paper. This file pertains to the real data analysis. The real data used in our analysis is loaded directly within the code using a URL.
* `ci-ma-fe-supp-figure-s1-s2.R`: An R script that generates Figure S-1 and S-2 from the supplementary material. 
* `ci-ma-fe-supp-figure-s3.R` through `ci-ma-fe-supp-figure-s6.R`: A series of R scripts used to generate Figures S-3 through S-6 found in the supplementary materials. These scripts generate their own `.Rdata` files and do not load any external files.