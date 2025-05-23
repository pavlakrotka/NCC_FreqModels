Statistical modeling to adjust for time trends in adaptive platform
trials utilizing non-concurrent controls
================

This is an accompanying repository for the paper *“Statistical modeling
to adjust for time trends in adaptive platform trials utilizing
non-concurrent controls’’* by Pavla Krotka, Martin Posch, Mohamed
Gewily, Günter Höglinger, and Marta Bofill Roig. It contains the code to
reproduce all simulations and figures, as well as the case study
presented in the paper.
<!-- "[Statistical modeling to adjust for time trends in adaptive platform trials utilizing non-concurrent controls](https://arxiv.org)". -->

The repository is structured as follows:

- Folder **simulations**:

  - *NCC_FreqModels_simscript.R*: This script contains the code to
    reproduce all simulations included in the paper. The results for the
    individual scenarios are then saved in the subfolder *results*.
  - *NCC_FreqModels_figures.Rmd*: This file contains the code to create
    all figures presented in the paper (Section 4) and the supplementary
    material (Section C). The figures are saved in the subfolder
    *figures* in a .pdf, .png and .tiff formats.

- Folder **case_studies**:

  - Subfolder **PSP**:

    - *NCC_FreqModels_case_study_PSP.Rmd*: This file includes the code
      used for illustrating the use of the methods from the paper on the
      data from the ABBV-8E12 trial, as shown in Section 5. It also
      contains the code to reproduce Figure 13, which is saved in the
      subfolder *figures*. For data privacy reasons the dataset used for
      the illustration could not be uploaded.
    - *NCC_FreqModels_case_study_PSP_synthetic.Rmd*: This file includes
      the code used for illustrating the use of the methods from the
      paper on synthetic data that mimic the data from the ABBV-8E12
      trial.

  - Subfolder **FLAIR**:

    - *NCC_FreqModels_case_study_PSP.Rmd*: This file includes the code
      used for illustrating the use of the methods from the paper on
      simulated data based on the FLAIR trial, as shown in the
      supplementary material (Section D). It also contains the code to
      reproduce Figure S12, which is saved in the subfolder *figures*.

## NCC R-package

The version of the `NCC` [R-package](https://pavlakrotka.github.io/NCC/)
that was used for the simulation study is included in this repository.

This version is labeled as Release 1.4 on
[GitHub](https://github.com/pavlakrotka/NCC) and can be installed by
running the following code:

``` r
# install.packages("devtools") 
devtools::install_github("pavlakrotka/NCC@v1.4", force = TRUE, build_vignettes = TRUE)
```

Please note that prior to installing the `NCC` package, the
[JAGS](https://mcmc-jags.sourceforge.io/) (\>=3.4.1) library needs to be
installed on your computer, as it is an external dependency of the
package. For more details, see
<https://pavlakrotka.github.io/NCC/articles/installation.html>.

## Working directories

The required working directory for each code file is the folder where
this file is located. E.g., for the files `NCC_FreqModels_simscript.R`
and `NCC_FreqModels_figures.Rmd`, the required working directory is the
folder **simulations**.

## Environment

See below the R version, operating system, and the versions of all R
packages used for running the code for this manuscript (simulations,
figures, as well as case studies).

``` r
> sessionInfo()
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Debian GNU/Linux 12 (bookworm)

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.11.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.11.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Vienna
tzcode source: system (glibc)

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridis_0.6.5     viridisLite_0.4.2 spaMM_4.5.0       lmerTest_3.1-3    lme4_1.1-36       Matrix_1.7-0      scales_1.3.0     
 [8] ggpubr_0.6.0      latex2exp_0.9.6   kableExtra_1.4.0  lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4      
[15] purrr_1.0.2       readr_2.1.5       tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.1     tidyverse_2.0.0   NCC_1.0          

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1    rjags_4-16          loo_2.8.0           fastmap_1.2.0       digest_0.6.37       BayesPPD_1.1.3      timechange_0.3.0   
 [8] lifecycle_1.0.4     StanHeaders_2.32.10 magrittr_2.0.3      compiler_4.4.2      rlang_1.1.5         tools_4.4.2         ggsignif_0.6.4     
[15] knitr_1.46          pkgbuild_1.4.6      curl_5.2.1          xml2_1.3.6          abind_1.4-8         registry_0.5-1      withr_3.0.2        
[22] numDeriv_2016.8-1.1 grid_4.4.2          stats4_4.4.2        colorspace_2.1-1    future_1.34.0       inline_0.3.21       globals_0.16.3     
[29] iterators_1.0.14    MASS_7.3-61         cli_3.6.4           mvtnorm_1.3-3       rmarkdown_2.26      crayon_1.5.3        reformulas_0.4.0   
[36] generics_0.1.3      RcppParallel_5.1.10 rstudioapi_0.16.0   future.apply_1.11.3 tzdb_0.4.0          minqa_1.2.8         pbapply_1.7-2      
[43] proxy_0.4-27        rstan_2.32.7        assertthat_0.2.1    parallel_4.4.2      matrixStats_1.5.0   vctrs_0.6.5         boot_1.3-31        
[50] carData_3.0-5       jsonlite_1.8.8      slam_0.1-55         car_3.1-2           hms_1.1.3           rstatix_0.7.2       Formula_1.2-5      
[57] listenv_0.9.1       systemfonts_1.0.6   foreach_1.5.2       glue_1.8.0          ROI_1.0-1           parallelly_1.42.0   nloptr_2.1.1       
[64] codetools_0.2-20    stringi_1.8.4       gtable_0.3.6        QuickJSR_1.6.0      munsell_0.5.1       doFuture_1.0.1      pillar_1.10.1      
[71] htmltools_0.5.8.1   R6_2.6.1            Rdpack_2.6.2        evaluate_0.24.0     lattice_0.22-6      rbibutils_2.3       backports_1.5.0    
[78] RBesT_1.8-1         broom_1.0.5         rstantools_2.4.0    Rcpp_1.0.14         svglite_2.1.3       coda_0.19-4.1       gridExtra_2.3      
[85] nlme_3.1-166        checkmate_2.3.2     xfun_0.46           mgcv_1.9-1          pkgconfig_2.0.3 
```

------------------------------------------------------------------------

**Funding**

This research was funded in whole or in part by the Austrian Science
Fund (FWF) \[ESP 442 ESPRIT-Programm\].

P.K., M.P., M.G., and G.H. received funding from the European Joint
Programme on Rare Diseases (Improve-PSP).
