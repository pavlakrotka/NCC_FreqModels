Statistical modeling to adjust for time trends in adaptive platform
trials utilizing non-concurrent controls
================

This is an accompanying repository for the paper *“Statistical modeling
to adjust for time trends in adaptive platform trials utilizing
non-concurrent controls’’* by Pavla Krotka, Martin Posch and Marta
Bofill Roig. It contains the code to reproduce all simulations and
figures, as well as the case study presented in the paper.
<!-- "[Statistical modeling to adjust for time trends in adaptive platform trials utilizing non-concurrent controls](https://arxiv.org)". -->

The repository is structured as follows:

- Folder **simulations**:

  - *NCC_FreqModels_simscript.R*: This script contains the code to
    reproduce all simulations included in the paper. The results for the
    individual scenarios are then saved in the subfolder *results*.
  - *NCC_FreqModels_figures.Rmd*: This file contains the code to create
    all figures presented in the paper (Section 4) and the supplementary
    material (Section C). The figures are saved in the subfolder
    *figures* in a .pdf and .png formats.

- Folder **case_study**:

  - *NCC_FreqModels_case_study.Rmd*: This file includes the code used
    for illustrating the use of the methods from the paper on the data
    from the ABBV-8E12 trial, as shown in Section 5. It also contains
    the code to reproduce Figure 11, which is saved in the subfolder
    *figures*. For data privacy reasons the dataset used for the
    illustration could not be uploaded.

## NCC R-package

The version of the `NCC` [R-package](https://pavlakrotka.github.io/NCC/)
that was used for the simulation study is included in this repository.

------------------------------------------------------------------------

**Funding**

This research was funded in whole or in part by the Austrian Science
Fund (FWF) \[ESP 442 ESPRIT-Programm\].

P.K. and M.P. received funding from the European Joint Programme on Rare
Diseases (Improve-PSP).
