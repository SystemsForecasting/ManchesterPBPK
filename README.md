# ManchesterPBPK
Basic physiologically based pharmacokinetics (PBPK) model for performing single subject bottom up predictions.

A running version of the app che be found here: https://manchester.shinyapps.io/pbpk/.

The last version of the software was updated on 03 September 2021.

Authors: *Nicola Melillo & Hitesh Mistry*

## Files
- `server.R`: Shiny server function, it contains the "logics" of the applications (e.g., what happens when clicking a button).
- `ui.R`: Shiny user interface function, it contains the visual aspect of the application.
- `main_pbpk_acat.R`: R script that can be used to run the model outside the Shiny app. This file is not used in `server.R` or `ui.R`.
- `functions/PBPK_model_rxode.R`: it contains the ordinary differential equations of the PBPK and compartmental absorption & transit (CAT) models, written in `RxODE` format.
- `functions/import_param.R`: it contains the functions used to import and elaborate the PBPK and CAT model parameters.
- `functions/functions_plot3.R`: it contains the functions used to generate the plots.
- `data/PBPK_parameters`: it contains `.xlsx` files of the physiological parameters for all the species supported by the Manchester PBPK app. All the references are included in the `.xlsx` files.
- `data/library_drugs`:it contains `.xlsx` files of the drug physicochemical parameters and of drugs pharmacokinetics. All the references are included in the `.xlsx` files and in `data/library_drugs/readme.txt`.

## How to install on a local machine

1. Install R (https://www.r-project.org/).
2. Install RStudio (https://www.rstudio.com/).
3. Open RStudio.
4. Install all libraries: copy, paste and run the following R code in the RStudio console or in a script.
```
install.packages(c("shiny","shinyBS","shinyjs","readxl","writexl","RxODE","dplyr","ggplot2","RColorBrewer","gridExtra","PKNCA","shinybusy"))
```
5. Restart RStudio.
6. Open server.R: double click on the file and it should open automatically.
7. Click on the green button "Run App" to run the app (It  may crash on first use after installing all libraries so please try again).
