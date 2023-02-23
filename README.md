# Educational PBPK, PBPKedu
This physiologically based pharmacokinetic (PBPK) model is designed as an educational tool, hence named PBPKedu. The aim is to provide a simple, open source, freely downloadable PBPK model for performing basic single subject bottom-up simulations in R. The model has not been tested against a large set of compounds. All parameters were taken from the literature. The PBPK model was implemented by using the [rxode2](https://cran.r-project.org/web/packages/rxode2/index.html) package and non-compartmental analysis (NCA) was performed with the [PKNCA](https://cran.r-project.org/web/packages/PKNCA/index.html) package.

A running version of the app (might be not updated) can be found here: https://manchester.shinyapps.io/pbpk/.

Additions, suggestions and amendments are welcome.

## Files
- `server.R`: Shiny server function, it contains the "logics" of the applications (e.g., what happens when clicking a button).
- `ui.R`: Shiny user interface function, it contains the visual aspect of the application.
- `main_pbpk_acat.R`: R script that can be used to run the model outside the Shiny app. This file is not used in `server.R` or `ui.R`.
- `functions/PBPK_model_rxode.R`: it contains the ordinary differential equations of the PBPK and compartmental absorption & transit (CAT) models, written in `RxODE` format.
- `functions/import_param.R`: it contains the functions used to import and elaborate the PBPK and CAT model parameters.
- `functions/functions_plot4.R`: it contains the functions used to generate the plots.
- `data/PBPK_parameters`: it contains `.xlsx` files of the physiological parameters for all the species supported by the PBPKedu app. All the references are included in the `.xlsx` files.
- `data/library_drugs`:it contains `.xlsx` files of the drug physicochemical parameters and of drugs pharmacokinetics. All the references are included in the `.xlsx` files and in `data/library_drugs/readme.txt`.

## How to install on a local machine

1. Install R (select your favourite CRAN mirror from https://cran.r-project.org/mirrors.html and dowload the latest version of R).
2. Install RStudio (https://posit.co/download/rstudio-desktop/).
3. If you have Windows OS, [install RTools](https://cran.r-project.org/bin/windows/Rtools/) (check your R version!). If you have Ubuntu, you don't need RTools, but probably you will need to install some libraries (i.e., [build-essential](https://github.com/SystemsForecasting/nlmixr2_on_AWS#install-r-and-rstudio), [make, cmake etc.](https://github.com/SystemsForecasting/nlmixr2_on_AWS#install-tidyverse-and-nlmixr2)). We haven't tested the installation on MacOS (if you do, please provide us some feedback!).
4. Open RStudio.
5. Install all libraries: copy, paste and run the following R code in the RStudio console or in a script (it may take a while...).
```
install.packages(c("shiny","shinyBS","shinyjs","readxl","writexl","rxode2","dplyr","ggplot2","RColorBrewer","gridExtra","PKNCA","shinybusy"))
```
6. Restart RStudio.
7. Download and unzip the ManchesterPBPK code. To download the code, press the green `<> code` button above and press `Download ZIP`.
8. Open the unzipped folder containing all the code. You should see the same files and folders present at the top of this page. Open server.R with RStudio: double click on the file and it should open automatically. If it does not work, right-click on the file and select open with RStudio.
9. Now, you should have the server.R code opened in RStudio. Next to the usual "Run" button, there should be a green button "Run App": click it to run the app (it  may crash on first use after installing all libraries so please close RStudio and try again point 9).


## Standalone version
An installer for a standalone version is made avaiable [here](https://drive.google.com/file/d/1tJAGFH0A9wbhUvd1z6KdHj7yJZJeBbGm/view?usp=sharing) (through Google Drive). Windows 10 and 11 are the only supported OS.  
With this standalone version you can play with the app also without installing R, RStudio and so on. To install the PBPKedu app, follow these steps.  

1. [Install RTools](https://cran.r-project.org/bin/windows/Rtools/) for the R 4.2 version. Unfortunately, this is still needed!
2. Dowload `installer_PBPKedu_windows10.exe` (despite the name it should work for Windows 11 as well) at [this link](https://drive.google.com/file/d/1tJAGFH0A9wbhUvd1z6KdHj7yJZJeBbGm/view?usp=sharing).
3. Open (double click) `installer_PBPKedu_windows10.exe`.
4. Allow the app to make changes on the device.
5. Select the app location, default should be `C:\Program Files (x86)\PBPKedu`.
6. Tick create desktop shortcut.
7. Continue and press install.

In a few minutes the app should be installed and a '90s style `PBPKedu` icon should appear on your desktop. Open the app: you should see PBPKedu opened in your default web browser.  
Enjoy!

