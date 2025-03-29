# Respiratory virus transmissiom models with socio-economic startification of disease states and social contacts

# Project overview
This modelling framework is part of Work Package 4 (WP4) of the Winter Pressures project. 
The project aims to quantify inequality in the health impacts and in the primary care
of respiratory infections in England during winter, using electronic health records accessed through OpenSAFELY.
https://www.opensafely.org/approved-projects/172/. 
WP4 aims to model winter epidemics of three main respiratory viruses impacting primary and secondary care, 
and to evaluate the potential of vaccination to ameliorate their unequal health impacts across 
socio-economic and ethnically distinct populations.

# Model framework
The dynamic transmission model has a single flexible compartmental structure parametersised for each 
of the respiratory viruses causing COVID-19, Influenza, and Respiratory Syncytial Virus (SRV) illness. 

The current model inputs (in /data/ ) include:
* literature-based model parameters for COVID-19, Influenza, and SRV.
* Social contact data stratified by age group and socio-economic strata
* Demographic data for England from the the ONS 2021 census

The pipeline is run by running the R file main.r, which:
* Sources R and c++ code in /code/
* Outputs model epidemic outcomes to /output/ for each respiratory virus

Output includes (currenlty set for the urban population):
* Plots of number of infectious people with clinical or sub-clinical infections
* Plots versions for overall population, or by socio-economic strata, or by age group
* Summary of parameter inputs and settings
* Diagnostics of performance of the model code

# Further work
Under development:
* Risk group definition through age and/or SES to refine specification of the vaccination programmes.
* Other refinements of vaccination implementation sepecific to each infcetion.

Implementation:
* Tests show the current model code is fast (median 20-25ms, depending on platform and running processes) and has a hight degree of accuracy.
* Further implementation may involve use of c++ numerical libraries; but requires more compact and more error-prone formulation  
 

| Folder    | Function  |
| :------------ | :----------------------------------------------------------------------------------------------------------------------------------------------- |
|  /            | Main R script for running the pipeline and figure with sructure of the trasnmission models.
|  /code/       | Code scripts sourced by the pipeline for settting up and running the model and plotting outputs.
|  /data/       | Demographic and social contact data inputs, including figures displaying the data. |
|  /output/     | Files with model outcomes and diagnostics for Influenza, COVID-19, and SRV. |
