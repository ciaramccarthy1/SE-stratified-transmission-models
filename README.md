# Respiratory-virus transmissiom model with socio-economic startification of disease states and social contacts

# Project overview
This modelling framework is part of Work Package 4 (WP4) of the Winter Pressures project. 
The project aims to quantify inequality in health impacts and primary care of respiratory infections 
in England during winter using electronic health records accessed through OpenSAFELY.
https://www.opensafely.org/approved-projects/172/. 
The WP4 aims to use transmission modelling approaches to evaluate how demographic factors and unequal vaccination coverage 
affect respiratory virus transmission in winter in different population groups and how these could be mitigated thorugh 
targetted vaccination. In the first stage of WP4, we have extended existing dynamic transmission models for influenza, COVID-19 
and RSV, used previously in the UK for vaccine decision making, to include socio-economic stratification.


# Model framework
The dynamic transmission modela have a single flexible compartmental structure parametersised for each 
of the three respiratory viruses. 

The model inputs (in /data/ ) include:
* literature-based model parameters for COVID-19, Influenza, and SRV.
* Social contact data for England stratified by age and socio-economic groups.
* Demographic data for England from the the ONS 2021 census.

The pipeline is run by running the R script 'main.r', which:
* Sources R and C++ code (in /code/).
* Outputs model epidemic outcomes to /output/ for each respiratory virus.

Output includes (currently set for the urban population):
* Plots of number of infectious people with clinical or sub-clinical infections.
* Plots for the overall population, or by socio-economic strata, or by age group.
* Summary of parameter inputs and settings.
* Diagnostics of performance of the model code.

# Further work
Under development:
* Risk-group definition through age and/or SES to target vaccination programmes.
* Other refinements of vaccination implementation sepecific to each infection.
* For influenza, a single infectious state (for both sub-clinical and clinical infections) is usually used in the model;
* it is under debate whether to split the sub-clinical infections, but if so literature parameter estimates may not be reliable.

Inclusion of hospitalisations:
* A further H compartment can be added to fit secondary care data (if available) or to estimate severe outcomes (if required); 
* robust estimates of clinical outcome parameters exists for COVID-19, but some may not be available for influenza and SRV. 

Fitting primary care disgnostic data or secondary care data, depending of availability:
* The model tracks the number of new clinical infections (all infections, for influenza) scaled by a reporting rate;
* this variable can be used to fit reported case data.

Implementation:
* Tests show the model code is very fast (median run 4-5 milliseconds, depending on platform and running processes) and has high accuracy.
* Further implementation may involve use of C++ numerical libraries; which could involve more compact formulation.
 

| Folder    | Function  |
| :------------ | :----------------------------------------------------------------------------------------------------------------------------------------------- |
|  /            | Main R script for running the pipeline, and diagrams with the structure of the transmission models.
|  /code/       | Code scripts sourced by the pipeline for settting up and running the model and plotting outputs.
|  /data/       | Demographic and social contact data inputs, including figures displaying the data. |
|  /output/     | Files with model outcomes and diagnostics for Influenza, COVID-19, and SRV. |
