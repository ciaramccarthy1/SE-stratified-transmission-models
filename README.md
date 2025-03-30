# Respiratory-virus transmissiom models with socio-economic startification of disease states and social contacts

# Project overview
This modelling framework is part of Work Package 4 (WP4) of the Winter Pressures project. 
The project aims to quantify inequality in health impacts and primary care of respiratory infections 
in England during winter using electronic health records accessed through OpenSAFELY.
https://www.opensafely.org/approved-projects/172/. 
WP4 aims to use transmission modelling approaches to evaluate how demographic factors and unequal vaccination coverage 
affect respiratory virus transmission in winter in different population groups and how these could be mitigated thorugh 
targetted vaccination. As the first stage, we have extended existing dynamic transmission models for influenza, COVID-19 and RSV, 
used previously in the UK for vaccine decision making, to include stratification by socioeconomic quintile.


# Model framework
The dynamic transmission model has a single flexible compartmental structure parametersised for each 
of the three respiratory viruses considered. 

The current model inputs (in /data/ ) include:
* literature-based model parameters for COVID-19, Influenza, and SRV.
* Social contact data stratified by age group and socio-economic strata
* Demographic data for England from the the ONS 2021 census

The pipeline is run by running the R file main.r, which:
* Sources R and c++ code in /code/
* Outputs model epidemic outcomes to /output/ for each respiratory virus

Output includes (currently set for the urban population):
* Plots of number of infectious people with clinical or sub-clinical infections
* Plots versions for overall population, or by socio-economic strata, or by age group
* Summary of parameter inputs and settings
* Diagnostics of performance of the model code

# Further work
Under development:
* Risk-group definition through age and/or SES to refine the specification of vaccination programmes.
* Other refinements of vaccination implementation sepecific to each infection.
* For influenza, a single infectious state is usually used in the model
* it's under debate whether to separate sub-clinical infections, but if so lit parameter estimates may not be reliable.

Inclusion of hospitalisations:
* A further H compartment can be added to fit secondary care data (if available) or estimate severe outcomes (if required); 
* robust estimates of clinical parameters exists for COVID-19, but some may not be available for influenza and SRV. 

Fitting primary care disgnostic data or secondary care data, depending of availability:
* The model tracks the number of new clinical infections (all infections, for influenza) scaled by a reporting rate;
* this variable can be used to fit reported case data 

Implementation:
* Tests show the current model code is very fast (median run 4-5ms, depending on platform and running processes) and has high accuracy.
* Further implementation may involve use of c++ numerical libraries; which requires more compact and error-prone formulation  
 

| Folder    | Function  |
| :------------ | :----------------------------------------------------------------------------------------------------------------------------------------------- |
|  /            | Main R script for running the pipeline and figure with sructure of the trasnmission models.
|  /code/       | Code scripts sourced by the pipeline for settting up and running the model and plotting outputs.
|  /data/       | Demographic and social contact data inputs, including figures displaying the data. |
|  /output/     | Files with model outcomes and diagnostics for Influenza, COVID-19, and SRV. |
