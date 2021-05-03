# COVID_US_2020

Code for Sen Pei, Teresa K. Yamana, Sasikiran Kandula, Marta Galanti, Jeffrey Shaman, "Overall burden and characteristics of COVID-19 in the United States during 2020", 2021.

Code to run inference for COVID-19 spread in the US at county level.

## Functions
1. infer.m: the main function to run inference
2. initialize.m: compute initial conditions for all subpopulations
3. seeding.m: set initial infection for all subpopulations
4. initializepara_eakf.m: initialize the ensemble of parameters for data assimilation
5. model_eakf.cpp: run the transmission model

## Data and data structure
1. commutedata.mat:
a. The commuting network structure is stored in two vectors: nl (neighbor list) and part (partition). The vector nl records all neighbors of counties connected in the commuting network, and the vector part specifies the index range for the neighbors of each county. To find the list of counties where commuters living in county i work, use nl(part(i):part(i+1)-1). For computational convenience, the first county in the list is county i itself.
b. The number of commuters is stored in the vector C. For instance, the number of commuters from county i to county j=nl(k) [where part(i)<=k<=part(i+1)-1] is C(k). Note the index k should be within the index range for county i, i.e., from part(i) to part(i+1)-1.
c. The average number of commuters between two counties is stored in the vector Cave. For instance, the average number of commuters between county i and county j=nl(k) [where part(i)<=k<=part(i+1)-1] is Cave(k). It is the average number of commuters in both directions.
2. countyfips.mat: the FIPS code of US counties
3. dailyincidence.mat: the number of confirmed cases in each county on each day. Row: county; column: day, starting from Feb 21 2021.
4. delaypara.mat: the parameters for the gamma distributions of reporting delay in each month from April to December 2020. Row: month from April; column: the shape parameter (first column) and scaling parameter (second column)
5. MI_inter.mat: county-level SafeGraph mobility data in 2020. The first column is county FIPS code. Other columns are the ratios of daily numbers of visitors from outside county to the number of inter-county visitors on March 1 2020.
6. parafit.mat: model parameters estimated using data until March 13. Row: beta, mu, Z, D, alpha, theta; column: 100 ensemble members for each parameter
7. popd.mat: population density of each county, persons per square miles
8. population.mat: county population

## How to run inference
To speed up the inference, the transmission model is programmed in C++. The inference function calls the C++ function through interface between MATLAB and C++.

1. Compile the C++ function in MATLAB using “mex model_eakf.cpp”. A C++ compiler (e.g., Xcode) needs to be installed on the computer before running this command.
2. Run infer.m

## Outputs
1. para_post: posterior parameters, each with 100 ensemble members
2. S_post: posterior susceptibility in each county on each day, 100 ensemble members
3. dailyIr_post_rec: posterior estimate of the daily newly infected documented infections in each county on each day, 100 ensemble members
4. dailyIu_post_rec: posterior estimate of the daily newly infected undocumented infections in each county on each day, 100 ensemble members
5. obs_temp: the posterior fitting of daily confirmed cases in each county, 100 ensemble members

