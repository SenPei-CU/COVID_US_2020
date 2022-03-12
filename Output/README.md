# COVID_US_2020

Outputs for Sen Pei, Teresa K. Yamana, Sasikiran Kandula, Marta Galanti, Jeffrey Shaman, "Overall burden and characteristics of COVID-19 in the United States during 2020", 2021.

1. ActiveContagiousPrevalence.csv: fraction of population that are contagious in each county on each day.
2. AscertainmentRate.csv: monthly ascertainment rate in each county from March to December 2020. Ascertainment rates for counties with monthly total reported cases <100 are not shown, marked with “NaN”.
3. CFR.csv: case fatality rate in each county on each day. For day t, if the average daily death per 100,000 people from day t-14 to day t+14 is below 0.5, CFR is not shown and marked with “NaN”.
4. CountyPopulation.csv: Population in each county.
5. Death_deconvoluted.csv: deconvoluted daily death for cases reported on each day in each county.
6. IFR.csv: infection fatality rate in each county on each day. For day t, if the average daily death per 100,000 people from day t-14 to day t+14 is below 0.5, CFR is not shown and marked with “NaN”.
7. Susceptibility.csv: daily susceptibility in each county on each day. A few counties with small population and abnormal reporting issues were estimated to have susceptibility below 25%. We ignore those potentially biased estimates and mark those counties with "NaN"
8. EstimatedDailyInfection.csv: estimated daily new infections (both reported and unreported) in each county on each day. Data are reported as 2.5, 25, 50, 75, and 97.5 percentiles.
