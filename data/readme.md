## How to compile the data

The dataset `ALL_NEIGHBORHOOD_COVARIATES_YEARLY_PHILLY.csv` has been created by merging together several datasets. 

### ACS 

Most of the data comes from the American Community Survey (ACS) and was downloaded from [data.census.gov](https://data.census.gov/cedsci/). 

To obtain the data aggregated at the correct level (Census Tracts in Philadelphia), use Advanced Search and select 
```
Geography > Tract > Pennsylvania > Philadelphia County, Pennsylvania > All Census Tracts within Philadelphia County, Pennsylvania.
```

The codes for the tables used for each variable and the tranformations used are reported in `File_names_and_calculations.pdf`.

### OpenDataPhilly

Some covariates were downloaded from [OpenDataPhilly.org](https://www.opendataphilly.org), in particular Housing Violations, Violent Crime Rate and Non-Violent Crime Rate.

- Housing violation (L&I Violations) was downloaded from [OpenDataPhilly](https://www.opendataphilly.org/dataset/licenses-and-inspections-violations) in July, 2017 (thus containing data from 2007 to 2017).

- Crime incidents were downloaded from [OpenDataPhilly](https://www.opendataphilly.org/dataset/crime-incidents) in Febrary 2019, aggregated into violent and nonviolent (property) incidents, using the definition of the Uniform Crime Reporting program of the FBI: violent crimes include homicides, rapes, robberies and aggravated assaults, while
nonviolent crimes include burglaries, thefts and motor vehicle thefts. Incidents were also aggregated into Census Tracts by year. Code and instructions to replicate this aggregation are provided [here](https://github.com/cecilia-balocchi/particle-optimization/tree/master/two_partitions/get_data). Violent and nonviolent crime *rates* were computed dividing the number of incidents in each area/year by the corresponding population and multiplying by 100,000.
