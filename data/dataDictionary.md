We analyzed two datasets in our article. Both are subsets of total precipitation rate (in m/s) files, produced by the Community Earth System Model (CESM) Large Ensemble Project: https://www.cesm.ucar.edu/community-projects/lens

Specifically, we considered log-transformed total precipitation rate in July in 98 consecutive years starting in the year 402. 
For ease of working with the data and reproducing our results, we provide the two datasets in this GitHub repository.

prec_days.RData contains three objects:

- lat: vector of length 2738 containing latitude values
- lon: vector of length 2738 containing longitude values
- precs: 3D array of size 2738x98x30 with dimensions corresponding to spatial locations, years, and days, respectively.

prec_all.RData contains three objects:

- lat: vector of length 55296 containing latitude values
- lon: vector of length 55296 containing longitude values
- prec: 2D array of size 55296x98 with dimensions corresponding to spatial locations and years, respectively.
